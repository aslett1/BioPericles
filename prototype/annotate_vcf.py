#!/usr/bin/env python

import argparse
import Bio.SeqIO
import logging
import sys
import vcf

from argparse import FileType
from biopericles.SNPFinder import SNPSitesReader
from Bio.Seq import Seq
from gt import GFF3InStream, FeatureIndexMemory, FeatureStream

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_arguments():
  parser=argparse.ArgumentParser()
  parser.add_argument('vcf_in', type=FileType('r'))
  parser.add_argument('gff_in', type=FileType('r'))
  parser.add_argument('fasta_in', type=FileType('r'))
  parser.add_argument('vcf_out', type=str)
  return parser.parse_args()

def get_file_handles(args):
  vcf_input_file = args.vcf_in
  gff_input_file = args.gff_in
  fasta_input_file = args.fasta_in
  if args.vcf_out == '-':
    vcf_output_file = sys.STDOUT
  else:
    vcf_output_file = open(args.vcf_out, 'w')
  return vcf_input_file, gff_input_file, fasta_input_file, vcf_output_file

def add_GT_format_header(vcf_input_reader):
  GT_header = vcf.parser._Format(id='GT',
                                 num=1,
                                 type='String',
                                 desc='Genotype')
  vcf_input_reader.formats['GT'] = GT_header

def remove_AB_info_header(vcf_input_reader):
  try:
    del(vcf_input_reader.infos['AB'])
  except KeyError:
    pass # Strange if AB is missing but not a problem

def add_consequences_info_header(vcf_input_reader):
  desc='Consequence type in a VEP like format. Format: Allele||||Consequence|||Protein_position|Amino_acids||||STRAND'
  vcf_input_reader.infos['CSQ_like'] = vcf.parser._Info(id='CSQ_like', num='.',
                                                        type='String',
                                                        desc=desc,
                                                        source=None,
                                                        version=None)

def build_feature_index(gff_input_file):
  gff_filename = gff_input_file.name
  gff_stream = GFF3InStream(gff_filename)
  feature_index = FeatureIndexMemory()
  feature_stream = FeatureStream(gff_stream, feature_index)

  gn = feature_stream.next_tree()
  while gn:
    gn = feature_stream.next_tree()

  return feature_index

def get_matching_CDS(record, feature_index,
                     chromosome_name_in_gff):
  features = feature_index.get_features_for_range(record.POS, record.POS, chromosome_name_in_gff)
  return [feature for feature in features if feature.get_type() == 'CDS']

def add_consequences_to_record(record, consequences):
  consequences_string = ','.join((consequence.format() for consequence in consequences))
  record.add_info('CSQ_like', consequences_string)

def get_strand_from_feature(feature):
  strand = feature.get_strand()
  if strand == '+':
    return '1'
  elif strand == '-':
    return '-1'
  else:
    raise ValueError("Expected strand for %s to be '+' or '-'; got %s" %
                     (feature, strand))

class Consequence(object):
  def __init__(self,  Allele, Consequence=None, Protein_position=None, Amino_acids=None, STRAND=None):
    self.Allele=Allele
    self.Consequence=Consequence
    self.Protein_position=Protein_position
    self.Amino_acids=Amino_acids
    self.STRAND=STRAND

  def format(self):
    none_to_empty_str = lambda s: '' if s == None else s
    values = map(none_to_empty_str, [self.Allele, self.Consequence,
                                     self.Protein_position, self.Amino_acids,
                                     self.STRAND])
    values = map(str, values)
    return "%s||||%s|||%s|%s||||%s" % tuple(values)

  __repr__ = format

def alter_base(sequence_obj, index, alternative_base):
  sequence_string = str(sequence_obj)
  sequence_list = [base for base in sequence_string]
  sequence_list[index] = alternative_base
  alternative_sequence_string = ''.join(sequence_list)
  return Seq(alternative_sequence_string, sequence_obj.alphabet)

def analyse_mutations(sequence_obj, vcf_record, cds_feature):
  # Makes some pretty naive assumptions about deletions in the sample;
  # probably wrong.

  # Local naming convention for sanity:
  # position_ => '1' indexed as in original vcf, gff files
  # index_ => '0' indexed as used in python objects

  snp_position = vcf_record.POS
  index_of_snp_in_protein = snp_position - cds_feature.start 
  index_of_amino_acid_in_protein = int(index_of_snp_in_protein/3)*3
  position_in_protein = int(index_of_snp_in_protein/3) + 1 # in amino acids, not bases
  position_of_amino_acid = index_of_amino_acid_in_protein + cds_feature.start
  amino_acid_sequence = sequence_obj.seq[position_of_amino_acid-1:position_of_amino_acid+2]
  original_amino_acid = amino_acid_sequence.translate()
  index_of_snp_in_amino_acid = snp_position - position_of_amino_acid 
  for alternative_base in (str(alt) for alt in vcf_record.ALT):
    alt_amino_acid_sequence = alter_base(amino_acid_sequence,
                                         index_of_snp_in_amino_acid,
                                         alternative_base)
    alt_amino_acid = alt_amino_acid_sequence.translate()
    yield (alternative_base, original_amino_acid, alt_amino_acid,
           position_in_protein)

def get_consequences(sequence_obj, vcf_record, cds_features):
  if len(cds_features) == 0:
    return [Consequence(str(alternative_base), 'intergenic_variant') for alternative_base in
            vcf_record.ALT]
  elif len(cds_features) > 1:
    logger.warn("Matched more than one CDS feature at POS '%s' on CHROM '%s'; skipping" %
                (record.POS, record.CHROM))
    return []
  else:
    consequences = []
    (cds_feature,) = cds_features # takes the first one, errors if there isn't just 1
    mutations = analyse_mutations(sequence_obj, vcf_record, cds_feature)
    strand = get_strand_from_feature(cds_feature)
    for alternative_base, original_amino_acid, new_amino_acid, position_in_protein in mutations:
      if str(original_amino_acid) == str(new_amino_acid):
        consequence_type = 'synonymous_variant'
        amino_acid_change = str(original_amino_acid)
      else:
        consequence_type = 'nonsynonymous_variant'
        amino_acid_change = "%s/%s" % (str(original_amino_acid),
                                       str(new_amino_acid))
      consequence = Consequence(alternative_base,
                                Consequence=consequence_type,
                                Protein_position=position_in_protein,
                                Amino_acids=amino_acid_change,
                                STRAND=strand)
      consequences.append(consequence)
  return consequences

if __name__ == '__main__':
  args = get_arguments()
  vcf_input_file, gff_input_file, fasta_input_file, vcf_output_file = get_file_handles(args)

  vcf_input_reader = SNPSitesReader(vcf_input_file)
  add_consequences_info_header(vcf_input_reader)
  add_GT_format_header(vcf_input_reader)
  remove_AB_info_header(vcf_input_reader)
  vcf_output_writer = vcf.Writer(vcf_output_file, vcf_input_reader)
  sequence = Bio.SeqIO.parse(fasta_input_file, 'fasta').next()

  feature_index = build_feature_index(gff_input_file)

  chromosome_name_in_vcf = '1'
  chromosome_name_in_gff = 'Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1|SC|contig000001'

  for record in vcf_input_reader:
    if record.CHROM == chromosome_name_in_vcf:
      matching_cds = get_matching_CDS(record, feature_index,
                                      chromosome_name_in_gff)
      consequences = get_consequences(sequence, record, matching_cds)
      add_consequences_to_record(record, consequences)
    vcf_output_writer.write_record(record)
