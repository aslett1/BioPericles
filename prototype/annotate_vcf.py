#!/usr/bin/env python

import argparse
import sys
import vcf

from argparse import FileType
from biopericles.SNPFinder import SNPSitesReader
from gt import GFF3InStream, FeatureIndexMemory, FeatureStream

def get_arguments():
  parser=argparse.ArgumentParser()
  parser.add_argument('vcf_in', type=FileType('r'))
  parser.add_argument('gff_in', type=FileType('r'))
  parser.add_argument('vcf_out', type=str)
  return parser.parse_args()

def get_file_handles(args):
  vcf_input_file = args.vcf_in
  gff_input_file = args.gff_in
  if args.vcf_out == '-':
    vcf_output_file = sys.STDOUT
  else:
    vcf_output_file = open(args.vcf_out, 'w')
  return vcf_input_file, gff_input_file, vcf_output_file

def add_consequences_info_header(vcf_input_reader):
  desc='Consequence type in a VEP like format. Format: Allele|Consequence|Protein_position|Amino_acids|STRAND">'
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
    return '|'.join(values)

  __repr__ = format

if __name__ == '__main__':
  args = get_arguments()
  vcf_input_file, gff_input_file, vcf_output_file = get_file_handles(args)

  vcf_input_reader = SNPSitesReader(vcf_input_file)
  add_consequences_info_header(vcf_input_reader)
  vcf_output_writer = vcf.Writer(vcf_output_file, vcf_input_reader)

  feature_index = build_feature_index(gff_input_file)

  chromosome_name_in_vcf = '1'
  chromosome_name_in_gff = 'Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1|SC|contig000001'

  for record in vcf_input_reader:
    if record.CHROM == chromosome_name_in_vcf:
      matching_cds = get_matching_CDS(record, feature_index,
                                      chromosome_name_in_gff)
      print record, len(matching_cds)
      if len(matching_cds) > 0:
        consequences = [Consequence(str(allele), 'intra_genic') for allele in record.ALT]
      else:
        consequences = [Consequence(str(allele), 'inter_genic') for allele in record.ALT]
      add_consequences_to_record(record, consequences)
    vcf_output_writer.write_record(record)
