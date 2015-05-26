#!/usr/bin/env python

import argparse
import sys
import vcf

from argparse import FileType
from biopericles.SNPFinder import SNPSitesReader
from gt import GFF3InStream, FeatureIndexMemory, FeatureStream

if __name__ == '__main__':
  parser=argparse.ArgumentParser()
  parser.add_argument('vcf_in', type=FileType('r'))
  parser.add_argument('gff_in', type=FileType('r'))
  parser.add_argument('vcf_out', type=str)

  args = parser.parse_args()

  vcf_input_file = args.vcf_in
  gff_input_file = args.gff_in
  if args.vcf_out == '-':
    vcf_output_file = sys.STDOUT
  else:
    vcf_output_file = open(args.vcf_out, 'w')

  vcf_input_reader = SNPSitesReader(vcf_input_file)
  vcf_output_writer = vcf.Writer(vcf_output_file, vcf_input_reader)

  gff_filename = args.gff_in.name
  gff_stream = GFF3InStream(gff_filename)
  feature_index = FeatureIndexMemory()
  feature_stream = FeatureStream(gff_stream, feature_index)

  gn = feature_stream.next_tree()
  while gn:
    gn = feature_stream.next_tree()

  chromosome_name_in_vcf = '1'
  chromosome_name_in_gff = 'Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1|SC|contig000001'

  for record in vcf_input_reader:
    get_type = lambda feature: feature.get_type()
    if record.CHROM == chromosome_name_in_vcf:
      features = feature_index.get_features_for_range(record.POS, record.POS, chromosome_name_in_gff)
      feature_types = [feature.get_type() for feature in features]
      if 'CDS' in feature_types:
        vcf_output_writer.write_record(record)
    vcf_output_writer.write_record(record)
