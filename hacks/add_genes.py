#!/usr/bin/env python

import sys

filename = sys.argv[1]
try:
  temp_output = sys.argv[2]
except IndexError:
  temp_output = filename + ".tmp"

with open(filename, 'r') as gff_in_file:
  with open(temp_output, 'w') as gff_out_file:
    for i,line in enumerate(gff_in_file):
      try:
        if line[:2] == '##':
          gff_out_file.write(line)
        else:
          elements = line.split('\t')
          feature_type = elements[2]
          if feature_type == 'CDS':
            attributes_list = (attribute.strip().split('=') for attribute in elements[-1].split(';'))
            attributes = {key: value for key,value in attributes_list}
            if 'gene' in attributes:
              fake_gene = elements[:2] + ['gene'] + elements[3:8] + ["ID={ID};gene={gene};locus_tag={locus_tag}".format(**attributes)]
            else:
              fake_gene = elements[:2] + ['gene'] + elements[3:8] + ["ID={ID};locus_tag={locus_tag}".format(**attributes)]
            gff_out_file.write("\t".join(fake_gene) + '\n')
          gff_out_file.write(line)
      except:
        print "Problem parsing line %s: %s" % (i, line)
        raise
