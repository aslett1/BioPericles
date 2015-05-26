#!/bin/bash

set -x
set -e

WORKING_DIRECTORY=/nfs/users/nfs_b/bt5/BioPericles
DATA_DIRECTORY=/lustre/scratch108/pathogen/bt5/pericles
VCF_INPUT=${DATA_DIRECTORY}/output/all-snp.vcf
GFF_INPUT=${DATA_DIRECTORY}/Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1.gff
VCF_OUTPUT=${DATA_DIRECTORY}/output/all-snp-with-GT.vcf

cd $WORKING_DIRECTORY
python hacks/annotate-vcf.py $VCF_INPUT $GFF_INPUT $VCF_OUTPUT
echo "Done"
