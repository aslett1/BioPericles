#!/bin/bash

set -x
set -e

WORKING_DIRECTORY=/lustre/scratch108/pathogen/bt5/pericles
VEP=/nfs/users/nfs_j/jl11/installations/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl 
VCF=${WORKING_DIRECTORY}/output/all-snp.vcf
OUTPUT_FILE=${WORKING_DIRECTORY}/output/all-snp-variations.txt

cd $WORKING_DIRECTORY
export COLUMNS=80
export LINES=100
perl $VEP -i $VCF --cache --species CT18 --cache_version 25 --offline --output_file $OUTPUT_FILE -coding_only
