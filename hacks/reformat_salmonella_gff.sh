#!/bin/bash

set -x
set -e

WORKING_DIR=/lustre/scratch108/pathogen/bt5/pericles
INPUT_FILE=${WORKING_DIR}/Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1.gff
OUTPUT_FILE=${WORKING_DIR}/salmonella_typhi.gff
TEMP_FILE=${OUTPUT_FILE}.tmp

rm $OUTPUT_FILE || true
sed 's/Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1|SC|contig000001/AL513382/g' $INPUT_FILE > $OUTPUT_FILE
grep -v '^Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1|SC|contig00000' $OUTPUT_FILE > $TEMP_FILE && mv $TEMP_FILE $OUTPUT_FILE
lines_before_fasta=$(($(cat -n salmonella_typhi.gff | awk '/FASTA/ {print $1}')-1))
head -n $lines_before_fasta $OUTPUT_FILE > $TEMP_FILE && mv $TEMP_FILE $OUTPUT_FILE
gt gff3 -tidy -sort -retainids $OUTPUT_FILE > $TEMP_FILE && mv $TEMP_FILE $OUTPUT_FILE
python add_genes.py $OUTPUT_FILE $TEMP_FILE && mv $TEMP_FILE $OUTPUT_FILE
