#!/bin/bash
#$ -N GenomicsDBImport
#$ -o GenomicsDBImport.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK/metadata"
OUTPUT="$HOME/jsb/RNAseq_GATK/analysis/GenomicsDBImport"

gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $OUTPUT/"my_database" \
       -L 1 \
       -L 2 \
       -L 3 \
       -L 4 \
       -L 5 \
       -L 6 \
       -L 7 \
       -L 8 \
       -L 9 \
      -L 10 \
      -L 11 \
      -L 12 \
      -L 13 \
      -L 14 \
      -L 15 \
      -L 16 \
      -L 17 \
      -L 18 \
      -L 19 \
      -L 20 \
      -L 21 \
      -L 22 \
      -L 23 \
      -L 24 \
      -L 25 \
      -L 26 \
      -L 27 \
      -L 28 \
      -L 29 \
       --sample-name-map $SAMPLES/"cohort_1.sample_map" \
       --tmp-dir $HOME/jsb/tmp \
       --reader-threads 5
conda deactivate

