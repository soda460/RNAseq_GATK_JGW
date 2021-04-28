#!/bin/bash
#$ -N BaseRecalibrator
#$ -o BaseQualityRecalibration_$TASK_ID.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK/metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK/analysis/SplitNCigarReads"
OUTPUT="$HOME/jsb/RNAseq_GATK/analysis/BaseQualityRecalibration"

input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk BaseRecalibrator \
-R ./refGenome.fasta \
-I $BAMPATH/$input"_SplitNCigarReads.bam" \
--known-sites ./known.vcf \
-O $OUTPUT/$input"_recal_data.table"

conda deactivate
