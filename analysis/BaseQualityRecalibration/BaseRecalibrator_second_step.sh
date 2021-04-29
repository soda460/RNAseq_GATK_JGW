#!/bin/bash
#$ -N BaseRecalibrator
#$ -o BaseQualityRecalibration_$TASK_ID.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK_JGW/metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK_JGW/analysis/BaseQualityRecalibration"
OUTPUT="$HOME/jsb/RNAseq_GATK_JGW/analysis/BaseQualityRecalibration"
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk BaseRecalibrator \
-R ../SplitNCigarReads/refGenome.fasta \
-I $BAMPATH/$input"_recal.bam" \
--known-sites ./snp50.vcf.gz \
-O $OUTPUT/$input"_recal_data2.table"

conda deactivate
