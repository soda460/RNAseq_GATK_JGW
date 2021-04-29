#!/bin/bash
#$ -N AnalyzeCovariates
#$ -o AnalyzeCovariates_$TASK_ID.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK_JGW/metadata/samples.list"
TABLE_PATH="$HOME/jsb/RNAseq_GATK_JGW/analysis/BaseQualityRecalibration"
OUTPUT="$HOME/jsb/RNAseq_GATK_JGW/analysis/BaseQualityRecalibration"


echo $SAMPLES
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk AnalyzeCovariates \
      -bqsr $TABLE_PATH/$input"_recal_data.table" \
      -plots $OUTPUT/$input"_AnalyzeCovariates.pdf"

conda deactivate
