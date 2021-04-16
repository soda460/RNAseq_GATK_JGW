#!/bin/bash

#$ -V
#$ -N AnalyzeCovariates
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e e
#$ -o AnalyzeCovariates_$TASK_ID.log
#$ -q all.q
#$ -t 1-24
#$ -pe smp 4

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/springer/metadata/samples.list"
TABLE_PATH="$HOME/jsb/springer/analysis/BaseQualityRecalibration"
OUTPUT="$HOME/jsb/springer/analysis/BaseQualityRecalibration"


echo $SAMPLES
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk AnalyzeCovariates \
      -bqsr $TABLE_PATH/$input"_recal_data.table" \
      -plots $OUTPUT/$input"_AnalyzeCovariates.pdf"

conda deactivate
