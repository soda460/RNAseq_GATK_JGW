#!/bin/bash
#$ -N ApplyBQSR
#$ -o ApplyBQSR_$TASK_ID.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK/metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK/analysis/SplitNCigarReads"
OUTPUT="$HOME/jsb/RNAseq_GATK/analysis/BaseQualityRecalibration"

echo $SAMPLES
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk ApplyBQSR \
-R ./refGenome.fasta \
-I $BAMPATH/$input"_SplitNCigarReads.bam" \
--bqsr-recal-file $OUTPUT/$input"_recal_data.table" \
-O $OUTPUT/$input"_recal.bam"

conda deactivate
