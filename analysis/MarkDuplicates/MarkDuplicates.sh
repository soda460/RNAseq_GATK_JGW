#!/bin/bash
#$ -N MarkDuplicates
#$ -o logfile.$TASK_ID.log

SAMPLES="$HOME/jsb/RNAseq_GATK_JGW/metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK_JGW/analysis/addReadGroups"
OUTPUT="$HOME/jsb/RNAseq_GATK_JGW/analysis/MarkDuplicates"

input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

java -jar $PICARD MarkDuplicates \
I=$BAMPATH/$input".bam" \
O=$OUTPUT/$input"_marked_duplicates.bam" \
M=$OUTPUT/$input"_marked_dup_metrics.txt"
