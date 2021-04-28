#!/bin/bash
#$ -N SplitNCigarReads
#$ -o SplitNCigarReads.$TASK_ID.log
	
eval "$(conda shell.bash hook)"
conda activate gatk4
	
SAMPLES="$HOME/jsb/RNAseq_GATK_JGW//metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK_JGW/analysis/MarkDuplicates"
OUTPUT="$HOME/jsb/RNAseq_GATK_JGW/analysis/SplitNCigarReads"
	
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk SplitNCigarReads \
-R refGenome.fasta \
-I $BAMPATH/$input"_marked_duplicates.bam" \
-O $OUTPUT/$input"_SplitNCigarReads.bam" \
--tmp-dir $OUTPUT/gatk_tmp
	
conda deactivate
