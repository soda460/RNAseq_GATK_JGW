#!/bin/bash
#$ -N GVCF
#$ -o singleSampleCallingGVCF_$TASK_ID.log

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/RNAseq_GATK_JGW/metadata/samples.list"
BAMPATH="$HOME/jsb/RNAseq_GATK_JGW/analysis/BaseQualityRecalibration"
OUTPUT="$HOME/jsb/RNAseq_GATK_JGW/analysis/GVCF"


echo $SAMPLES
input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

gatk --java-options "-Xmx4g" HaplotypeCaller \
	 --dont-use-soft-clipped-bases true \
     -R ../SplitNCigarReads/refGenome.fasta \
     -I $BAMPATH/$input"_recal.bam" \
     -O $OUTPUT/$input".g.vcf.gz" \
     -ERC GVCF

conda deactivate

