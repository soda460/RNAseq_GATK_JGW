#!/bin/bash
#$ -N GenotypeGVCF
#$ -o GenotypeGVCF.log


eval "$(conda shell.bash hook)"
conda activate gatk4

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R ../SplitNCigarReads/refGenome.fasta \
   -V gendb://my_database \
   -O output.vcf.gz

conda deactivate

