#!/bin/bash
#$ -S /bin/bash
#$ -cwd

for i in `ls SRR5487400Aligned*.bam | xargs basename -s '.bam'`; do
samtools sort $i.bam -o $i"_sorted.bam"
done 
