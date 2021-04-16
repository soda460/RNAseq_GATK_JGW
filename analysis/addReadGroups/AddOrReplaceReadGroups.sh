#!/bin/bash

#$ -V
#$ -N AddOrReplaceReadGroups
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e e
#$ -o logfile.txt
#$ -q all.q
#$ -t 1-24
#$ -pe smp 1

eval "$(conda shell.bash hook)"
conda activate gatk4

SAMPLES="$HOME/jsb/springer/metadata/samples.list"
BAMPATH="$HOME/jsb/springer/analysis/STAR/output_2nd_pass"
OUTPUT="$HOME/jsb/springer/analysis/addReadGroups"

readarray -t RGLB < ./RGLB.txt
readarray -t RGSM < ./RGSM.txt
readarray -t RGPU < ./RGPU.txt

input=$(head -n $SGE_TASK_ID $SAMPLES | tail -n 1)

java -jar $PICARD AddOrReplaceReadGroups \
       I=$BAMPATH/$input"Aligned.sortedByCoord.out.bam" \
       O=$OUTPUT/$input".bam" \
       RGLB=$RGLB \
       RGPL=ILLUMINA \
       RGPU=$RGPU \
       RGSM=$input

conda deactivate
