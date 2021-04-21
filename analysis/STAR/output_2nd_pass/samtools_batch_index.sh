#!/bin/bash
#$ -V
#$ -N samtools_sort
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e e
#$ -o o
#$ -q all.q
#$ -t 1-24
#$ -pe smp 2

input=$(head -n $SGE_TASK_ID bam.list | tail -n 1)

samtools sort $input

