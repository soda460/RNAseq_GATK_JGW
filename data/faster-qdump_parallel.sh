#!/bin/bash
#$ -V
#$ -N fasterq-dump
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e ./qsub_err.txt
#$ -o ./qsub_log.txt
#$ -q all.q
#$ -t 1-24
#$ -pe smp 4

input=$(head -n $SGE_TASK_ID SRR.list | tail -n 1)

cd $input
fasterq-dump --split-files $input.sra
cd ..
