#!/bin/bash
#$ -N samtools_index
#$ -o samtools_index.$TASK_ID.log

input=$(head -n $SGE_TASK_ID bam.list | tail -n 1)

samtools index $input
