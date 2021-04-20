#!/bin/bash
#$ -V
#$ -N 'STAR_aligner'
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e e
#$ -o STAR_aligner_1st_pass.$TASK_ID.log
#$ -q all.q
#$ -t 1-24
#$ -pe smp 4  


input=$(head -n $SGE_TASK_ID SRR.list | tail -n 1)


mkdir -p ../analysis/STAR/output/$input

STAR --genomeDir ./Bos_taurus.UMD3.1.87_index \
--runThreadN 12 \
--readFilesIn  ./$input/$input"_1.fastq" ./$input/$input"_2.fastq" \
--outFileNamePrefix ../analysis/STAR/output/$input \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard






