#!/bin/bash
#$ -V
#$ -N 'STAR_aligner_2nd_pass'
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -b n
#$ -e e
#$ -o STAR_aligner_2nd_pass.$TASK_ID.log
#$ -q all.q
#$ -t 1-24
#$ -pe smp 4


input=$(head -n $SGE_TASK_ID SRR.list | tail -n 1)

mkdir -p ../analysis/STAR/output_2nd_pass/$input

STAR --genomeDir ./Bos_taurus.UMD3.1.87_index \
--runThreadN 4 \
--readFilesIn  ./$input/$input"_1.fastq" ./$input/$input"_2.fastq" \
--outFileNamePrefix ../analysis/STAR/output_2nd_pass/$input \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outTmpDir ../analysis/STAR/output_2nd_pass/_STARtmp_$SGE_TASK_ID \
--sjdbFileChrStartEnd \
../analysis/STAR/output/SRR5487372SJ.out.tab \
../analysis/STAR/output/SRR5487384SJ.out.tab \
../analysis/STAR/output/SRR5487396SJ.out.tab \
../analysis/STAR/output/SRR5487408SJ.out.tab \
../analysis/STAR/output/SRR5487420SJ.out.tab \
../analysis/STAR/output/SRR5487432SJ.out.tab \
../analysis/STAR/output/SRR5487376SJ.out.tab \
../analysis/STAR/output/SRR5487388SJ.out.tab \
../analysis/STAR/output/SRR5487400SJ.out.tab \
../analysis/STAR/output/SRR5487412SJ.out.tab \
../analysis/STAR/output/SRR5487424SJ.out.tab \
../analysis/STAR/output/SRR5487436SJ.out.tab \
../analysis/STAR/output/SRR5487378SJ.out.tab \
../analysis/STAR/output/SRR5487390SJ.out.tab \
../analysis/STAR/output/SRR5487402SJ.out.tab \
../analysis/STAR/output/SRR5487414SJ.out.tab \
../analysis/STAR/output/SRR5487426SJ.out.tab \
../analysis/STAR/output/SRR5487438SJ.out.tab \
../analysis/STAR/output/SRR5487382SJ.out.tab \
../analysis/STAR/output/SRR5487394SJ.out.tab \
../analysis/STAR/output/SRR5487406SJ.out.tab \
../analysis/STAR/output/SRR5487418SJ.out.tab \
../analysis/STAR/output/SRR5487430SJ.out.tab \
../analysis/STAR/output/SRR5487442SJ.out.tab



