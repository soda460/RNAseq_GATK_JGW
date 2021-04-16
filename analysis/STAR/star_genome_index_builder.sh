#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N 'STAR_genome_builder'
#$ -pe smp 6
#$ -o ./qsub_log.txt
#$ -e ./qsub_err.txt 

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir Bos_taurus.UMD3.1.87_index \
--genomeFastaFiles ./Bos_taurus.UMD3.1.dna.toplevel.fa \
--sjdbGTFfile ./Bos_taurus.UMD3.1.87.gtf \
--sjdbOverhang 99
