#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N 'fasterq-dump'
#$ -pe smp 12
#$ -o ./qsub_log.txt
#$ -e ./qsub_err.txt
	
for i in `ls -d SRR*`; do
	cd $i
	fasterq-dump --split-files $i.sra
	cd ..
done
