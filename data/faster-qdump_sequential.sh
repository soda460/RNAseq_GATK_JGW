#! /bin/bash
#$ -N 'fasterq-dump'
#$ -o ./faster-qdump_log.txt
	
for i in `ls -d SRR*`; do
	cd $i
	fasterq-dump --split-files $i.sra
	cd ..
done
