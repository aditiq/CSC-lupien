#!/bin/bash

## Called from /mnt/work1/users/lupiengroup/People/qamraa99/Pancancer.CSC

## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 

for f in /mnt/work1/users/lupiengroup/Projects/RawData/JEDickProjects/AML_Fractions_ATAC/AML_Fastq/*.gz
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' ) ; 
	mkdir -p data/bams/LSC/"$name"
	qsub -l mem_free=8G -l h_rt=08:00:00 -N Align.$name -e stdout/align/ -o stdout/align/ -q light.q \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/LSC/"$name"/ "$name"  ;
done
