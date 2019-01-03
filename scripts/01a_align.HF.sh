#!/bin/bash

## Called from /mnt/work1/users/lupiengroup/People/qamraa99/Pancancer.CSC

## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 

## Running part of the samples.. Awaiting fastqs files for others

files="/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160916_D00331_0203_AC9NC2ANXX_Lupien_Paul/Sample_HF6562_P9_SF_AT_1/HF6562_P9_SF_AT_1_TCCTGA_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160916_D00331_0203_AC9NC2ANXX_Lupien_Paul/Sample_HF5205_P10_SF_AT_1/HF5205_P10_SF_AT_1_AGGCAG_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160916_D00331_0203_AC9NC2ANXX_Lupien_Paul/Sample_HF7450_P9_SF_AT_1/HF7450_P9_SF_AT_1_GGACTC_L007_R1.fastq.gz" ;


for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ; 
	mkdir -p data/bams/HF/"$name"
	qsub -l mem_free=8G -l h_rt=08:00:00 -N Align.$name -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/HF/"$name"/ "$name"  ;
done
