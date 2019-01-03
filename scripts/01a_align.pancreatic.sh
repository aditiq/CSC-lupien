#!/bin/bash


## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 

## Running part of the samples.. Awaiting fastqs files for others

files="/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_29b_S5_L001_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_29b_S5_L002_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_29b_S5_L003_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_29b_S5_L004_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_30b_S6_L001_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_30b_S6_L002_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_30b_S6_L003_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_30b_S6_L004_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_46a_S7_L001_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_46a_S7_L002_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_46a_S7_L003_R1_001.fastq.gz,/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/PPTO_46a_S7_L004_R1_001.fastq.gz"

for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1,2 ) ; 
	mkdir -p data/bams/PPTO/"$name"
	qsub -l mem_free=8G -l h_rt=08:00:00 -N Align.$name -q all.q -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/PPTO/"$name"/ "$name"  ;
done
