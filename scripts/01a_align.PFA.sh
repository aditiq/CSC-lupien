#!/bin/bash

## Called from /mnt/work1/users/lupiengroup/People/qamraa99/Pancancer.CSC

## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 

## Running part of the samples.. Awaiting fastqs files for others

files="/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA8_SF_AT_3/PFA8_SF_AT_3_GGACTCCT_L006_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA6_SF_AT_3/PFA6_SF_AT_3_CTCTCTAC_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA3_SF_AT_3/PFA3_SF_AT_3_GGACTCCT_L005_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA2_SF_AT_3/PFA2_SF_AT_3_AGGCAGAA_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA7_SF_AT_3/PFA7_SF_AT_3_CAGAGAGG_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA5_SF_AT_3/PFA5_SF_AT_3_TAAGGCGA_L007_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_PFA4_SF_AT_4/PFA4_SF_AT_4_TAGGCATG_L006_R1.fastq.gz" ;


for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ; 
	mkdir -p data/bams/PFA/"$name"
	qsub -l mem_free=8G -l h_rt=08:00:00 -N Align.$name -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/PFA/"$name"/ "$name"  ;
done
