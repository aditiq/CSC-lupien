#!/bin/bash


## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 

## Running part of the samples.. Awaiting fastqs files for others

files="/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_13b_S27_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_21C_S26_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_36b_S28_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_64a_S32_L007_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_65c_S33_L007_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_80a_S30_L007_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_85c_S29_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180913_D00165_0222_ACCG93ANXX_Lupien_kip_Aditi.corrected/PPTO_90a_S31_L007_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_26A_S21_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_29b_S17_L005_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_2A_S24_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_30b_S18_L005_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_46a_S19_L005_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_49A_S20_L005_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_6A_S22_L006_R1_001.fastq.gz
/mnt/work1/users/lupiengroup/People/qamraa99/rawdata/180921_D00331_0345_BCCFYVANXX_Lupien_Kip_Aditi/PPTO_76A_S23_L006_R1_001.fastq.gz" 


for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1,2 ) ; 
	mkdir -p data/bams/PPTO.hiseq/"$name"
	qsub -q lupiengroup -l mem_free=8G -l h_rt=08:00:00 -N Align."$name".hiseq -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/PPTO.hiseq/"$name"/ "$name"  ;
done

