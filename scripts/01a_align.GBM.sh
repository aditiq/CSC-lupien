#!/bin/bash


## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 
files="/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G361/G361_v2.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G411/G411.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G551/G551.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G564/G564.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G566/G566.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G719/G719.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G729/G729.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G799/G799.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/151211_SN1068_0259_BC883RACXX_Paul/Sample_G489_P14_SF_AT_1/G489_P14_SF_AT_1_AGGCAGAA_L003_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160330_SN1080_0266_AC8U8JACXX_Paul/Sample_G549_P8_SF_AT_1/G549_P8_SF_AT_1_AGGCAGAA_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/151211_SN1068_0259_BC883RACXX_Paul/Sample_G571_P11_SF_AT_1/G571_P11_SF_AT_1_TAGGCATG_L003_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160122_SN1080_0260_BC8MM3ACXX_Paul/Sample_G577_P10_SF_AT_1/G577_P10_SF_AT_1_GGACTCCT_L002_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160330_SN1080_0266_AC8U8JACXX_Paul/Sample_G584_P8_SF_AT_1/G584_P8_SF_AT_1_TAGGCATG_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160419_SN1080_0270_BC98HFACXX_Lupien_Paul/Sample_G594_P8_SF_AT_1/G594_P8_SF_AT_1_AGGCAG_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160930_D00343_0140_AC9NGFANXX_Lupien_Paul/Sample_G613_P12_SF_AT_2/G613_P12_SF_AT_2_GGACTCCT_L005_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_G620_P8_SF_AT_2/G620_P8_SF_AT_2_TCCTGAGC_L006_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/151211_SN1068_0259_BC883RACXX_Paul/Sample_G702_P11_SF_AT_1/G702_P11_SF_AT_1_GGACTCCT_L004_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160122_SN1080_0260_BC8MM3ACXX_Paul/Sample_G705_P9_SF_AT_1/G705_P9_SF_AT_1_AGGCAGAA_L003_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160419_SN1080_0270_BC98HFACXX_Lupien_Paul/Sample_G706_P8_SF_AT_1/G706_P8_SF_AT_1_GGACTC_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160122_SN1080_0260_BC8MM3ACXX_Paul/Sample_G744_P8_SF_AT_1/G744_P8_SF_AT_1_TAGGCATG_L002_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160122_SN1080_0260_BC8MM3ACXX_Paul/Sample_G797_P9_SF_AT_2/G797_P9_SF_AT_2_TAGGCATG_L003_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160419_SN1080_0270_BC98HFACXX_Lupien_Paul/Sample_G800_P8_SF_AT_1/G800_P8_SF_AT_1_TAGGCA_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160122_SN1080_0260_BC8MM3ACXX_Paul/Sample_G498_P8_SF_AT_1/G498_P8_SF_AT_1_AGGCAGAA_L002_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_G523_P9_SF_AT_2/G523_P9_SF_AT_2_TCCTGAGC_L005_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160330_SN1080_0266_AC8U8JACXX_Paul/Sample_G567_P8_SF_AT_1/G567_P8_SF_AT_1_GGACTCCT_L008_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/160930_D00343_0140_AC9NGFANXX_Lupien_Paul/Sample_G583_P11_SF_AT_3/G583_P11_SF_AT_3_TCCTGAGC_L005_R1.fastq.gz
/mnt/work1/users/lupiengroup/People/Paul/FASTQ/151211_SN1068_0259_BC883RACXX_Paul/Sample_G648_P14_SF_AT_1/G648_P14_SF_AT_1_AGGCAGAA_L004_R1.fastq.gz" ;


for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ; 
	mkdir -p data/bams/GBM/"$name"
	qsub -l mem_free=8G -l h_rt=08:00:00 -N Align.$name -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/GBM/"$name"/ "$name"  ;
done
