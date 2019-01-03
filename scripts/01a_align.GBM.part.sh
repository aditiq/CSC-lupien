#!/bin/bash

## Rerunning from different steps for different  files


## Arguments
#1 Original.fastq.gz file
#2 Output dir 
#3 Output name
 
files="/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G361/G361_v2.fastq.gz"


for  f in $files ;
do
	name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ; 
	mkdir -p data/bams/GBM/"$name"
	qsub -l mem_free=15G -l h_rt=24:00:00 -q all.q -N Rerun.Align.$name -e stdout/align/ -o stdout/align/ \
	/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.sh $f data/bams/GBM/"$name"/ "$name"  ;
done



files="/mnt/work1/users/lupiengroup/People/Paul/GBM/ATAC/G719/G719.fastq.gz"

for  f in $files ;
do
        name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ;
        mkdir -p data/bams/GBM/"$name"
        qsub -l mem_free=15G -l h_rt=24:00:00 -q all.q -N Rerun.Align.$name -e stdout/align/ -o stdout/align/ \
        /mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.part1.sh $f data/bams/GBM/"$name"/ "$name"  ;
done



files="/mnt/work1/users/lupiengroup/People/Paul/FASTQ/161221_D00343_0154_ACAAE8ANXX_Lupien_Paul/Sample_G523_P9_SF_AT_2/G523_P9_SF_AT_2_TCCTGAGC_L005_R1.fastq.gz" ;

for  f in $files ;
do
        name=$( echo $f | awk -F"/" '{print $NF}' | sed -e 's/.fastq.gz//g' | cut -d"_" -f1 ) ;
        mkdir -p data/bams/GBM/"$name"
        qsub -l mem_free=15G -l h_rt=24:00:00 -q all.q -N Rerun.Align.$name -e stdout/align/ -o stdout/align/ \
        /mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-align.v1.hg38.part2.sh $f data/bams/GBM/"$name"/ "$name"  ;
done

