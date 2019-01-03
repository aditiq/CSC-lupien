#!/bin/bash

for f in data/bams/PPTO/PPTO_*trim//*.filt.srt.dedup.q30.bam
do
	name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.filt.srt.dedup.q30.bam//g') ;
	qsub -N Macs2."$name" -e stdout/macs/ -o stdout/macs/  /mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-macs2.v1.hg38.keepdup1.sh $f data/peaks/$name $name ;
done

