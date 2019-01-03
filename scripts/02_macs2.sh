#!/bin/bash

for f in data/bams/*/*/*.filt.srt.dedup.q30.bam
do
	name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.filt.srt.dedup.q30.bam//g') ;
	qsub -N Macs2."$name" -e stdout/ -o stdout/ /mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/atac-seq-macs2.v1.hg38.sh $f data/peaks/$name $name ;
done

# end
