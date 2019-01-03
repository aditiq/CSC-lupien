#!/bin/bash

module load deeptools/2.4.2    
for f in  data/bams/*/keep/*/*.filt.srt.dedup.q30.bam ; 
do
    path=$(echo $f)
    name=$(echo $f | awk -F"/" '{ print $NF}' | awk -F"." '{ print $1}' )
    
    cmd1="module load deeptools/2.4.2 ; bamCoverage -b \"$path\" -o \"$path\".bw -of bigwig --binSize 1 --blackListFileName /mnt/work1/users/lupiengroup/Projects/CommonData/hg19.blacklist.bed --normalizeUsingRPKM -ignore chrM chrY chrX"
    
    echo $cmd1 > scripts/bamcoverage/BamCoverage."$name".sh
    sed -i 's/"//g' scripts/bamcoverage/BamCoverage."$name".sh

    qsub -l mem_free=30G -l h_rt=08:00:00 -q highmem.q -N BamCoverage.$name -e stdout/bamcoverage -o stdout/bamcoverage -cwd  scripts/bamcoverage/BamCoverage."$name".sh

done




