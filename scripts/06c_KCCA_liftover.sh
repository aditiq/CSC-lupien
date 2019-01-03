#!/bin/bash

## Liftover for great analysis


for f in results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/*Promotergrp.txt ;
do
    name=$( echo $f  | awk -F"/" '{ print $NF}' | sed -e 's/.txt//g')
    mkdir -p results/PCSC1/Cluster//KCCA.Flexclust/ExtractedClusters//liftover/
    liftOver $f ../common.data/hg38ToHg19.over.chain.gz results/PCSC1/Cluster//KCCA.Flexclust/ExtractedClusters//liftover/$name.hg19.txt results/PCSC1/Cluster/ExtractClusters.Flexclust/liftover/$name.hg19.unmap.txt 
done


for f in results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/*Enhancergrp.txt ;
do
    name=$( echo $f  | awk -F"/" '{ print $NF}' | sed -e 's/.txt//g')
    mkdir -p results/PCSC1/Cluster//KCCA.Flexclust/ExtractedClusters//liftover/
    liftOver $f ../common.data/hg38ToHg19.over.chain.gz results/PCSC1/Cluster//KCCA.Flexclust/ExtractedClusters//liftover/$name.hg19.txt results/PCSC1/Cluster/ExtractClusters.Flexclust/liftover/$name.hg19.unmap.txt 
done

for f in results/PCSC1/Cluster//KCCA.Flexclust/ExtractedClusters//liftover/*.hg19.txt ;
do 
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' $f > "$f".bed
done