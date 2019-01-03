#!/bin/bash

## Split clusters

for f in results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/FlexClust.Enhancer.Groups.txt  ;
do
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.txt//g' )
    awk ' BEGIN {FS="_"} { if (NR>1) print $1"\t"$2"\t"$3}' $f | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3,$5}' | sort -k4,4 | awk -F"\t" '{print > $4".Enhancergrp.txt"}' ;
done

for f in results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/FlexClust.Promoter.Groups.txt ;
do
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.txt//g' )
    awk ' BEGIN {FS="_"} { if (NR>1) print $1"\t"$2"\t"$3}' $f | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3,$5}' | sort -k4,4 | awk -F"\t" '{print > $4".Promotergrp.txt"}' ;
done
