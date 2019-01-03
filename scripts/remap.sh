#!/bin/bash
## This scripts needs cleaning up
## Intersect with Remap

module load bedtools/2.26.0

for f in results/PCSC1/Cluster/ExtractClusters.Flexclust/*Enhancergrp.txt ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.Enhancergrp.txt//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Enhancer/Support.$name.txt
    done
done


for f in results/PCSC1/Cluster/ExtractClusters.Flexclust/*Promotergrp.txt ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.Promotergrp.txt//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Promoter/Support.$name.txt
    done
done

for f in data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancers.bed  ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' )
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Enhancer/Support.$name.txt
    done
done


for f in data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoters.bed  ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' )
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Promoter/Support.$name.txt
    done
done

## scABC

for f in results/PCSC1/Cluster/scABC/*Enhancer.p0.05.bed ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.Enhancer.p0.05.bed//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/Cluster/scABC/remap/Enhancer/Support.$name.txt
    done
done


for f in results/PCSC1/Cluster/scABC/*Promoter.p0.05.bed ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.Promoter.p0.05.bed//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/Cluster/scABC/remap/Promoter/Support.$name.txt
    done
done

##  Cicero -

for f in results/PCSC1/cicero/run2.readcount/*coaccess*bed
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.bed//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/cicero//run2.readcount/remap/Support.$name.txt
    done
done

for f in data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Enhancers.bed  ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' )
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/cicero//run2.readcount/remap//Support.$name.txt
    done
done


for f in data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Promoters.bed  ;
do
    name=$(echo $f |awk -F"/" '{print $NF}' )
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/cicero//run2.readcount/remap/Support.$name.txt
    done
done
