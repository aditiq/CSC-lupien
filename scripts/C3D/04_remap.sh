#!/bin/bash
## Remap enrichment of common/shared/notshared regions identified by C3D
## mordor

module load bedtools/2.23.0

for f in results/PCSC1/C3D/finalresults/*coaccess.regions.bed
do
    name=$(echo $f |awk -F"/" '{print $NF}' | sed -e 's/.bed//g')
    echo $name $f
    for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
    do
        echo $k
        count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
        tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
        echo $tf $count >> results/PCSC1/C3D/finalresults/remap/Support.$name.txt
    done
done
