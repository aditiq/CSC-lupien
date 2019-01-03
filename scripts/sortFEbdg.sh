#!/bin/bash


for f in data/peaks/*/*_FE.bedgraph
do
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/_FE.bedgraph//g' ) ;
    qsub -q light.q -e stdout/sortfebdg/ -o stdout/sortfebdg/ -N sortFEbdg."$name" /mnt/work1/users/lupiengroup/People/qamraa99/common.scripts/sortbdg.sh $f data/peaks/$name/"$name"_sortedFE.bdg
done
