#!/bin/bash
# not run

module load deeptools

for f in data/peaks/*/*.bigwig ;
do 
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/_FE.bdg.bigwig//g')
    qsub -cwd -N plotheatmap."$name" -q all.q -e stdout/deeptools -o stdout/deeptools \
    scripts/deeptools.plotheatmap.sh $name.matrix.gz $name
done 

