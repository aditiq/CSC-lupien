#!/bin/bash
# not run

module load deeptools

TSSfile="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/gencode.v24.annotation.transcript.wdgeneinfo.bed"
for f in data/peaks/*/*.bigwig ;
do 
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/_FE.bdg.bigwig//g')
    qsub -cwd -N computematrix."$name" -q all.q -e stdout/deeptools -o stdout/deeptools \
    scripts/deeptools.computematrix.sh $TSSfile $f $name 3000 3000 50
done 

