#!/bin/bash

for f in /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/LSCp/*.gz ;
do 
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.fastq.gz//g' )
    cmd="module load fastqc/0.11.5  ; fastqc $f" ;
    echo $cmd > scripts/rnaseq/"$name".fastqc.sh
    sed -i 's/"//g' scripts/rnaseq/"$name".fastqc.sh
    qsub -q all.q -N LSC.fastqc.$name -e stdout/rnaseq/ -o stdout/rnaseq/ -cwd scripts/rnaseq/"$name".fastqc.sh
done

