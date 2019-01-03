#!/bin/bash
# mordor
# GFFcompare to identify novel transcripts in the merged gtfs

mergedfiles="./data/RNAseq/LSCp/LSCp.merged.gtf
./data/RNAseq/GBM.PFA.HF/GBM.merged.gtf
./data/RNAseq/GBM.PFA.HF/PFA.merged.gtf
./data/RNAseq/LSCp.GBM.PFA.merged.gtf" 

for f in $mergedfiles ;
do
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.merged.gtf//g' )
    echo $f > tmp."$name".txt
    ~/gffcompare-0.10.4.Linux_x86_64/gffcompare -o $name -p $name -V -r /mnt/work1/users/lupiengroup/People/qamraa99/common.data/gencode.v24.annotation.gtf -i tmp."$name".txt
done
    

