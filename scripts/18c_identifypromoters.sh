#!/bin/bash

#============================================================================================================================================================
# Objective: Identify promoters from co-accessible regions and get list of deduped pairs from Cicero
#============================================================================================================================================================

files="GBM
PFA
LSCp
PCSC1"

module load bedtools/2.23.0
for f in $files;
do
    rm tmp1.txt  2> /dev/null || echo > /dev/null
    rm tmp2.txt  2> /dev/null || echo > /dev/null
    awk 'BEGIN {FS="\t" } { if (NR>1) print $1}' results/PCSC1/cicero/"$f".connsfull.bed > tmp1.txt
    awk 'BEGIN {FS="\t" } { if (NR>1) print $2}' results/PCSC1/cicero/"$f".connsfull.bed > tmp2.txt

    cat tmp1.txt tmp2.txt | awk '!x[$0]++' | awk 'BEGIN {FS="_"} { print $1"\t"$2"\t"$3}' | \
    closestBed -a stdin  -b ../common.data/gencode.v24.annotation.TSS.wdgeneinfo.sorted.bed -d | \
    awk 'BEGIN {FS=OFS="\t" }{ if ($NF<=500) print $1,$2,$3,$8,$9,$10,$7,$11, $1"_"$2"_"$3}' > results/PCSC1/cicero/"$f".promoters.bed

    awk ' BEGIN {FS=OFS="\t"} { if(NR>1) print $1,$2}' results/PCSC1/cicero/"$f".connsfull.bed | \
    awk '!seen[$1>$2 ? $1 FS $2 : $2 FS $1]++' > results/PCSC1/cicero/"$f".deduped.connsfull.pairs.bed

done

files="GBM
PFA
LSCp"

module load bedtools/2.23.0
for f in $files;
do
    rm tmp1.txt  2> /dev/null || echo > /dev/null
    rm tmp2.txt  2> /dev/null || echo > /dev/null
    awk 'BEGIN {FS="\t" } { if (NR>1) print $1}' results/PCSC1/cicero/run2.readcount/"$f".connsfull.bed > tmp1.txt
    awk 'BEGIN {FS="\t" } { if (NR>1) print $2}' results/PCSC1/cicero/run2.readcount/"$f".connsfull.bed > tmp2.txt

    cat tmp1.txt tmp2.txt | awk '!x[$0]++' | awk 'BEGIN {FS="_"} { print $1"\t"$2"\t"$3}' | sortBed -i stdin | \
    closestBed -a stdin  -b ../common.data/gencode.v24.annotation.TSS.wdgeneinfo.sorted.bed -d | \
    awk 'BEGIN {FS=OFS="\t" }{ if ($NF<=500) print $1,$2,$3,$8,$9,$10,$7,$11, $1"_"$2"_"$3}' > results/PCSC1/cicero/run2.readcount/"$f".promoters.bed

    awk ' BEGIN {FS=OFS="\t"} { if(NR>1) print $1,$2}' results/PCSC1/cicero/run2.readcount/"$f".connsfull.bed | \
    awk '!seen[$1>$2 ? $1 FS $2 : $2 FS $1]++' > results/PCSC1/cicero/run2.readcount/"$f".deduped.connsfull.pairs.bed

done
