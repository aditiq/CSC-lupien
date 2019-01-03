#!/bin/bash

module load bedops/2.4.14   

$1=inputfile
$2=samplenum
$3=opdir
$4=opname

## sort file
sort-bed $inputfile > elements.bed

## Split file by last column and count no. of lines 

for state in `cut -f4 elements.bed | sort -u`; 
do 
    awk 'BEGIN {FS=OFS="\t"} { if ($4=="'${state}'") print $0}' elements.bed > elements.$state.bed; 
    count=`wc -l elements.$state.bed | cut -d' ' -f1`; 
    echo -e "$state\t$count" >> counts.txt; 
done 

## total no. of records and % of each chr in that to get proportion
sum=`cut -f2 counts.txt | perl -nle '$sum += $_ } END { print $sum' -`
awk -v sum=$sum -v sampleSize=$samplenum '{print $0"\t"($2/sum)"\t"int(sampleSize*$2/sum)+1}' counts.txt > proportions.txt

iter=0

while [ "$iter" -lt 100 ]
do

echo $iter

awk '{ \
    perChrName=$1; \
    perChrN=$4; \
    cmd="shuf elements."perChrName".bed | head -n "perChrN; \
    system(cmd); \
}' proportions.txt \
| sort-bed - > $opdir/"$opname".sample."$iter".bed

iter=`expr $iter + 1`

done

rm elements*bed counts.txt  proportions.txt

