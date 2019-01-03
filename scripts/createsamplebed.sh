#!/bin/bash

module load bedops/2.4.14   

$1=inputfile
$2=samplenum
$3=opdir
$4=opname



sort-bed $inputfile > elements.bed
 for chr in `bedextract --list-chr elements.bed`; 
 do 
    bedextract $chr elements.bed > elements.$chr.bed; 
    count=`wc -l elements.$chr.bed | cut -d' ' -f1`; 
    echo -e "$chr\t$count" >> counts.txt; 
done 

sum=`cut -f2 counts.txt | perl -nle '$sum += $_ } END { print $sum' -`
awk -v sum=$sum -v sampleSize=$samplenum '{print $0"\t"($2/sum)"\t"int(sampleSize*$2/sum)}' counts.txt > proportions.txt

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

