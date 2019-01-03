#!/bin/sh

## Array job to run jaccard index claculation for repeat element enrichment
## Will be running 11 Array jobs - 10 with 100 jobs each and 1 with 25 jobs ..
## So each job within the array will be for one repeat file. All repeat files end with a unique number between 1 and 1544
## In each array job, jaccard index will be calculated for SGE_TASK_ID repeat files on 10K permutations

module load  bedtools/2.23.0  ## Don't use 2.26 because shuffle and jaccard bed are much slower

analysisname=$1
repeatfilesdir=$2
workingdir=$3

name=$(echo $repeatfilesdir/"$SGE_TASK_ID"_* | awk -F"/" '{ print $NF}' |  sed -e 's/.bed.sorted.bed//g' )

for k in `seq 1 1000` ;
do
    echo " Running permutation number \"$k\" for \"$name\" using Bg.\"$k\" ";
    sortBed -i $workingdir/results/$analysisname/bgfiles.onlyexcl/Bg."$k" | bedtools jaccard -a $repeatfilesdir/"$SGE_TASK_ID"_* -b stdin | awk -F"\t" '{ if (NR>1) print $3, "'${k}'"}' >> $workingdir/results/$analysisname/jaccardfiles/Bg.Jaccard."$name".bed
done
