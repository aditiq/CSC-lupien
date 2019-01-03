#!/bin/bash

module load  bedtools/2.23.0

workingdir=$1
genomesizefile=$2
analysisname=$3
k=$SGE_TASK_ID

echo $k

for f in $workingdir/results/$analysisname/queryfile.anno.* ;
do

    anno=$(echo $f | awk -F/ '{ print $NF}' | sed -e 's/queryfile.anno.//g' )

    bedtools complement -i $workingdir/results/$analysisname/backgroundfile.anno.$anno -g $genomesizefile | shuffleBed -i $f  -f 0.0 -excl stdin -g $genomesizefile -maxTries 10000 -chrom -seed $k -noOverlapping >> $workingdir/results/$analysisname/bgfiles.onlyexcl/Bg.$k  ;

done;
