#!/bin/bash

###=================================================================================================================
### Load dependencies
###=================================================================================================================
module load homer/4.7
module load bedtools/2.23.0
module load bedops/2.4.14
module load R/3.4.1

###=================================================================================================================
### Read in config file
###=================================================================================================================

backgroundfile=$( awk -F"\t" '{if($1=="backgroundfile") print $2}' $1 )
queryfile=$( awk -F"\t" '{if($1=="queryfile") print $2}' $1)
analysisname=$( awk -F"\t" '{if($1=="analysisname") print $2}' $1)
nperm=$( awk -F"\t" '{if($1=="nperm") print $2}' $1)
scriptdir=$( awk -F"\t" '{if($1=="scriptdir") print $2}' $1)
workingdir=$( awk -F"\t" '{if($1=="workingdir") print $2}' $1)
repeatfilesdir=$( awk -F"\t" '{if($1=="repeatfilesdir") print $2}' $1)
genomesizefile=$( awk -F"\t" '{if($1=="genomesizefile") print $2}' $1)


###=================================================================================================================
### Compute Jaccard index for background files
###=================================================================================================================


for f in `seq 1 100 1500` ;
do
    k=$(($f+99));
    qsub -cwd -N jaccard."$f"."$k" -t $f-$k -q all.q -e $workingdir/stdout/$analysisname/array.jaccard/ -o $workingdir/stdout/$analysisname/array.jaccard/ $scriptdir/jaccard.arrayjob.sh $analysisname $repeatfilesdir $workingdir ;
done

qsub -cwd -N jaccard.1501.1544 -t 1501-1544 -q all.q -e $workingdir/stdout/$analysisname/array.jaccard/ -o $workingdir/stdout/$analysisname/array.jaccard/ $scriptdir/jaccard.arrayjob.sh $analysisname $repeatfilesdir $workingdir ;
