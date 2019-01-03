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
### Shuffle background keeping the chromosome and genomic location same
###=================================================================================================================

echo "#-----------------------------------------------------------------------"
echo " Scheduling jobs for generating background files"
echo "#-----------------------------------------------------------------------"
echo " "

nperm2=$(($nperm-99))

for f in `seq 1 100 $nperm2` ;
do
    k=$(($f+99));
    qsub -N bgfiles."$analysisname"."$f"."$k" -t $f-$k -q all.q -e $workingdir/stdout/$analysisname/bgfiles/ -o $workingdir/stdout/$analysisname/bgfiles/ -cwd $scriptdir/genbg.onlyexcl.sh $workingdir $genomesizefile $analysisname
done
