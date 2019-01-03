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


echo "#-----------------------------------------------------------------------"
echo "Starting repeat enrichment analysis"
echo "#-----------------------------------------------------------------------"


echo " "
echo "#-----------------------------------------------------------------------"
echo "Parameters being used are -- "
echo "#-----------------------------------------------------------------------"

echo " "

echo "Background file being used is " $backgroundfile
echo "Query file being used is " $queryfile
echo "Genome build for annotation is " $genomesizefile
echo "All results and stdout files will be generated at " $workingdir
echo "All repeat files are being read from " $repeatfilesdir
echo "Results will be found in " $analysisname " folder in results directory of the specified working directory "
echo $nperm " no of permutations will be run for enrichment analysis "
echo "All scripts for analysis are kept at " $scriptdir
echo " "

###=================================================================================================================
### Make directories
###=================================================================================================================

echo "#-----------------------------------------------------------------------"
echo "Generating results and stdout directories at " $workingdir
echo "#-----------------------------------------------------------------------"
echo " "

mkdir -p $workingdir/results/$analysisname/jaccardfiles $workingdir/results/$analysisname/bgfiles.onlyexcl
mkdir -p $workingdir/stdout/$analysisname/array.jaccard $workingdir/stdout/$analysisname/bgfiles

##=================================================================================================================
# Query file jaccard index calculation
##=================================================================================================================

sortBed -i $queryfile  > $queryfile.sorted 

for f in `seq 1 1544` ;
do
    echo $f
    bedtools jaccard -a $repeatfilesdir/"$f"_* -b $queryfile.sorted | \
    awk -F"\t" '{ if (NR>1) print $3, "'${f}'"}' >> $workingdir/results/$analysisname/jaccardfiles/queryfile.bed
done

###=================================================================================================================
### Annotate and split files
###=================================================================================================================

echo "#-----------------------------------------------------------------------------------------------------"
echo "Annotating and splitting background and query files by genomic location at " $workingdir"/results/"
echo "#-----------------------------------------------------------------------------------------------------"
echo " "

Rscript $scriptdir/annotate.R $backgroundfile $workingdir/results/"$analysisname"/backgroundfile.annotated
Rscript $scriptdir/annotate.R $queryfile $workingdir/results/"$analysisname"/queryfile.annotated

## Split by annotation and chromosome
grep -v chrM $workingdir/results/$analysisname/backgroundfile.annotated | grep -v random | grep -v chrUn  | grep -v alt | grep -v KI27 | grep -v chrGL | awk -F'\t' '{ if (NR>1) print $1,$2,$3, $5}'  | awk 'BEGIN { FS="(" ; OFS="\t"} { print $1}' | awk -F' '  '{ print $1"\t"$2"\t"$3"\t"$1"."$4} ' | sed -e 's/'\''/prime/g' | awk -F"\t" '{print >  "'${workingdir}'""/results/""'${analysisname}'""/backgroundfile.anno."$4}'
grep -v chrM $workingdir/results/$analysisname/queryfile.annotated | grep -v random | grep -v chrUn  |   grep -v alt | grep -v KI27 |   grep -v chrGL | awk -F'\t' '{ if (NR>1) print $1,$2,$3, $5}'  | awk 'BEGIN { FS="(" ; OFS="\t"} { print $1}' | awk -F' '  '{ print $1"\t"$2"\t"$3"\t"$1"."$4} ' | sed -e 's/'\''/prime/g' | awk -F"\t" '{print >  "'${workingdir}'""/results/""'${analysisname}'""/queryfile.anno."$4}'

###=================================================================================================================
### Shuffle background keeping the chromosome and genomic location same
###=================================================================================================================

# echo "#-----------------------------------------------------------------------"
# echo " Scheduling jobs for generating background files"
# echo "#-----------------------------------------------------------------------"
# echo " "

# nperm2=$(($nperm-99))

# for f in `seq 1 100 $nperm2` ;
# do
#     k=$(($f+99));
#     qsub -N bgfiles."$f"."$k" -t $f-$k -q all.q -e $workingdir/stdout/$analysisname/bgfiles/ -o $workingdir/stdout/$analysisname/bgfiles/ -cwd $scriptdir/genbg.onlyexcl.sh $workingdir $genomesizefile $analysisname;
# done
