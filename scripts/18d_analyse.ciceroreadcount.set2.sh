#============================================================================================================================================================
# Objective: Second set of analysis from Cicero results
#============================================================================================================================================================
    # 1. Do LSC, ds PFA share co-accessibility hubs ?
        # A. How many of them map to the same gene -- shared chromosomal structure
        # B. How many map to different genes -- different model of regulation
#============================================================================================================================================================

### mordor

#============================================================================================================================================================
## Identify E-P or P-E interactions and how many are common between CSCs
#============================================================================================================================================================

names="GBM
PFA
LSCp"

for f in $names ;
do
    awk ' BEGIN {FS=OFS="\t"} {  if (substr($5,1,4)=="ENSG") print $0, "P" ; else print $0, "E" }' results/PCSC1/cicero/run2.readcount/"$f".Deduplicated.coaccess.Annotated.txt | \
    awk ' BEGIN {FS=OFS="\t"} {  if (substr($7 ,1,4)=="ENSG") print $0, "P" ; else print $0, "E"}' | awk ' BEGIN {FS=OFS="\t"} {  print $0, $9"_"$10 }' | \
    awk ' BEGIN {FS=OFS="\t"} { if ($11=="E_P" || $11=="P_E") print $0}' |   \
    awk ' BEGIN {FS=OFS="\t"} { if ($11=="E_P") print $1,$2, $3,$7 ; else print $2,$1,$3,$5}' | awk '!seen[$0]++'  > results/PCSC1/cicero/run2.readcount/"$f".E-P.txt
done

## identify common genes present in all three or any two of the populations
nums="0
0.2
0.5
0.9
0.6
0.7
0.8"

for threshold in $nums
do
  echo $threshold

  names="GBM
  PFA
  LSCp"

  for f in $names ;
  do
      echo $f
      awk '{ if ($3 > "'${threshold}'") print $4, "'${f}'"}' results/PCSC1/cicero/run2.readcount/"$f".E-P.txt | awk '!seen[$0]++'  > "$f".cicero.tmp.txt
  done

  echo "generating files for " $threshold

  cat *.cicero.tmp.txt | awk -F" " '{a[$1]=a[$1]?a[$1]":"$2:$2;}END{for (i in a)print i, a[i];}' OFS="\t" | awk -F ":" ' { print NF-1"\t"$0 } '  > results/PCSC1/cicero/run2.readcount/E-P.coaccessgt"$threshold".genefreq.txt
  awk 'BEGIN {FS=OFS="\t"} { if ($1>1)  print $2,$3}'  results/PCSC1/cicero/run2.readcount/E-P.coaccessgt"$threshold".genefreq.txt > results/PCSC1/cicero/run2.readcount/common.coaccessgt"$threshold".E-P.genes.txt
  awk 'BEGIN {FS=OFS="\t"} { if ($1>0)  print $2,$3}'  results/PCSC1/cicero/run2.readcount/E-P.coaccessgt"$threshold".genefreq.txt > results/PCSC1/cicero/run2.readcount/shared.coaccessgt"$threshold".E-P.genes.txt
  awk 'BEGIN {FS=OFS="\t"} { if ($1==0)  print $2,$3}'  results/PCSC1/cicero/run2.readcount/E-P.coaccessgt"$threshold".genefreq.txt > results/PCSC1/cicero/run2.readcount/notshared.coaccessgt"$threshold".E-P.genes.txt

  rm *.cicero.tmp.txt
done

#============================================================================================================================================================
## Identify regions from GBM, PFA and LSC for common, shared and unique regions
#============================================================================================================================================================


#=============================================================================================
## Common genes and Shared genes
#=============================================================================================

identify_regions()
{
    name=$1
    threshold=$2
    genefile=$3
    opname=$4

    rm tmp.txt  2> /dev/null || echo > /dev/null
    rm tmp1.txt  2> /dev/null || echo > /dev/null

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } { if($3 > '${threshold}' ) print $0}' results/PCSC1/cicero/run2.readcount/"$name".E-P.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmp.txt
    cat <(cut -d" " -f2 tmp.txt ) <(cut -d" " -f3 tmp.txt ) |  sort -k1,1 -u | cut -d"_" -f1,2,3 --output-delimiter=$'\t' > tmp1.txt

    ## identify enhancers and promoters
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/run2.readcount/"$name".promoters.bed -u -f 1 > results/PCSC1/cicero/run2.readcount/"$name"."$opname".promoters.bed
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/run2.readcount/"$name".promoters.bed -v -f 1 > results/PCSC1/cicero/run2.readcount/"$name"."$opname".enhancer.bed
}

filenames="GBM
PFA
LSCp"

for j in $filenames ;
do

    identify_regions $j 0 results/PCSC1/cicero/run2.readcount/common.coaccessgt0.E-P.genes.txt "common.coaccessgt0"
    identify_regions $j 0.2 results/PCSC1/cicero/run2.readcount/common.coaccessgt0.2.E-P.genes.txt "common.coaccessgt0.2"
    identify_regions $j 0.5 results/PCSC1/cicero/run2.readcount/common.coaccessgt0.5.E-P.genes.txt "common.coaccessgt0.5"

    identify_regions $j 0 results/PCSC1/cicero/run2.readcount/shared.coaccessgt0.E-P.genes.txt "shared.coaccessgt0"
    identify_regions $j 0.2 results/PCSC1/cicero/run2.readcount/shared.coaccessgt0.2.E-P.genes.txt "shared.coaccessgt0.2"
    identify_regions $j 0.5 results/PCSC1/cicero/run2.readcount/shared.coaccessgt0.5.E-P.genes.txt "shared.coaccessgt0.5"
done


#=============================================================================================
## Regions mapping to not shared genes
#=============================================================================================

## This could also include regions that are same but map to different genes

awk -F"\t" '{print > "results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0."$2".E-P.genes.txt"}' results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.E-P.genes.txt
awk -F"\t" '{print > "results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.2."$2".E-P.genes.txt"}' results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.2.E-P.genes.txt
awk -F"\t" '{print > "results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.5."$2".E-P.genes.coaccessgt0.5.txt"}' results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.5.E-P.genes.txt


identify_unique_regions()
{
    name=$1
    threshold=$2
    genefile=$3
    opname=$4

    rm tmpns.txt  2> /dev/null || echo > /dev/null
    rm tmp1ns.txt  2> /dev/null || echo > /dev/null

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } { if($3 > '${threshold}' ) print $0}' results/PCSC1/cicero/run2.readcount/"$name".E-P.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmpns.txt
    cat <(cut -d" " -f2 tmpns.txt ) <(cut -d" " -f3 tmpns.txt ) |  sort -k1,1 -u | cut -d"_" -f1,2,3 --output-delimiter=$'\t' > tmp1ns.txt

    ## identify enhancers and promoters
    intersectBed -a tmp1ns.txt -b results/PCSC1/cicero/run2.readcount/"$name".promoters.bed -u -f 1 > results/PCSC1/cicero/run2.readcount/"$name"."$opname".promoters.bed
    intersectBed -a tmp1ns.txt -b results/PCSC1/cicero/run2.readcount/"$name".promoters.bed -v -f 1 > results/PCSC1/cicero/run2.readcount/"$name"."$opname".enhancer.bed
}

filenames="GBM
PFA
LSCp"

for j in $filenames ;
do
    echo $j
    identify_unique_regions $j 0 results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.E-P.genes.txt "notshared.coaccessgt0"
    echo "2"
    identify_unique_regions $j 0.2 results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.2.E-P.genes.txt "notshared.coaccessgt0.2"
    echo "3"
    identify_unique_regions $j 0.5 results/PCSC1/cicero/run2.readcount/notshared.coaccessgt0.5.E-P.genes.txt "notshared.coaccessgt0.5"
done
