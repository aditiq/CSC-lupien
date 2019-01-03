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
    awk ' BEGIN {FS=OFS="\t"} {  if (substr($5,1,4)=="ENSG") print $0, "P" ; else print $0, "E" }' results/PCSC1/cicero/"$f".Deduplicated.coaccess.Annotated.txt | \
    awk ' BEGIN {FS=OFS="\t"} {  if (substr($7 ,1,4)=="ENSG") print $0, "P" ; else print $0, "E"}' | awk ' BEGIN {FS=OFS="\t"} {  print $0, $9"_"$10 }' | \
    awk ' BEGIN {FS=OFS="\t"} { if ($11=="E_P" || $11=="P_E") print $0}' |   \
    awk ' BEGIN {FS=OFS="\t"} { if ($11=="E_P") print $1,$2, $3,$7 ; else print $2,$1,$3,$5}' | awk '!seen[$0]++'  > results/PCSC1/cicero/"$f".E-P.txt
done

## identify common genes present in all three or any two of the populations
for f in $names ;
do
    awk '{print $4, "'${f}'"}' results/PCSC1/cicero/"$f".E-P.txt | awk '!seen[$0]++'  > "$f".cicero.tmp.txt
done

cat *.cicero.tmp.txt | awk -F" " '{a[$1]=a[$1]?a[$1]":"$2:$2;}END{for (i in a)print i, a[i];}' OFS="\t" | awk -F ":" ' { print NF-1"\t"$0 } '  > results/PCSC1/cicero/E-P.genefreq.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>1)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.txt > results/PCSC1/cicero/common.E-P.genes.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.txt > results/PCSC1/cicero/shared.E-P.genes.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1==0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.txt > results/PCSC1/cicero/notshared.E-P.genes.txt

rm *.cicero.tmp.txt

## Higher co-accessibility
for f in $names ;
do
    awk '{ if ($3 >0.2) print $4, "'${f}'"}' results/PCSC1/cicero/"$f".E-P.txt | awk '!seen[$0]++' > "$f".cicero.tmp.txt
done

cat *.cicero.tmp.txt | awk -F" " '{a[$1]=a[$1]?a[$1]":"$2:$2;}END{for (i in a)print i, a[i];}' OFS="\t" | awk -F ":" ' { print NF-1"\t"$0 }'  > results/PCSC1/cicero/E-P.genefreq.coaccessgt0.2.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>1)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.2.txt > results/PCSC1/cicero/common.E-P.genes.coaccessgt0.2.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.2.txt > results/PCSC1/cicero/shared.E-P.genes.coaccessgt0.2.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1==0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.2.txt > results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.2.txt

rm *.cicero.tmp.txt


for f in $names ;
do
    awk '{ if ($3 >0.5) print $4, "'${f}'"}' results/PCSC1/cicero/"$f".E-P.txt | sort -k1,1 -u > "$f".cicero.tmp.txt
done

cat *.cicero.tmp.txt | awk -F" " '{a[$1]=a[$1]?a[$1]":"$2:$2;}END{for (i in a)print i, a[i];}' OFS="\t" | awk -F ":" ' { print NF-1"\t"$0 }'  > results/PCSC1/cicero/E-P.genefreq.coaccessgt0.5.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>1)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.5.txt > results/PCSC1/cicero/common.E-P.genes.coaccessgt0.5.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.5.txt > results/PCSC1/cicero/shared.E-P.genes.coaccessgt0.5.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1==0)  print $2,$3}'  results/PCSC1/cicero/E-P.genefreq.coaccessgt0.5.txt > results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.5.txt

rm *.cicero.tmp.txt

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

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } { if($3 > '${threshold}' ) print $0}' results/PCSC1/cicero/"$name".E-P.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmp.txt
    cat <(cut -d" " -f2 tmp.txt ) <(cut -d" " -f3 tmp.txt ) |  sort -k1,1 -u | cut -d"_" -f1,2,3 --output-delimiter=$'\t' > tmp1.txt

    ## identify enhancers and promoters
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/"$name".promoters.bed -u -f 1 > results/PCSC1/cicero/"$name"."$opname".promoters.bed
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/"$name".promoters.bed -v -f 1 > results/PCSC1/cicero/"$name"."$opname".enhancer.bed
}

filenames="GBM
PFA
LSCp"

for j in $filenames ;
do
    identify_regions $j 0 results/PCSC1/cicero/common.E-P.genes.txt "common.coaccess"
    identify_regions $j 0.2 results/PCSC1/cicero/common.E-P.genes.coaccessgt0.2.txt "common.coaccessgt0.2"
    identify_regions $j 0.5 results/PCSC1/cicero/common.E-P.genes.coaccessgt0.5.txt "common.coaccessgt0.5"

    identify_regions $j 0 results/PCSC1/cicero/shared.E-P.genes.txt "shared.coaccess"
    identify_regions $j 0.2 results/PCSC1/cicero/shared.E-P.genes.coaccessgt0.2.txt "shared.coaccessgt0.2"
    identify_regions $j 0.5 results/PCSC1/cicero/shared.E-P.genes.coaccessgt0.5.txt "shared.coaccessgt0.5"
done


#=============================================================================================
## Regions mapping to not shared genes
#=============================================================================================

## This could also include regions that are same but map to different genes

awk -F"\t" '{print > "results/PCSC1/cicero/notshared."$2".E-P.genes.txt"}' results/PCSC1/cicero/notshared.E-P.genes.txt

awk -F"\t" '{print > "results/PCSC1/cicero/notshared."$2".E-P.genes.coaccessgt0.2.txt"}' results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.2.txt

awk -F"\t" '{print > "results/PCSC1/cicero/notshared."$2".E-P.genes.coaccessgt0.5.txt"}' results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.5.txt



identify_unique_regions()
{
    name=$1
    threshold=$2
    genefile=$3
    opname=$4

    rm tmp.txt  2> /dev/null || echo > /dev/null
    rm tmp1.txt  2> /dev/null || echo > /dev/null

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } { if($3 > '${threshold}' ) print $0}' results/PCSC1/cicero/"$name".E-P.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmp.txt
    cat <(cut -d" " -f2 tmp.txt ) <(cut -d" " -f3 tmp.txt ) |  sort -k1,1 -u | cut -d"_" -f1,2,3 --output-delimiter=$'\t' > tmp1.txt

    ## identify enhancers and promoters
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/"$name".promoters.bed -u -f 1 > results/PCSC1/cicero/"$name"."$opname".promoters.bed
    intersectBed -a tmp1.txt -b results/PCSC1/cicero/"$name".promoters.bed -v -f 1 > results/PCSC1/cicero/"$name"."$opname".enhancer.bed
}

for j in $filenames ;
do
    identify_unique_regions $j 0 results/PCSC1/cicero/notshared.E-P.genes.txt "notshared.coaccess"
    identify_unique_regions $j 0.2 results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.2.txt "notshared.coaccessgt0.2"
    identify_unique_regions $j 0.5 results/PCSC1/cicero/notshared.E-P.genes.coaccessgt0.5.txt "notshared.coaccessgt0.5"
done
