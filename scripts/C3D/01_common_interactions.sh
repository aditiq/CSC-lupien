module load bedtools/2.23.0

# Identify no. of E-P interactions mapping to common genes and different genes

for f in results/PCSC1/C3D/finalresults/*combined.results.txt ;
do
  name=$(echo $f | awk -F"/" '{ print $NF}'  | sed -e 's/.combined.results.txt//g')
  awk 'BEGIN  {FS=OFS="\t"} { if($NF <0.01 && $7 >= 0.7 && $1!="COORD_1") \
  print $1,$4,$7,$8,$9,$10}' $f > results/PCSC1/C3D/finalresults/$name.sig.results.txt
done

## identify common genes present in all three or any two of the populations
names="GBM
PFA
LSCp" ;

for f in $names ;
do
    awk '{print $4, "'${f}'"}' results/PCSC1/C3D/finalresults/$f.sig.results.txt | \
    awk '!seen[$0]++'  > \
    "$f".c3d.tmp.txt
done

cat *.c3d.tmp.txt | \
awk -F" " '{a[$1]=a[$1]?a[$1]":"$2:$2;}END{for (i in a)print i, a[i];}' OFS="\t" | \
awk -F ":" ' { print NF-1"\t"$0 } '  > results/PCSC1/C3D/finalresults/SigResults.genefreq.txt

awk 'BEGIN {FS=OFS="\t"} { if ($1>1)  print $2,$3}'  results/PCSC1/C3D/finalresults/SigResults.genefreq.txt > results/PCSC1/C3D/finalresults/common.interactions.genes.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1>0)  print $2,$3}'  results/PCSC1/C3D/finalresults/SigResults.genefreq.txt > results/PCSC1/C3D/finalresults/shared.interactions.genes.txt
awk 'BEGIN {FS=OFS="\t"} { if ($1==0)  print $2,$3}' results/PCSC1/C3D/finalresults/SigResults.genefreq.txt > results/PCSC1/C3D/finalresults/notshared.interactions.genes.txt

#============================================================================================================================================================
## Identify regions from GBM, PFA and LSC for common, shared and unique regions
#============================================================================================================================================================

#=============================================================================================
## Common genes and Shared genes
#=============================================================================================

identify_regions()
{
    name=$1
    genefile=$2
    opname=$3
    anchor="../common.data/c3d.anchors.gencodev24.500bparoundTSS.bed"

    rm tmp.txt  2> /dev/null || echo > /dev/null
    rm tmp1.txt  2> /dev/null || echo > /dev/null

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } { print $0}' \
    results/PCSC1/C3D/finalresults/"$name".sig.results.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmp.txt

    cat <(cut -d" " -f2 tmp.txt ) <(cut -d" " -f3 tmp.txt ) | \
    sort -k1,1 -u | cut -d":" -f1,2 --output-delimiter=$'\t' | \
    cut -d"-" -f1,2 --output-delimiter=$'\t'  | \
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3}' > tmp1.txt

    ## split into TSS and ATAC Regions
    intersectBed -a tmp1.txt -b $anchor -u -f 1 > results/PCSC1/C3D/finalresults/"$name"."$opname".anchors.bed
    intersectBed -a tmp1.txt -b $anchor -v -f 1 > results/PCSC1/C3D/finalresults/"$name"."$opname".regions.bed
}



filenames="GBM
PFA
LSCp"

for j in $filenames ;
do
    identify_regions $j results/PCSC1/C3D/finalresults/common.interactions.genes.txt "common.coaccess"
    identify_regions $j results/PCSC1/C3D/finalresults/shared.interactions.genes.txt "shared.coaccess"
done


#=============================================================================================
## Regions mapping to not shared genes
#=============================================================================================

## This could also include regions that are same but map to different genes

awk -F"\t" '{print > "results/PCSC1/C3D/finalresults/notshared."$2".interactions.genes.txt"}' results/PCSC1/C3D/finalresults/notshared.interactions.genes.txt

identify_unique_regions_ns()
{
    name=$1
    genefile=$2
    opname=$3
    anchor="../common.data/c3d.anchors.gencodev24.500bparoundTSS.bed"

    rm tmp.txt  2> /dev/null || echo > /dev/null
    rm tmp1.txt  2> /dev/null || echo > /dev/null

    join -1 4 -2 1 <(awk ' BEGIN {FS=OFS="\t" } {  print $0}' results/PCSC1/C3D/finalresults/"$name".sig.results.txt | \
    sort -k4,4 ) <( cut -f1 $genefile | sort -k1,1 -u)  > tmp.txt
    cat <(cut -d" " -f2 tmp.txt ) <(cut -d" " -f3 tmp.txt ) |  \
    sort -k1,1 -u | cut -d":" -f1,2 --output-delimiter=$'\t' | \
    cut -d"-" -f1,2 --output-delimiter=$'\t'  | \
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3}' > tmp1.txt

    ## split into TSS and ATAC Regions
    intersectBed -a tmp1.txt -b $anchor -u -f 1 > results/PCSC1/C3D/finalresults/"$name"."$opname".anchors.bed
    intersectBed -a tmp1.txt -b $anchor -v -f 1 > results/PCSC1/C3D/finalresults/"$name"."$opname".regions.bed

}

for j in $filenames ;
do
    identify_unique_regions_ns $j results/PCSC1/C3D/finalresults/notshared.interactions.genes.txt "notshared.coaccess"
done
