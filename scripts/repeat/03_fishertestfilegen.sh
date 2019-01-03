#!/bin/bash

# For each query file, calculate jaccard index for --

#a = Jaccard index between query and repeat
#b = Jaccard index between bg-gbm and repeat
#c = Jaccard index between gbm and bgrepeat-repeat
#d = Jaccard index bg-gbm and bgrepeat-repeat

#Fisher test
#fisher.test(matrix(c(a, b, c, d), 2, 2), alternative='greater')$p.value;

module load bedtools/2.23.0
queryfile=$1
bgfile=$2
opdir=$3
qname=$4

for f in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/repeats/hg38/repname/*.bed.sorted ;
do
  name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/.bed.sorted//g')
  queryname=$(echo $queryfile |  awk -F"/" '{ print $NF}' | sed -e 's/.sorted/.bed/g')
  bedtools jaccard -a $f -b $queryfile | awk -F"\t" '{ if (NR>1) print $3, "'${name}'"}' >> "$opdir"/a."$qname".jaccard.bed
  bedtools jaccard -a $f -b data/ConsensusSet/Repeatanalysis/Bg.subtracted.$queryname | awk -F"\t" '{ if (NR>1) print $3, "'${name}'"}' >> "$opdir"/b."$qname".jaccard.bed
  bedtools jaccard -a /mnt/work1/users/lupiengroup/People/qamraa99/common.data/repeats/hg38/repname.bgsubtracted/Bg.subtracted.$name.bed  -b $queryfile | awk -F"\t" '{ if (NR>1) print $3, "'${name}'"}' >> "$opdir"/c."$qname".jaccard.bed
  bedtools jaccard -a /mnt/work1/users/lupiengroup/People/qamraa99/common.data/repeats/hg38/repname.bgsubtracted/Bg.subtracted.$name.bed  -b data/ConsensusSet/Repeatanalysis/Bg.subtracted.$queryname  | awk -F"\t" '{ if (NR>1) print $3, "'${name}'"}' >> "$opdir"/d."$qname".jaccard.bed
done
