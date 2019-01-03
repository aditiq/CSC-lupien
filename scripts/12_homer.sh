#!/bin/bash

#======================================
## scABC clusters
#======================================

for f in results/PCSC1/Cluster/scABC/*Enhancer.p0.05.bed ;
do
    name=$( echo $f  | awk -F"/" '{ print $NF}' | sed -e 's/.p0.05.bed//g')
    mkdir -p results/PCSC1/Cluster/scABC/motif/${name}.bg
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $f > results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bed
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancers.bed > results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bg.bed
    qsub -q light.q -cwd -N Homer.Bg.${name} -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh \
    results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bed results/PCSC1/Cluster/scABC/motif/${name}.bg results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bg.bed
done

for f in results/PCSC1/Cluster/scABC/*Promoter.p0.05.bed ;
do
    name=$( echo $f  | awk -F"/" '{ print $NF}' | sed -e 's/.p0.05.bed//g')
    mkdir -p results/PCSC1/Cluster/scABC/motif/${name}.bg
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $f > results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bed
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoters.bed > results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bg.bed
    qsub -q light.q -cwd -N Homer.Bg.${name} -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh \
    results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bed results/PCSC1/Cluster/scABC/motif/${name}.bg results/PCSC1/Cluster/scABC/motif/${name}.bg/tmp.bg.bed
done

#======================================
## Monocle difftest regions
#======================================

for f in results/PCSC1/monocle.difftest/*beta.betapatch.qval.0.05.bed;
do
  name=$( echo $f  | awk -F"/" '{ print $NF}' | sed -e 's/.bed//g')
  name2=$(echo $name |  sed -e 's/negativebeta.betapatch.qval.0.05//g' | sed -e 's/positivebeta.betapatch.qval.0.05//g'  )

  mkdir -p results/PCSC1/monocle.difftest/homer/${name}.bg

  awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $f > results/PCSC1/monocle.difftest/homer/${name}.bg/tmp.bed
  awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' data/ConsensusSet/PCSC1/"$name2".Consensus.Catalogue.narrowPeak > results/PCSC1/monocle.difftest/homer/${name}.bg/tmp.bg.bed
  qsub -q light.q -cwd -N Homer.Bg.${name} -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh \
  results/PCSC1/monocle.difftest/homer/${name}.bg/tmp.bed results/PCSC1/monocle.difftest/homer/${name}.bg/ results/PCSC1/monocle.difftest/homer/${name}.bg/tmp.bg.bed

done

#======================================
## KCCA
#======================================

## Motif enrichment on CSC specific regions from KCCA clustering on brain tissues
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Brain.CSC.sp.txt > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg/tmp.bed
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Brain.txt  > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg/tmp.bg.bed
qsub -q light.q -cwd -N Homer.Bg.KCCA.Brain.CSC.sp -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg/tmp.bed results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg/ results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg/tmp.bg.bed

## Motif enrichment on CSC common regions from KCCA clustering not overlapping KS1 ubiquitous regions
cd temp/
intersectBed -a Common.PCSC1.bed -b K100.KitchenSink1.Flexclust.Enhancer.Promoter.Common.bed -v > Common.PCSC1.notHK.bed
mv Common.PCSC1.notHK.bed ../
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' Common.PCSC1.notHK.bed > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/common.nothk/tmp.bed
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak  > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/common.nothk/tmp.bg.bed
qsub -q light.q -cwd -N Homer.Bg.KCCA.commonPCSC1.nothk -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/common.nothk/tmp.bed results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/common.nothk/ results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/common.nothk/tmp.bg.bed

## Motif enrichment on CSC shared regions from KCCA clustering not overlapping KS1 ubiquitous regions
DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters" ;
cat $DIR/Shared.Promotergrp.txt $DIR/Shared.Enhancergrp.txt > $DIR/Shared.PCSC1.txt
intersectBed -a $DIR/Shared.PCSC1.txt -b temp/K100.KitchenSink1.Flexclust.Enhancer.Promoter.Common.bed -v > $DIR/Shared.PCSC1.notHK.bed

awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/Shared.PCSC1.notHK.bed> results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/shared.nothk/tmp.bed
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak  > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/shared.nothk/tmp.bg.bed
qsub -q light.q -cwd -N Homer.Bg.KCCA.sharedPCSC1.nothk -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/shared.nothk/tmp.bed results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/shared.nothk/ results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/shared.nothk/tmp.bg.bed

#======================================
## CREAM
#======================================

## Motif enrichment in CREAM regions
DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/CREAM/" ;

files="GBM
PFA
LSCp"

for f in $files ;
do
  awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/"$f".CREAM.bed > $DIR/homer/"$f"/tmp.bed
  cat $DIR/GBM.CREAM.bed $DIR/PFA.CREAM.bed $DIR/LSCp.CREAM.bed | sortBed -i stdin | mergeBed -i stdin | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }'  > $DIR/homer/"$f"/tmp.bg.bed
  qsub -q light.q -cwd -N Homer.CREAM."$f" -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh $DIR/homer/"$f"/tmp.bed $DIR/homer/"$f"/ $DIR/homer/"$f"/tmp.bg.bed
done

## Motif enrichment in CREAM regions
DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/CREAM/" ;

files="GBM
PFA
LSCp"

for f in $files ;
do
  awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/"$f".CREAM.bed > $DIR/homer/"$f".v2/tmp.bed
  awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' results/PCSC1/CREAM/KitchenSink1.CREAM.bed.sorted > $DIR/homer/"$f".v2/tmp.bg.bed
  qsub -q light.q -cwd -N Homer.CREAM."$f" -e stdout/homer -o stdout/homer/ ../common.scripts/callhomer.hg38.bg.sh $DIR/homer/"$f".v2/tmp.bed $DIR/homer/"$f".v2/ $DIR/homer/"$f".v2/tmp.bg.bed
done

#======================================
## CICERO
#======================================

DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/cicero/" ;

files="GBM.common.coaccess.promoters.bed
GBM.common.coaccess.enhancer.bed
GBM.common.coaccessgt0.2.promoters.bed
GBM.common.coaccessgt0.2.enhancer.bed
GBM.common.coaccessgt0.5.promoters.bed
GBM.common.coaccessgt0.5.enhancer.bed
GBM.shared.coaccess.promoters.bed
GBM.shared.coaccess.enhancer.bed
GBM.shared.coaccessgt0.2.promoters.bed
GBM.shared.coaccessgt0.2.enhancer.bed
GBM.shared.coaccessgt0.5.promoters.bed
GBM.shared.coaccessgt0.5.enhancer.bed
PFA.common.coaccess.promoters.bed
PFA.common.coaccess.enhancer.bed
PFA.common.coaccessgt0.2.promoters.bed
PFA.common.coaccessgt0.2.enhancer.bed
PFA.common.coaccessgt0.5.promoters.bed
PFA.common.coaccessgt0.5.enhancer.bed
PFA.shared.coaccess.promoters.bed
PFA.shared.coaccess.enhancer.bed
PFA.shared.coaccessgt0.2.promoters.bed
PFA.shared.coaccessgt0.2.enhancer.bed
PFA.shared.coaccessgt0.5.promoters.bed
PFA.shared.coaccessgt0.5.enhancer.bed
LSCp.common.coaccess.promoters.bed
LSCp.common.coaccess.enhancer.bed
LSCp.common.coaccessgt0.2.promoters.bed
LSCp.common.coaccessgt0.2.enhancer.bed
LSCp.common.coaccessgt0.5.promoters.bed
LSCp.common.coaccessgt0.5.enhancer.bed
LSCp.shared.coaccess.promoters.bed
LSCp.shared.coaccess.enhancer.bed
LSCp.shared.coaccessgt0.2.promoters.bed
LSCp.shared.coaccessgt0.2.enhancer.bed
LSCp.shared.coaccessgt0.5.promoters.bed
LSCp.shared.coaccessgt0.5.enhancer.bed"

bgfiles="data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak
data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak"

for f in $files ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.bed//g')

  for j in $bgfiles ;
  do
      bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.narrowPeak//g')

      mkdir $DIR/homer/$dirname.$bgname

      awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname/tmp.bg.bed
      awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/$f  > $DIR/homer/$dirname.$bgname/tmp.bed

      qsub -q highmem.q -cwd -N Homer."$dirname".$bgname  -e stdout/homer -o stdout/homer/ \
      ../common.scripts/callhomer.hg38.bg.sh \
      $DIR/homer/$dirname.$bgname/tmp.bed $DIR/homer/$dirname.$bgname/ $DIR/homer/$dirname.$bgname/tmp.bg.bed

  done
done


notsharedfiles="GBM.notshared.coaccess.enhancer.bed
LSCp.notshared.coaccessgt0.5.enhancer.bed
GBM.notshared.coaccessgt0.2.enhancer.bed
LSCp.notshared.coaccessgt0.5.promoters.bed
GBM.notshared.coaccessgt0.2.promoters.bed
LSCp.notshared.coaccess.promoters.bed
GBM.notshared.coaccessgt0.5.enhancer.bed
PFA.notshared.coaccess.enhancer.bed
GBM.notshared.coaccessgt0.5.promoters.bed
PFA.notshared.coaccessgt0.2.enhancer.bed
GBM.notshared.coaccess.promoters.bed
PFA.notshared.coaccessgt0.2.promoters.bed
LSCp.notshared.coaccess.enhancer.bed
PFA.notshared.coaccessgt0.5.enhancer.bed
LSCp.notshared.coaccessgt0.2.enhancer.bed
PFA.notshared.coaccessgt0.5.promoters.bed
LSCp.notshared.coaccessgt0.2.promoters.bed
PFA.notshared.coaccess.promoters.bed"

bgfiles="data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak
data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak"

for f in $notsharedfiles ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.bed//g')
  name2=$(echo $f | sed -e 's/notshared/shared/g')

  for j in $bgfiles ;
  do
      bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.narrowPeak//g')

      mkdir $DIR/homer/$dirname.$bgname

      awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname/tmp.bg.bed
      intersectBed -a $DIR/$f -b $DIR/$name2 -v | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' > $DIR/homer/$dirname.$bgname/tmp.bed

      qsub -q highmem.q -cwd -N Homer."$dirname".$bgname  -e stdout/homer -o stdout/homer/ \
      ../common.scripts/callhomer.hg38.bg.sh \
      $DIR/homer/$dirname.$bgname/tmp.bed $DIR/homer/$dirname.$bgname/ $DIR/homer/$dirname.$bgname/tmp.bg.bed
  done
done

for f in $notsharedfiles ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.bed//g')
  name2=$(echo $f | sed -e 's/notshared/shared/g')

  flag=`echo $f|awk '{print match($0,"enhancer")}'`;

  if [ $flag -gt 0 ];then

    for j in $bgfilesenh ;
    do
        bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.Enhancers.bed//g')
        mkdir $DIR/homer/$dirname.$bgname.run2

        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
        intersectBed -a $DIR/$f -b $DIR/$name2 -v | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' > $DIR/homer/$dirname.$bgname.run2/tmp.bed

        qsub -q light.q -cwd -N Homer."$dirname".$bgname.run2  -e stdout/homer -o stdout/homer/ \
        ../common.scripts/callhomer.hg38.bg.sh \
        $DIR/homer/$dirname.$bgname.run2/tmp.bed $DIR/homer/$dirname.$bgname.run2/ $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
    done

  else
    for j in $bgfilesprom ;
    do
        bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.Promoters.bed//g')
        mkdir $DIR/homer/$dirname.$bgname.run2

        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
        intersectBed -a $DIR/$f -b $DIR/$name2 -v | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' > $DIR/homer/$dirname.$bgname.run2/tmp.bed

        qsub -q light.q -cwd -N Homer."$dirname".$bgname.run2  -e stdout/homer -o stdout/homer/ \
        ../common.scripts/callhomer.hg38.bg.sh \
        $DIR/homer/$dirname.$bgname.run2/tmp.bed $DIR/homer/$dirname.$bgname.run2/ $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
    done
  fi
done


## run2 - Cicero

DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/cicero/run2.readcount/" ;

bgfilesenh="data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Enhancers.bed
data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Enhancers.bed
data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancers.bed"

files="GBM.common.coaccessgt0.promoters.bed
GBM.common.coaccessgt0.enhancer.bed
GBM.common.coaccessgt0.2.promoters.bed
GBM.common.coaccessgt0.2.enhancer.bed
GBM.common.coaccessgt0.5.promoters.bed
GBM.common.coaccessgt0.5.enhancer.bed
GBM.shared.coaccessgt0.promoters.bed
GBM.shared.coaccessgt0.enhancer.bed
GBM.shared.coaccessgt0.2.promoters.bed
GBM.shared.coaccessgt0.2.enhancer.bed
GBM.shared.coaccessgt0.5.promoters.bed
GBM.shared.coaccessgt0.5.enhancer.bed
PFA.common.coaccessgt0.promoters.bed
PFA.common.coaccessgt0.enhancer.bed
PFA.common.coaccessgt0.2.promoters.bed
PFA.common.coaccessgt0.2.enhancer.bed
PFA.common.coaccessgt0.5.promoters.bed
PFA.common.coaccessgt0.5.enhancer.bed
PFA.shared.coaccessgt0.promoters.bed
PFA.shared.coaccessgt0.enhancer.bed
PFA.shared.coaccessgt0.2.promoters.bed
PFA.shared.coaccessgt0.2.enhancer.bed
PFA.shared.coaccessgt0.5.promoters.bed
PFA.shared.coaccessgt0.5.enhancer.bed
LSCp.common.coaccessgt0.promoters.bed
LSCp.common.coaccessgt0.enhancer.bed
LSCp.common.coaccessgt0.2.promoters.bed
LSCp.common.coaccessgt0.2.enhancer.bed
LSCp.common.coaccessgt0.5.promoters.bed
LSCp.common.coaccessgt0.5.enhancer.bed
LSCp.shared.coaccessgt0.promoters.bed
LSCp.shared.coaccessgt0.enhancer.bed
LSCp.shared.coaccessgt0.2.promoters.bed
LSCp.shared.coaccessgt0.2.enhancer.bed
LSCp.shared.coaccessgt0.5.promoters.bed
LSCp.shared.coaccessgt0.5.enhancer.bed"

for f in $files ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.bed//g')
  flag=`echo $f | awk '{print match($0,"enhancer")}'`;
  if [ $flag -gt 0 ];
  then

    for j in $bgfilesenh ;
    do
        bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.Enhancers.bed//g')
        mkdir $DIR/homer/$dirname.$bgname.run2

        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/$f  > $DIR/homer/$dirname.$bgname.run2/tmp.bed

        qsub -q lupiengroup -cwd -N Homer."$dirname".$bgname.run2  -e stdout/homer/cicero/ -o stdout/homer/cicero/ \
        ../common.scripts/callhomer.hg38.bg.sh \
        $DIR/homer/$dirname.$bgname.run2/tmp.bed $DIR/homer/$dirname.$bgname.run2/ $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
    done
  fi
done


notsharedfiles="GBM.notshared.coaccessgt0.enhancer.bed
LSCp.notshared.coaccessgt0.5.enhancer.bed
GBM.notshared.coaccessgt0.2.enhancer.bed
LSCp.notshared.coaccessgt0.5.promoters.bed
GBM.notshared.coaccessgt0.2.promoters.bed
LSCp.notshared.coaccessgt0.promoters.bed
GBM.notshared.coaccessgt0.5.enhancer.bed
PFA.notshared.coaccessgt0.enhancer.bed
GBM.notshared.coaccessgt0.5.promoters.bed
PFA.notshared.coaccessgt0.2.enhancer.bed
GBM.notshared.coaccessgt0.promoters.bed
PFA.notshared.coaccessgt0.2.promoters.bed
LSCp.notshared.coaccessgt0.enhancer.bed
PFA.notshared.coaccessgt0.5.enhancer.bed
LSCp.notshared.coaccessgt0.2.enhancer.bed
PFA.notshared.coaccessgt0.5.promoters.bed
LSCp.notshared.coaccessgt0.2.promoters.bed
PFA.notshared.coaccessgt0.promoters.bed"

for f in $notsharedfiles ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.bed//g')
  name2=$(echo $f | sed -e 's/notshared/shared/g')

  flag=`echo $f|awk '{print match($0,"enhancer")}'`;

  if [ $flag -gt 0 ];then

    for j in $bgfilesenh ;
    do
        bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.Enhancers.bed//g')
        mkdir $DIR/homer/$dirname.$bgname.run2

        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
        intersectBed -a $DIR/$f -b $DIR/$name2 -v | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' > $DIR/homer/$dirname.$bgname.run2/tmp.bed

        qsub -q light.q -cwd -N Homer."$dirname".$bgname.run2  -e stdout/homer/cicero/ -o stdout/homer/cicero/ \
        ../common.scripts/callhomer.hg38.bg.sh \
        $DIR/homer/$dirname.$bgname.run2/tmp.bed $DIR/homer/$dirname.$bgname.run2/ $DIR/homer/$dirname.$bgname.run2/tmp.bg.bed
    done
  fi
done
