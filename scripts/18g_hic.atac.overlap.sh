#!/bin/bash

# objective : identify  which HiC loops lie in ATAC regions and  which cicero loops lie in HiC

## LSCp with OCIAML2 HiC

loopfile="data/HiC/LSC/All_Chr_Merged_Loops.bedpe"
atacregions="data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak"
cicerofile="results/PCSC1/cicero/LSCp.deduped.connsfull.readcount.pairs.bed"

## identify which HiC loops lie in ATAC regions
awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$1,$2,$3, "chr"$4,$5,$6}' $loopfile | intersectBed -a stdin -b $atacregions -u >> data/HiC/LSC/All_Chr_Merged_Loops.LSCatac.bedpe
awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$4,$5,$6, "chr"$1,$2,$3}' $loopfile | intersectBed -a stdin -b $atacregions -u | \
awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print $4,$5,$6, $1,$2,$3}'   >> data/HiC/LSC/All_Chr_Merged_Loops.LSCatac.bedpe

sort -u data/HiC/LSC/All_Chr_Merged_Loops.LSCatac.bedpe > data/HiC/LSC/All_Chr_Merged_Loops.LSCatac.dedup.bedpe

## Which atac regions lie within 1kb of a Hic Loops
awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$1,$2,$3}' $loopfile | windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/LSCp.regions.in.HiC.bed
awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$4,$5,$6}' $loopfile | windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/LSCp.regions.in.HiC.bed


## GBM with GSC HiC

loopfile="data/HiC/GBM/merged_loops_for_G523_aggr_ChIP_with_motifs.bedpe
data/HiC/GBM/merged_loops_for_G583_aggr_ChIP_with_motifs.bedpe
data/HiC/GBM/merged_loops_for_G567_aggr_ChIP_with_motifs.bedpe"

atacregions="data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak"
cicerofile="results/PCSC1/cicero/run2.readcount/GBM.deduped.connsfull.pairs.bed"

## Which atac regions lie within 1kb of a Hic Loops
for f in $loopfile ;
do
  name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/merged_loops_for_//g' | sed -e 's/_aggr_ChIP_with_motifs.bedpe//g')

  awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$1,$2,$3}' $f | \
  windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/GBM.regions.in.$name.HiC.bed

  awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$4,$5,$6}' $f | \
  windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/GBM.regions.in.$name.HiC.bed
done



cat data/HiC/GBM/merged_loops_for_G523_aggr_ChIP_with_motifs.bedpe \
data/HiC/GBM/merged_loops_for_G583_aggr_ChIP_with_motifs.bedpe \
data/HiC/GBM/merged_loops_for_G567_aggr_ChIP_with_motifs.bedpe > tmp.txt

awk 'BEGIN {FS=OFS="\t"} {  print "chr"$1,$2,$3}' tmp.txt | grep -v x1 | \
windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' > results/PCSC1/cicero/GBM.regions.in.HiC.bed

awk 'BEGIN {FS=OFS="\t"} {  print "chr"$4,$5,$6}' tmp.txt | grep -v y1 | \
windowBed -a $atacregions -b stdin -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/GBM.regions.in.HiC.bed

## Which HiC regions overlap ATAC peaks
for f in $loopfile ;
do
  name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/merged_loops_for_//g' | sed -e 's/_aggr_ChIP_with_motifs.bedpe//g')

  awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$1,$2,$3}' $f | \
  windowBed -a stdin -b data/peaks/$name/*Peak -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/$name.HiC.regions.in.$name.atac.bed

  awk 'BEGIN {FS=OFS="\t"} {  if (NR>1) print "chr"$4,$5,$6}' $f | \
  windowBed -a stdin -b data/peaks/$name/*Peak -w 1000 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' >> results/PCSC1/cicero/$name.HiC.regions.in.$name.atac.bed
done
