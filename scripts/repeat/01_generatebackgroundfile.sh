#!/bin/bash

## Comine all ATAC-seq peaks to generate background for repeat enrichment analysis

module load bedtools/2.23.0

cat /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/peaks/*/*.narrowPeak \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/peaks/*.narrowPeak \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.atlas/hg38/*-includes-blacklisted-regions_peaks.hg38.bed \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/*bed >> data/ConsensusSet/Repeatanalysis/Background.bed

awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}'  data/ConsensusSet/Repeatanalysis/Background.bed  | sortBed -i stdin | mergeBed -i stdin > data/ConsensusSet/Repeatanalysis/Background.merged.bed
