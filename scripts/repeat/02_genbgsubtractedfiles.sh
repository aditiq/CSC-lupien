#!/bin/bash


queryfiles="data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/LSC.neg.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/LSC.bulk.Consensus.Catalogue.narrowPeak
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/HSC.ConsensusSet.bed
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.Progenitor.ConsensusSet.bed
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.differentiated.ConsensusSet.bed
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/Consensus.ESC.bed
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.Consensus/Brain.Consensus.sortBed"

bgfile="data/ConsensusSet/Repeatanalysis/Background.merged.bed"

for f in $queryfiles ;
do
  name=$(echo $f | awk -F"/" '{ print $NF}' )
  cmd1="module load bedtools/2.23.0 ;
  intersectBed -a \"$bgfile\" -b \"$f\" -v | sortBed -i stdin > data/ConsensusSet/Repeatanalysis/Bg.subtracted.\"$name\".bed"

  echo $cmd1 > scripts/repeat/"$name".genbgsubtracted.sh
  sed -i 's/"//g' scripts/repeat/"$name".genbgsubtracted.sh
  qsub -q lupiengroup -cwd -e stdout/repeat/ -o stdout/repeat/ -N "$name".bgsubtraction scripts/repeat/"$name".genbgsubtracted.sh

done
