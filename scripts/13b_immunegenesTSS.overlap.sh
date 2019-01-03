#!/bin/bash

#-----------------------------------------------------------
## Objective : Intersect atac peaks with immune genes TSS
#-----------------------------------------------------------

module load bedtools/2.26.0


## In house ATAC consensus set
intersectBed -a data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak.bed -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/PCSC1ConsensusCatalogue.ImmuneTSS.intersectu.bed

## GBM
intersectBed -a data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/GBM.ConsensusCatalogue.ImmuneTSS.intersectu.bed

## PFA
intersectBed -a data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/PFA.ConsensusCatalogue.ImmuneTSS.intersectu.bed

## LSCp
intersectBed -a data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/LSCp.ConsensusCatalogue.ImmuneTSS.intersectu.bed

## LSC-negative
intersectBed -a data/ConsensusSet/PCSC1/LSC.neg.Consensus.Catalogue.narrowPeak -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/LSCneg.ConsensusCatalogue.ImmuneTSS.intersectu.bed

## HF
intersectBed -a data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/HF.ConsensusCatalogue.ImmuneTSS.intersectu.bed

## ATAC from https://www.encodeproject.org/matrix/?type=Experiment&assay_title=ATAC-seq
for f in  /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-ENCODE/peaks/ENCFF*Peak ;
do
    name=$(echo $f | awk -F"/" '{print $NF}'  | sed -e 's/_peaks.narrowPeak//g' )
    intersectBed -a $f -b data/immunelandscape/TSSofImmuneGenes.bed -wo > data/immunelandscape/$name.ImmuneTSS.intersectu.bed
done

