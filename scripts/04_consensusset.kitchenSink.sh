#!/bin/bash

## Script Run in directory /mnt/work1/users/lupiengroup/People/qamraa99/HG38.HG38.Pancancer.CSC

module load R/3.4.1
module load bedtools/2.26.0  
module load crossmap/0.1.5 

DATAPATH="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/"
OPDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/"
OPNAME="KitchenSink1"
SCRIPTDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts"
COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"
ALIGNMENTSTATS="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/Alignment.Stats.txt"
GENCODETSS=$COMMONDIR"/gencode.v24.annotation.TSS.wdgeneinfo.bed"
mkdir -p $OPDIR -p $OPDIR/$OPNAME $OPDIR/$OPNAME/temp.peaks 


#-------------------------------------------------------------------------
## Consensus Sets
#-------------------------------------------------------------------------

BRAIN="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.Consensus/Brain.Consensus.sortBed"
ESC="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/Consensus.ESC.bed"
HEMATDIFF="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.differentiated.ConsensusSet.bed"
HEMATHSC="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/HSC.ConsensusSet.bed"
HEMATProg="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.Progenitor.ConsensusSet.bed"
GBM="data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak"
PFA="data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak"
HF="data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak"
LSCp="data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak"
LSCn="data/ConsensusSet/PCSC1/LSC.neg.Consensus.Catalogue.narrowPeak"
LSCb="data/ConsensusSet/PCSC1/LSC.bulk.Consensus.Catalogue.narrowPeak"
ENCODE="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/consensus.set/Consensus.bed"

#-------------------------------------------------------------------------
## Peak files
#-------------------------------------------------------------------------

gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)
lscpfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
hffiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 &&  $3=="HF" && $10==1) print $1}' $ALIGNMENTSTATS)
lscnegfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="negative" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscbulkfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="Bulk" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)



#-------------------------------------------------------------------------
## Get consensus file
#-------------------------------------------------------------------------

cat $BRAIN $ESC $HEMATDIFF $HEMATHSC $HEMATProg $GBM $PFA $HF $LSCp $LSCn $LSCb $ENCODE | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
grep -v chrM | grep -v Un | grep -v random | grep -v alt | grep -v KI27 | grep -v chrGL | grep -v chrEBV | \
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/KitchenSink1.Consensus.Catalogue.narrowPeak

awk 'BEGIN {FS=OFS="\t"} { print $1, $2-2500, $3+1000}' $GENCODETSS  | grep -v chrEBV | grep -v "chrM" | grep -v "Un" | grep -v "random" | grep -v "alt" | grep -v KI27 | grep -v chrGL | \
intersectBed -a $OPDIR/$OPNAME/KitchenSink1.Consensus.Catalogue.narrowPeak -b stdin -u > $OPDIR/$OPNAME/KitchenSink1.Consensus.Catalogue.Promoters.bed

awk 'BEGIN {FS=OFS="\t"} { print $1, $2-2500, $3+1000}' $GENCODETSS | grep -v chrEBV | grep -v "chrM" | grep -v "Un" | grep -v "random" | grep -v "alt" | grep -v KI27 | grep -v chrGL | \
intersectBed -a $OPDIR/$OPNAME/KitchenSink1.Consensus.Catalogue.narrowPeak -b stdin -v  > $OPDIR/$OPNAME/KitchenSink1.Consensus.Catalogue.Enhancers.bed

## Generate Binary matrix
cp /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.atlas/hg38/*bed $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/ESC/hg38/*.bed $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/LMPP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CMP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/MPP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CLP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/GMP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/MEP/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/HSC/*.narrowPeak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD4/*Peak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD8/*Peak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Bcell/*Peak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/NK/*Peak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Mono/*Peak $OPDIR/$OPNAME/temp.peaks
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Ery/*Peak  $OPDIR/$OPNAME/temp.peaks
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks; done
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $hffiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $lscbulkfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks; done
for f in $lscnegfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/*narrowPeak $OPDIR/$OPNAME/temp.peaks

rename .bed .narrowPeak $OPDIR/$OPNAME/temp.peaks/*.bed
rm $OPDIR/$OPNAME/temp.peaks/Consensus*

## Remove random chromosomes from the temp.peaks as well

for f in $OPDIR/$OPNAME/temp.peaks/*Peak ;
do  
    grep -v chrM $f | grep -v Un | grep -v random | grep -v alt | grep -v KI27 | grep -v chrGL | grep -v chrEBV > "$f".removed
done

rm $OPDIR/$OPNAME/temp.peaks/*Peak
rename .removed .narrowPeak $OPDIR/$OPNAME/temp.peaks/*.removed
rename .narrowPeak.narrowPeak .narrowPeak $OPDIR/$OPNAME/temp.peaks/*.narrowPeak

Rscript $SCRIPTDIR/createbinarymat.R data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.peaks ".narrowPeak" $OPDIR/$OPNAME/ KitchenSink1.Consensus.Catalogue.Binarymat.txt

