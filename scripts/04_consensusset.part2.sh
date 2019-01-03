#!/bin/bash
#-------------------------------------------------------------------------
## Script Run in directory /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC
#-------------------------------------------------------------------------

module load R/3.4.1
module load bedtools/2.26.0
module load crossmap/0.1.5

DATAPATH="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/"
OPDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/"
OPNAME="PCSC1"
SCRIPTDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts"
COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"
ALIGNMENTSTATS="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/Alignment.Stats.txt"
GENCODETSS=$COMMONDIR"/gencode.v24.annotation.TSS.wdgeneinfo.bed"
ATACDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue"
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

gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)
lscpfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscnegfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="negative" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscbulkfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="Bulk" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
hffiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 &&  $3=="HF" && $10==1) print $1}' $ALIGNMENTSTATS)

#-------------------------------------------------------------------------
## Create consensus
#-------------------------------------------------------------------------

## LSCp.ESC
cat $LSCp $ESC  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/LSCp.ESC.Consensus.Catalogue.narrowPeak

## LSCp.HSC
cat $LSCp $HEMATHSC  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/LSCp.HEMATHSC.Consensus.Catalogue.narrowPeak

## LSCp.HEMATDIFF
cat $LSCp $HEMATDIFF  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/LSCp.HEMATDIFF.Consensus.Catalogue.narrowPeak

## LSCp.LSCbulk
cat $LSCp $LSCb  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/LSCp.LSCb.Consensus.Catalogue.narrowPeak

## GBM.ESC
cat $GBM $ESC  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/GBM.ESC.Consensus.Catalogue.narrowPeak

## PFA.ESC
cat $PFA $ESC  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/PFA.ESC.Consensus.Catalogue.narrowPeak

## GBM.BrainAtlas
cat $GBM $BRAIN $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF682RLN_peaks.narrowPeak $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF114WRE_peaks.narrowPeak | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/GBM.BRAIN.Consensus.Catalogue.narrowPeak

## PFA.BrainAtlas
cat $PFA $BRAIN $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF682RLN_peaks.narrowPeak $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF114WRE_peaks.narrowPeak | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/PFA.BRAIN.Consensus.Catalogue.narrowPeak

## CSC.ESC
cat $PFA $GBM $LSCp $ESC  | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' |\
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/CSC.ESC.Consensus.Catalogue.narrowPeak

#-------------------------------------------------------------------------
## Copy peaks to temp folder for binary matrix construction
#-------------------------------------------------------------------------

mkdir -p $OPDIR/$OPNAME/temp.LSCp.ESC
mkdir -p $OPDIR/$OPNAME/temp.LSCp.HSC
mkdir -p $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF
mkdir -p $OPDIR/$OPNAME/temp.LSCp.LSCb
mkdir -p $OPDIR/$OPNAME/temp.GBM.ESC
mkdir -p $OPDIR/$OPNAME/temp.PFA.ESC
mkdir -p $OPDIR/$OPNAME/temp.GBM.BRAIN
mkdir -p $OPDIR/$OPNAME/temp.PFA.BRAIN
mkdir -p $OPDIR/$OPNAME/temp.CSC.ESC

## LSCp.ESC
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.ESC/ ; done
cp $ATACDIR/ESC/hg38/*.bed $OPDIR/$OPNAME/temp.LSCp.ESC/

## LSCp.HSC
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.HSC/ ; done
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/HSC/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.HSC/

## LSCp.HEMATDIFF
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/ ; done
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD4/*Peak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD8/*Peak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Bcell/*Peak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/NK/*Peak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Mono/*Peak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Ery/*Peak  $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/

## LSCp.LSCbulk
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.LSCb/ ; done
for f in $lscbulkfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.LSCp.LSCb/ ; done

## GBM.ESC
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.GBM.ESC/ ; done
cp $ATACDIR/ESC/hg38/*.bed $OPDIR/$OPNAME/temp.GBM.ESC/

## PFA.ESC
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.PFA.ESC/ ; done
cp $ATACDIR/ESC/hg38/*.bed $OPDIR/$OPNAME/temp.PFA.ESC/

## GBM.BrainAtlas
cp /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.atlas/hg38/*bed $OPDIR/$OPNAME/temp.GBM.BRAIN/
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.GBM.BRAIN/ ; done
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF682RLN_peaks.narrowPeak  $OPDIR/$OPNAME/temp.GBM.BRAIN/
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF114WRE_peaks.narrowPeak $OPDIR/$OPNAME/temp.GBM.BRAIN/

## PFA.BrainAtlas
cp /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.atlas/hg38/*bed $OPDIR/$OPNAME/temp.PFA.BRAIN/
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.PFA.BRAIN/ ; done
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF682RLN_peaks.narrowPeak  $OPDIR/$OPNAME/temp.PFA.BRAIN/
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF114WRE_peaks.narrowPeak $OPDIR/$OPNAME/temp.PFA.BRAIN/

## CSCs and ESCs
cp $ATACDIR/ESC/hg38/*.bed $OPDIR/$OPNAME/temp.CSC.ESC/
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.CSC.ESC/ ; done
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.CSC.ESC/ ; done
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.CSC.ESC/ ; done

## Rename all .bed file to .narrowPeak in temp folders
rename .bed .narrowPeak $OPDIR/$OPNAME/temp*/*.bed

## Generate binary matrix
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCp.ESC.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.LSCp.ESC/ ".narrowPeak" $OPDIR/$OPNAME/ LSCp.ESC.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCp.HEMATHSC.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.LSCp.HSC/ ".narrowPeak" $OPDIR/$OPNAME/ LSCp.HEMATHSC.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCp.LSCb.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.LSCp.LSCb/ ".narrowPeak" $OPDIR/$OPNAME/ LSCp.LSCb.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/GBM.ESC.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.GBM.ESC/ ".narrowPeak" $OPDIR/$OPNAME/ GBM.ESC.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/PFA.ESC.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.PFA.ESC/ ".narrowPeak" $OPDIR/$OPNAME/ PFA.ESC.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/GBM.BRAIN.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.GBM.BRAIN/ ".narrowPeak" $OPDIR/$OPNAME/ GBM.BRAIN.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/PFA.BRAIN.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.PFA.BRAIN/ ".narrowPeak" $OPDIR/$OPNAME/ PFA.BRAIN.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/CSC.ESC.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.CSC.ESC/ ".narrowPeak" $OPDIR/$OPNAME/ CSC.ESC.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCp.HEMATDIFF.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.LSCp.HEMATDIFF/ ".narrowPeak" $OPDIR/$OPNAME/ LSCp.HEMATDIFF.Consensus.Catalogue.Binarymat.txt


#-------------------------------------------------------------------------
### 9th August 2018
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

gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)
lscpfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscnegfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="negative" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscbulkfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="Bulk" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
hffiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 &&  $3=="HF" && $10==1) print $1}' $ALIGNMENTSTATS)


mkdir -p $OPDIR/$OPNAME/temp.blood
mkdir -p $OPDIR/$OPNAME/temp.brain


## All blood related data
for f in $lscpfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.blood/ ; done
for f in $lscnegfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.blood/ ; done
for f in $lscbulkfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.blood/ ; done
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/HSC/*Peak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/LMPP/*Peak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CMP/*Peak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/MPP/*Peak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CLP/*Peak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/GMP/*Peak  $OPDIR/$OPNAME/temp.LSCp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/MEP/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD4/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/CD8/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Bcell/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/NK/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Mono/*.narrowPeak $OPDIR/$OPNAME/temp.blood/
cp $ATACDIR/RyanCorces.Hemat/hg38/peakCalls/Ery/*.narrowPeak $OPDIR/$OPNAME/temp.blood/


cat $OPDIR/$OPNAME/temp.blood/* | \
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' | grep -v "chrM" | grep -v "chrEBV" | grep -v "random" |  grep -v "chrUn" |  grep -v "alt" | \
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/Combined.Blood.Consensus.Catalogue.narrowPeak

## All brain related data
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.brain/ ; done
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.brain/ ; done
for f in $hffiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.brain/ ; done
cp /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.atlas/hg38/*bed $OPDIR/$OPNAME/temp.brain/
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF682RLN_peaks.narrowPeak  $OPDIR/$OPNAME/temp.brain/
cp $ATACDIR/ATAC-ENCODE/hg38/peaks/ENCFF114WRE_peaks.narrowPeak $OPDIR/$OPNAME/temp.brain/

cat $OPDIR/$OPNAME/temp.brain/* | \
awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' | grep -v "chrM" | grep -v "chrEBV" | grep -v "random" |  grep -v "chrUn" |  grep -v "alt" | \
sortBed -i stdin | mergeBed -i stdin | sortBed -i stdin > $OPDIR/$OPNAME/Combined.Brain.Consensus.Catalogue.narrowPeak

for f in $OPDIR/$OPNAME/temp.brain/* ;
do
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' $f | grep -v "chrM" | grep -v "chrEBV" | grep -v "random" |  grep -v "chrUn" |  grep -v "alt" > "$f".v2 ;
done

for f in $OPDIR/$OPNAME/temp.blood/* ;
do
    awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' $f | grep -v "chrM" | grep -v "chrEBV" | grep -v "random" |  grep -v "chrUn" |  grep -v "alt" > "$f".v2 ;
done


Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/Combined.Brain.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.brain/ ".narrowPeak.v2" $OPDIR/$OPNAME/ Combined.Brain.Consensus.Catalogue.Binarymat.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/Combined.Blood.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.blood/ ".narrowPeak.v2" $OPDIR/$OPNAME/ Combined.Blood.Consensus.Catalogue.Binarymat.txt

## Binary catalogue with remap TFs annotation
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/Combined.Blood.Consensus.Catalogue.narrowPeak /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/ ".bed" $OPDIR/$OPNAME/ Combined.Blood.Consensus.Catalogue.Binarymat.with.remap.txt
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/Combined.Brain.Consensus.Catalogue.narrowPeak /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/ ".bed" $OPDIR/$OPNAME/ Combined.Brain.Consensus.Catalogue.Binarymat.with.remap.txt

## Binary catalogue for shared PCSC1 CSC regions
DIR="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters" ;
Rscript $SCRIPTDIR/createbinarymat.R $DIR/Shared.PCSC1.notHK.bed data/ConsensusSet/KitchenSink2/temp.peaks/ ".narrowPeak" ./ Shared.PCSC1.notHK.KS2.Binarymat.txt
