#!/bin/bash

## Script Run in directory /mnt/work1/users/lupiengroup/People/qamraa99/HG38.HG38.Pancancer.CSC

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
mkdir -p $OPDIR -p $OPDIR/$OPNAME $OPDIR/$OPNAME/temp.peaks 

##########################################################################################################################
## Generate catalogues
##########################################################################################################################

gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)
lscfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)

## Concatenate
for f in $gbmfiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/GBM.Catalogue.narrowPeak ; done
for f in $pfafiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/PFA.Catalogue.narrowPeak ; done
for f in $lscfiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/LSCp.Catalogue.narrowPeak ; done

cut -d"_" -f1 $OPDIR/$OPNAME/PFA.Catalogue.narrowPeak | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/PFA.Consensus.Catalogue.narrowPeak
cut -d"_" -f1 $OPDIR/$OPNAME/GBM.Catalogue.narrowPeak | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/GBM.Consensus.Catalogue.narrowPeak
sed -e 's/_peak/;peak/g' $OPDIR/$OPNAME/LSCp.Catalogue.narrowPeak | cut -d";" -f1  | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/LSCp.Consensus.Catalogue.narrowPeak

##################################
## Combine and merge peaks
##################################

cat $OPDIR/$OPNAME/GBM.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/PFA.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/LSCp.Consensus.Catalogue.narrowPeak | \
sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct | sortBed -i stdin > $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.narrowPeak


##########################################################################################################################
## Split into Promoters and Enhancers 
##########################################################################################################################

## Split into enhancers and promoters

awk 'BEGIN {FS=OFS="\t"} { print $1, $2-2500, $3+1000}' $GENCODETSS  | grep -v "chrM" |\
intersectBed -a $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.narrowPeak -b stdin -u > $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.Promoters.bed

awk 'BEGIN {FS=OFS="\t"} { print $1, $2-2500, $3+1000}' $GENCODETSS | grep -v "chrM" |\
intersectBed -a $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.narrowPeak -b stdin -v > $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.Enhancers.bed

##########################################################################################################################
## Copy peaks to temp folder for binary matrix construction
##########################################################################################################################
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $lscfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done

Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/$OPNAME.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.peaks ".narrowPeak" $OPDIR/$OPNAME/ $OPNAME.Consensus.Catalogue.Binarymat.txt
rm $OPDIR/$OPNAME/temp.peaks/*


##########################################################################################################################
## Create consensus for HF, LSC Bulk and LSC negative
##########################################################################################################################
hffiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 &&  $3=="HF" && $10==1) print $1}' $ALIGNMENTSTATS)
lscnegfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="negative" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscbulkfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="Bulk" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)

## Concatenate
for f in $hffiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/HF.Catalogue.narrowPeak ; done
for f in $lscnegfiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/LSC.neg.Catalogue.narrowPeak ; done
for f in $lscbulkfiles ; do cat $DATAPATH/peaks/$f/*.narrowPeak >> $OPDIR/$OPNAME/LSC.bulk.Catalogue.narrowPeak ; done

cut -d"_" -f1 $OPDIR/$OPNAME/HF.Catalogue.narrowPeak | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/HF.Consensus.Catalogue.narrowPeak
sed -e 's/_peak/;peak/g' $OPDIR/$OPNAME/LSC.neg.Catalogue.narrowPeak | cut -d";" -f1  | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/LSC.neg.Consensus.Catalogue.narrowPeak
sed -e 's/_peak/;peak/g' $OPDIR/$OPNAME/LSC.bulk.Catalogue.narrowPeak | cut -d";" -f1  | sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct > $OPDIR/$OPNAME/LSC.bulk.Consensus.Catalogue.narrowPeak

## gbm, pfa, hf

cat $OPDIR/$OPNAME/GBM.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/PFA.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/HF.Consensus.Catalogue.narrowPeak | \
sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct | sortBed -i stdin > $OPDIR/$OPNAME/GBM.PFA.HF.Consensus.Catalogue.narrowPeak


## gbm, hf
cat $OPDIR/$OPNAME/GBM.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/HF.Consensus.Catalogue.narrowPeak | \
sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct | sortBed -i stdin > $OPDIR/$OPNAME/GBM.HF.Consensus.Catalogue.narrowPeak

##lsc positive and lsc negative
cat $OPDIR/$OPNAME/LSCp.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/LSC.neg.Consensus.Catalogue.narrowPeak  | \
sortBed -i stdin | mergeBed -i stdin -c 4 -o count_distinct | sortBed -i stdin > $OPDIR/$OPNAME/LSCposandneg.Consensus.Catalogue.narrowPeak


##########################################################################################################################
## Copy peaks to temp folder for binary matrix construction
##########################################################################################################################

mkdir -p $OPDIR/$OPNAME/temp.gbm.pfa.hf
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm.pfa.hf/ ; done
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm.pfa.hf/ ; done
for f in $hffiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm.pfa.hf/ ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/GBM.PFA.HF.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.gbm.pfa.hf ".narrowPeak" $OPDIR/$OPNAME/ GBM.PFA.HF.Consensus.Catalogue.Binarymat.txt

mkdir -p $OPDIR/$OPNAME/temp.gbm.hf
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm.hf/ ; done
for f in $hffiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm.hf/ ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/GBM.HF.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.gbm.hf ".narrowPeak" $OPDIR/$OPNAME/ GBM.HF.Consensus.Catalogue.Binarymat.txt

mkdir -p $OPDIR/$OPNAME/temp.lscp.lscneg/ 
for f in $lscnegfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.lscp.lscneg/  ; done
for f in $lscfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.lscp.lscneg/  ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCposandneg.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.lscp.lscneg/  ".narrowPeak" $OPDIR/$OPNAME/ LSCposandneg.Consensus.Catalogue.Binarymat.txt


mkdir -p $OPDIR/$OPNAME/temp.gbm
for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.gbm/  ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/GBM.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.gbm/  ".narrowPeak" $OPDIR/$OPNAME/ GBM.Consensus.Catalogue.Binarymat.txt


mkdir -p $OPDIR/$OPNAME/temp.pfa
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.pfa/  ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/PFA.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.pfa/  ".narrowPeak" $OPDIR/$OPNAME/ PFA.Consensus.Catalogue.Binarymat.txt


mkdir -p $OPDIR/$OPNAME/temp.lscp
for f in $lscfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.lscp  ; done
Rscript $SCRIPTDIR/createbinarymat.R $OPDIR/$OPNAME/LSCp.Consensus.Catalogue.narrowPeak $OPDIR/$OPNAME/temp.lscp  ".narrowPeak" $OPDIR/$OPNAME/ LSCp.Consensus.Catalogue.Binarymat.txt
