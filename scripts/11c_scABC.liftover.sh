#!/bin/bash


## Convert scABC peaks to hg19 for use with GREAT analysis

COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"

for f in results/PCSC1/Cluster/scABC/*p0.05.bed ;
do 
     liftOver $f $COMMONDIR/hg38ToHg19.over.chain.gz "$f".hg19 "$f".hg19.unmap ; 
done


names="scABC.Combined.LSCnegvsLSCpos.LSCneg.p0.05
scABC.Combined.LSCnegvsLSCpos.LSCpos.p0.05"

for f in $names ;
do
  liftOver "$f".bed /mnt/work1/users/lupiengroup/People/qamraa99/common.data/hg38ToHg19.over.chain "$f".hg19.bed "$f".unmap
done

liftOver LSCposandneg.Consensus.Catalogue.narrowPeak /mnt/work1/users/lupiengroup/People/qamraa99/common.data/hg38ToHg19.over.chain LSCposandneg.Consensus.Catalogue.narrowPeak.hg19 LSCposandneg.Consensus.Catalogue.narrowPeak.unmap 


