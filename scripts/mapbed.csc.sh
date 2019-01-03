#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Get signal matrix for GBM , LSC and PFA separately
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


module load R/3.4.1
module load bedtools/2.26.0

DATAPATH="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/"
OPDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/"
OPNAME="PCSC1"
SCRIPTDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts"
COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"
ALIGNMENTSTATS="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/Alignment.Stats.txt"

lscfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)

gbmcons="data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak.sorted"
lscpcons="data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak.sorted"
pfacons="data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak.sorted"
pcsc1cons="data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak"

for f in $lscfiles ;
do
    qsub -cwd -N Mapbed.LSCp.Consensus.$f -e stdout/mapbed/lscp -o stdout/mapbed/lscp  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $lscpcons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/lscp/ LSCp."$f".Consensus.mapbed
done

for f in $gbmfiles ;
do
    qsub -cwd -N Mapbed.GBM.Consensus.$f -e stdout/mapbed/gbm -o stdout/mapbed/gbm  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $gbmcons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/gbm/ GBM."$f".Consensus.mapbed
done

for f in $pfafiles ;
do
    qsub -cwd -N Mapbed.PFA.Consensus.$f -e stdout/mapbed/pfa -o stdout/mapbed/pfa  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $pfacons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/pfa/ PFA."$f".Consensus.mapbed
done


for f in $lscfiles ;
do
    qsub -cwd -N Mapbed.PCSC1.Consensus.$f -e stdout/mapbed/pcsc1 -o stdout/mapbed/pcsc1  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $pcsc1cons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/pcsc1/ PCSC1."$f".Consensus.mapbed
done

for f in $gbmfiles ;
do
    qsub -cwd -N Mapbed.PCSC1.Consensus.$f -e stdout/mapbed/pcsc1 -o stdout/mapbed/pcsc1  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $pcsc1cons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/pcsc1/ PCSC1."$f".Consensus.mapbed
done

for f in $pfafiles ;
do
    qsub -cwd -N Mapbed.PCSC1.Consensus.$f -e stdout/mapbed/pcsc1 -o stdout/mapbed/pcsc1  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $pcsc1cons $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/pcsc1/ PCSC1."$f".Consensus.mapbed
done
