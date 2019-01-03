#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify LSCp and LSCn samples that don't behave similar to their group i.e LSCp that show ATAC signal like LSCn and vice versa
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

module load R/3.4.1
module load bedtools/2.26.0
module load crossmap/0.1.5

DATAPATH="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/"
OPDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/"
OPNAME="PCSC1"
SCRIPTDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.scripts"
COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"
ALIGNMENTSTATS="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/Alignment.Stats.txt"


lscfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)
lscnegfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="negative" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)

#--------------------------------------------------------------------------------------------------------------------------------
## Get Coverage for LSC pos and neg consensus set and run diagnostics
#--------------------------------------------------------------------------------------------------------------------------------


for f in $lscfiles ;
do
    qsub -cwd -N Mapbed.All.LSC.Consensus.$f -e stdout/mapbed/lscn.lscp/ -o stdout/mapbed/lscn.lscp/  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $OPDIR/$OPNAME/LSCposandneg.Consensus.Catalogue.narrowPeak $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/LSCn.LSCp LSCp."$f".All.LSC.Consensus.mapbed
done

for f in $lscnegfiles ;
do
    qsub -cwd -N Mapbed.All.LSC.Consensus.$f -e stdout/mapbed/lscn.lscp/ -o stdout/mapbed/lscn.lscp/  -q light.q $SCRIPTDIR/generic_mapbed.sh \
    $OPDIR/$OPNAME/LSCposandneg.Consensus.Catalogue.narrowPeak $DATAPATH/peaks/"$f"/"$f"_sortedFE.bdg $OPDIR/$OPNAME/mapbed/LSCn.LSCp LSCn."$f".All.LSC.Consensus.mapbed
done
