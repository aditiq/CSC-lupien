#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Get readcounts overlappign consensus peaks set of GBM , PFA and LSC
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


module load R/3.4.1
module load bedtools/2.23.0

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
    qsub -cwd -N CoverageBed.LSCp.Consensus.$f -e stdout/covbed/lscp -o stdout/covbed/lscp  -q light.q $SCRIPTDIR/generic_coverageBed.sh \
    $DATAPATH/bams/LSC/"$f"/"$f".filt.srt.dedup.q30.bam $lscpcons $OPDIR/$OPNAME/readcounts/lscp/ LSCp."$f".Consensus.coveragebed
done

for f in $gbmfiles ;
do
    qsub -cwd -N CoverageBed.GBM.Consensus.$f -e stdout/covbed/gbm -o stdout/covbed/gbm  -q light.q $SCRIPTDIR/generic_coverageBed.sh \
    $DATAPATH/bams/GBM/"$f"/"$f".filt.srt.dedup.q30.bam $gbmcons $OPDIR/$OPNAME/readcounts/gbm/ GBM."$f".Consensus.coveragebed
done

for f in $pfafiles ;
do
    qsub -cwd -N CoverageBed.PFA.Consensus.$f -e stdout/covbed/pfa -o stdout/covbed/pfa  -q light.q $SCRIPTDIR/generic_coverageBed.sh \
    $DATAPATH/bams/PFA/"$f"/"$f".filt.srt.dedup.q30.bam $pfacons $OPDIR/$OPNAME/readcounts/pfa/ PFA."$f".Consensus.coveragebed
done
