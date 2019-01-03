#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.5.0
# mordor
# Objective : Run Cicero to identify co-accessible peaks in each CSC
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

module load R/3.5.0

files="GBM
PFA
LSCp"

for f in $files ;
do
    qsub -cwd -N cicero."$f" -q lupiengroup \
    -e stdout/cicero -o stdout/cicero \
    -b y /mnt/work1/software/R/3.5.0/Rscript scripts/18_cicero.R data/ConsensusSet/PCSC1/"$f".Consensus.Mapbed.Qnorm.txt "$f"
done

files="PCSC1"
for f in $files ;
do
    qsub -cwd -N cicero."$f" -q lupiengroup \
    -e stdout/cicero -o stdout/cicero \
    -b y /mnt/work1/software/R/3.5.0/Rscript scripts/18_cicero.R data/ConsensusSet/PCSC1/"$f".Consensus.Mapbed.Qnorm.txt "$f"
done
