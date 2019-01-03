



ALIGNMENTSTATS="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/Alignment.Stats.txt"
OPDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/"
OPNAME="PCSC1"
DATAPATH="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/"

gbmfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="GBM" && $10==1) print $1}' $ALIGNMENTSTATS)
pfafiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="PFA" && $10==1) print $1}' $ALIGNMENTSTATS)
lscfiles=$( awk 'BEGIN {FS=OFS="\t"} { if (NR>1 && $2=="positive" && $3=="LSC" && $10==1) print $1}' $ALIGNMENTSTATS)

for f in $gbmfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak  $OPDIR/$OPNAME/temp.peaks ; done
for f in $pfafiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done
for f in $lscfiles ; do cp $DATAPATH/peaks/$f/*.narrowPeak $OPDIR/$OPNAME/temp.peaks ; done

## Generate binary matrix
qsub -cwd -N binarymat.KS2.CREAM -q lupiengroup -e stdout/ -o stdout/ -b y Rscript scripts/createbinarymat.R results/PCSC1/CREAM/KitchenSink2.CREAM.bed.sorted data/ConsensusSet/KitchenSink2/temp.peaks ".narrowPeak" results/PCSC1/CREAM/ KitchenSink2.CREAM.Binarymat.txt

qsub -cwd -N binarymat.PCSC1.CREAM -q lupiengroup -e stdout/ -o stdout/ -b y Rscript scripts/createbinarymat.R results/PCSC1/CREAM/PCSC1.CREAM.bed.sorted data/ConsensusSet/PCSC1/temp.peaks ".narrowPeak" results/PCSC1/CREAM/ PCSC1.CREAM.Binarymat.txt

## Run kcca
f=100
qsub -cwd -N KCCA.Kitchensink1.CREAM."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "results/PCSC1/CREAM/KitchenSink1.CREAM.Binarymat.txt" "KitchenSink1.CREAM" $f ;
qsub -cwd -N KCCA.Kitchensink2.CREAM."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "results/PCSC1/CREAM/KitchenSink2.CREAM.Binarymat.txt" "KitchenSink2.CREAM" $f ;
qsub -cwd -N KCCA.PCSC1.CREAM."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "results/PCSC1/CREAM/PCSC1.CREAM.Binarymat.txt" "PCSC1.CREAM" $f ;


