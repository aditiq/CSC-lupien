# Objective : Run homer on common, shared and not common regions identified by C3D


DIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/finalresults/" ;

bgfiles="data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak
data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.narrowPeak
data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak"

files="GBM.common.coaccess.regions.bed
GBM.notshared.coaccess.regions.bed
GBM.shared.coaccess.regions.bed
LSCp.common.coaccess.regions.bed
LSCp.notshared.coaccess.regions.bed
LSCp.shared.coaccess.regions.bed
PFA.common.coaccess.regions.bed
PFA.notshared.coaccess.regions.bed
PFA.shared.coaccess.regions.bed"

for f in $files ;
do

  name=$(echo $f | cut -d"." -f1)
  dirname=$(echo $f | sed -e 's/.coaccess.regions.bed//g')

    for j in $bgfiles ;
    do
        bgname=$(echo $j | awk -F"/" '{print $NF}' | sed -e 's/.Consensus.Catalogue.narrowPeak //g')
        mkdir $DIR/homer/$dirname.$bgname

        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $j  > $DIR/homer/$dirname.$bgname/tmp.bg.bed
        awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3, "Peak."NR }' $DIR/$f  > $DIR/homer/$dirname.$bgname/tmp.bed

        qsub -q lupiengroup -cwd -N Homer."$dirname".$bgname  -e stdout/homer/c3d/ -o stdout/homer/c3d/ \
        ../common.scripts/callhomer.hg38.bg.sh \
        $DIR/homer/$dirname.$bgname/tmp.bed $DIR/homer/$dirname.$bgname/ $DIR/homer/$dirname.$bgname/tmp.bg.bed
    done
done
