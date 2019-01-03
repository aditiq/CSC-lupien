#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
#mordor
#Objective: Alternative promoter usage in CSC compared to differentiated tissues
#           For a list of genes-txts in gencode v24, do diff vs CSC show atac peaks at same txt or different txts of the same gene ?
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

module load bedtools/2.23.0

COMMONDIR="/mnt/work1/users/lupiengroup/People/qamraa99/common.data"
GENCODEFILE=$COMMONDIR"/c3d.anchors.gencodev24.500bparoundTSS.bed"

SCRIPTDIR="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts";
Rscript $SCRIPTDIR/createbinarymat.R $GENCODEFILE data/ConsensusSet/KitchenSink2/temp.peaks ".narrowPeak" results/KitchenSink2/alternatepromoters/ KitchenSink2.GencodeTSS.Binarymat.txt
