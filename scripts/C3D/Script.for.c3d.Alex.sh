## Main Script to Kick off C3D

#--------------------------
## README FIRST
#--------------------------

#-- This script uses
  #-- c3d.sh
  #-- c3d.R
  #-- c3d.full.r
      #-- Make sure the options specified in c3d.full.r  are what you generated ( lines 34 - 50)
      #-- e.g. genome build for tracks, window size, correlation thresholds etc

# Bedgraph files for hg38 DNAse data -- /mnt/work1/users/lupiengroup/Projects/CommonData/ENCODE/hg38/dnase.hg38.bedgraph.txt
# for atac data hg38 -- /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt

# Anchor file for gencode v24 -- /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3d.anchors.gencodev24.500bparoundTSS.bed


#--------------------------
## Load Dependencies
#--------------------------
module load bedtools/2.23.0
module load R/3.4.1

## Step 1 -- Generate and combine signal matrix + anchors.bed file in a format suitable for running c3d
## Don't run this if you already have a combined signal matrix file
qsub -cwd -q light.q -e stdout/c3d -o stdout/c3d -N $scriptname -b y sh c3d.sh $referencebedfile $databasebedgraphfiles $anchorfile $outputdirectory "";


## Step 2 -- Run using the generated signal matrix Files

timestamp() {
  date +"%Y-%m-%d_%H-%M-%S"
}

qsub -cwd -q light.q -e stdout/c3d -o stdout/c3d -N Run.C3D."$f"."$name".c3d -b y \
Rscript c3d.full.r \
$outputdirectoryfromstep1 \ ## Can be anything if you already have the signalmatrixfile
$opdirforresults \
$fullpathforanchorfile \ ## if you ran step1 then it will be the $outputdirectory specified in step1
$databasebedgraphfiles \ ## same as step1 -- also specified in README
$(timestamp) \
$fullpathtosignalmatrixfile
