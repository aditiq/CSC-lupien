
module load bedtools/2.23.0
module load R/3.4.1

## Generate signal matrix
files="G523
G583
G567";

for f in $files ;
do
  for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr1.bed;
  do
      name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
      mkdir -p results/PCSC1/C3D/$f.$name
      qsub -cwd -q light.q -e stdout/c3d -o stdout/c3d -N "$f"."$name".c3d -b y sh scripts/C3D/c3d.sh \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/peaks/$f/*Peak" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt" "$k" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/" "";
  done
done


## Combine signal matrix -- commented out everything else in the R script
files="G523
G583
G567";

for f in $files ;
do
  for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr*.bed;
  do
      name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
      qsub -cwd -q all.q -e stdout/c3d -o stdout/c3d -N "$f"."$name".c3d -b y sh temp.sh \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/$f.Consensus.Catalogue.narrowPeak.sorted" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt" "$k" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/" "";
  done
done


## Run using signal matrix files
## commented out everything else

files="G523
G583
G567";

for f in $files ;
do
  for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr*;
  do

    name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
    timestamp() {
      date +"%Y-%m-%d_%H-%M-%S"
    }

    qsub -cwd -q all.q -e stdout/c3d -o stdout/c3d -N Run.C3D."$f"."$name".c3d -b y \
    /mnt/work1/software/R/3.4.1/Rscript /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts/C3D/c3d.full.r \
    "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/" \
    "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/" \
    "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/anchors.bed" \
    "/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt" \
    $(timestamp) \
    "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/"$f".chr1/signalMatrix.txt";
  done
done
