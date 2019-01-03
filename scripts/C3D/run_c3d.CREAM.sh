
module load bedtools/2.23.0
module load R/3.4.1

## Generate signal matrix
files="PCSC1
KitchenSink2
KitchenSink1" ;

for f in $files ;
do
  for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr1.bed;
  do
      name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
      qsub -cwd -q light.q -e stdout/c3d -o stdout/c3d -N "$f"."$name".c3d -b y sh scripts/C3D/c3d.sh \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/CREAM/$f.CREAM.bed.sorted" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt" "$k" \
      "/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/$f.$name/" "";
  done
done
