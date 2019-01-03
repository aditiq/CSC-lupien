#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Analysis of common/ shared regions between CSCs and other tissues.
              # Are common CSCs regions really common to CSCs or are they also found in other tissues ?

# Conclusion : Yes, more than 95% of the enhancers and promoters identified as common across Cancer stem cells actually are accessible in related and unrelated tissues.
#              When analysed with just ENCODE DNAse data - the % was around 50% but the remaining regions did not show any clear signals. This clears that confusion.
#              Maybe because of differences in DNAse and ATAC assay
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


# Files required
#  - KCCA - common to CSC - PCSC1
#  - KCCA - common to All - KitchenSink1
#  - KCCA - common to CSCs and related tissues - KitchenSink2
#  - cisTopic - common to CSCs and related tissues ( kind of ?) - KitchenSink2

#awk 'BEGIN {FS=OFS="\t"} { if ($4==93) print $1,$2,$3}' results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink2.CREAM.Flexclust.txt > results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink2.CREAM.Flexclust.Common.ClusterNo93.bed

#awk 'BEGIN {FS=OFS="\t"} { if ($4==9) print $1,$2,$3}' results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.CREAM.Flexclust.txt > results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.CREAM.Flexclust.Common.ClusterNo9.bed

#awk 'BEGIN {FS=OFS="\t"} { if ($4==89) print $1,$2,$3}' results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/PCSC1.CREAM.Flexclust.txt > results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/PCSC1.CREAM.Flexclust.Common.ClusterNo89.bed


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
library(UpSetR)

#----------------------------------------------------------------
# upset plot
#----------------------------------------------------------------

# kitchensink2="results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink2.CREAM.Flexclust.Common.ClusterNo93.bed"
# kitchensink1="results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.CREAM.Flexclust.Common.ClusterNo9.bed"
# pcsc1="results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/PCSC1.CREAM.Flexclust.Common.ClusterNo89.bed"

# cat $kitchensink2 $kitchensink1 $pcsc1 | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' - | sortBed -i stdin | mergeBed -i stdin >> temp.cream/Common.bed
# cp $kitchensink2  temp.cream/
# cp $kitchensink1  temp.cream/
# cp $pcsc1  temp.cream/

# Rscript scripts/createbinarymat.R temp.cream/Common.bed temp.cream/ ".bed" temp.cream/ Common.binarymat.txt

binarymat <- fread("temp.cream/Common.binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
binarymat$Common <- NULL
colnames(binarymat)[4:ncol(binarymat)] <- c("Common.KitchenSink1", "Common.KitchenSink2", "Common.PCSC1")

pdf("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/UpsetPlot.Common.CREAM.pdf", useDingbats=F)
upset(binarymat[,c(4:ncol(binarymat))],nsets=3, order.by = "freq")
dev.off()
