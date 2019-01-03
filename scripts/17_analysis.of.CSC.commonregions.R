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

#awk 'BEGIN {FS=OFS="\t"} { if ($4==79) print $1,$2,$3}' results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Promoter.txt > results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Promoter.Common.ClusterNo79.bed

#awk 'BEGIN {FS=OFS="\t"} { if ($4==12) print $1,$2,$3}' results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Enhancer.txt > results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Enhancer.Common.ClusterNo12.bed

#awk 'BEGIN {FS=OFS="\t"} { if ($4==1) print $1,$2,$3}' results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Promoter.txt > results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Promoter.Common.ClusterNo1.bed

#awk 'BEGIN {FS=OFS="\t"} { if ($4==94) print $1,$2,$3}' results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Enhancer.txt > results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Enhancer.Common.ClusterNo94.bed

#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
library(UpSetR)

#----------------------------------------------------------------
# upset plot
#----------------------------------------------------------------

# kitchensink1p="results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Promoter.Common.ClusterNo1.bed"
# kitchensink1e="results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink1.Flexclust.Enhancer.Common.ClusterNo94.bed"
# kitchensink2p="results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Promoter.Common.ClusterNo79.bed"
# kitchensink2e="results/KitchenSink2/Cluster/KCCA.Flexclust/ExtractedClusters/K100.KitchenSink2.Flexclust.Enhancer.Common.ClusterNo12.bed"
# pcsc1e="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Enhancergrp.txt"
# pcsc1p="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Promotergrp.txt"
# pcsc1pe="temp/Common.PCSC1.bed"
# ks1pe="temp/K100.KitchenSink1.Flexclust.Enhancer.Promoter.Common.bed"

# cat $kitchensink1p $kitchensink1e $kitchensink2p $kitchensink2e $pcsc1e $pcsc1p | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' - | sortBed -i stdin | mergeBed -i stdin >> Common.bed
# cp $kitchensink1p  temp/
# cp $kitchensink1e  temp/
# cp $kitchensink2p  temp/
# cp $kitchensink2e  temp/
# cp $pcsc1e  temp/
# cp $pcsc1p  temp/
# rename txt bed temp/*txt
# Rscript scripts/createbinarymat.R Common.bed temp/ ".bed" ./ Common.binarymat.txt

#cat Common.Enhancergrp.bed Common.Promotergrp.bed > Common.PCSC1.bed
#cat K100.KitchenSink1.Flexclust.Promoter.Common.ClusterNo1.bed  K100.KitchenSink1.Flexclust.Enhancer.Common.ClusterNo94.bed > K100.KitchenSink1.Flexclust.Enhancer.Promoter.Common.bed

binarymat <- fread("Common.binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
colnames(binarymat)[4:11] <- c("Common.CSC.E", "Common.CSC","Common.CSC.P",
                                "Common.All.E","Common.All", "Common.All.P", "KS2.E", "KS2.P")

pdf("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/UpsetPlot.CommonEnhancer.wd.KitchenSink.pdf", useDingbats=F)
upset(binarymat[,c(4,7)],nsets=6, order.by = "freq")
dev.off()

pdf("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/UpsetPlot.CommonPromoter.wd.KitchenSink.pdf", useDingbats=F)
upset(binarymat[,c(6,9)],nsets=6, order.by = "freq")
dev.off()

## common for enhancer and promoters
pdf("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/UpsetPlot.CommonPCSC1.wd.KitchenSink.pdf", useDingbats=F)
upset(binarymat[,c(5,8)],nsets=6, order.by = "freq")
dev.off()
