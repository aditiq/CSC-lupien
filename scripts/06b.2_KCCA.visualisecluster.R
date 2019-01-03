#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Plot clusters from KCCA
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
library(gplots)
library(pheatmap)

bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")
#-----------------------------------
# Load binary matrix
#-----------------------------------
lscmat <- read.table("data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.Binarymat.txt",check.names=F, stringsAsFactors = F, header=T, sep="\t")
rownames(lscmat) <- paste(lscmat[,1], lscmat[,2], lscmat[,3], sep="_")

#-----------------------------------
# Load KCCA dataset
#-----------------------------------
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.LSC.10.Rdata")
flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
rownames(flexclust.clus) <- paste(lscmat$seqnames, lscmat$start, lscmat$end,sep="_")

cc <- as.matrix(kcca.cl$centers)
colnames(cc) <- gsub("_peaks","",colnames(lscmat)[4:ncol(lscmat)])
rownames(cc) <- seq(1, max(kcca.cl$cluster),1)


pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/Figures/K10.LSC.Flexclust.pdf"))
pheatmap(as.matrix(cc), 
          col=hmcols, scale="none",
          cluster_rows=TRUE, cluster_cols=TRUE, 
          clustering_method="complete",
          clustering_distance_rows = "euclidean",
          trace="none")
dev.off()