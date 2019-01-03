#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify cluster specific peaks using scABC's function on GBM and HF
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
source("scripts/scABC.cluster.R") #source("https://raw.githubusercontent.com/timydaley/scABC/master/R/cluster.R")
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)

#----------------------------------------------------------------
# Define function to run scABC on different pairwise comparisons
#----------------------------------------------------------------

#--------------------------------
## Read in binary matrix
#--------------------------------
binarymatds="data/ConsensusSet/PCSC1/GBM.HF.Consensus.Catalogue.Binarymat.txt"
binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

#--------------------------------
## Annotations
#--------------------------------

bg <- c(rep(1,26))
gbm.hf.cl <- c(rep(1,23),rep(2,3))

names(bg) <- colnames(ncounts)
names(gbm.hf.cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
gbm.hf.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =gbm.hf.cl, background_medians = bg)
save(gbm.hf.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.GBMvsHF.peakselection.Rdata"))

#--------------------------------
## Extract P values
#--------------------------------

gbm.hf.scabc.pval = as.data.frame(gbm.hf.scabc$pvalue)
rownames(gbm.hf.scabc.pval) <- rownames(ncounts)
colnames(gbm.hf.scabc.pval) <- c("GBM","HF")
gbm.hf.scabc.pval$id <- rownames(gbm.hf.scabc.pval)
gbm.hf.scabc.pval <- gbm.hf.scabc.pval[order(gbm.hf.scabc.pval$id),]

#--------------------------------
## Combine dataset
#--------------------------------
opname="GBMvsHF"
write.table(gbm.hf.scabc.pval, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

#--------------------------------
## Plot images
#--------------------------------
gbmsp <- subset(gbm.hf.scabc.pval$id, gbm.hf.scabc.pval$GBM < 0.05)
hfsp <- subset(gbm.hf.scabc.pval$id, gbm.hf.scabc.pval$HF < 0.05)

pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
image(t(apply(as.matrix(ncounts[c(gbmsp,hfsp),]),2,rev)), col=scalered)
dev.off()
