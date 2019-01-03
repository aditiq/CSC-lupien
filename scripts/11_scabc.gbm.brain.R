#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify cluster specific peaks using scABC's function on different comparisons
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

binarymatds="data/ConsensusSet/PCSC1/GBM.BRAIN.Consensus.Catalogue.Binarymat.txt"
binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

#--------------------------------
## Annotations
#--------------------------------
bg <- c(rep(1,53))
cl <- c(rep(1,8), rep(2, 23), rep(1, 22))
names(bg) <- colnames(ncounts)
names(cl) <- colnames(ncounts)

#--------------------------------
## Run scabc
#--------------------------------
gbm.brain.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =cl, background_medians = bg)
save(gbm.brain.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.GBM.BRAIN.peakselection.Rdata"))

#--------------------------------
## Extract P values
#--------------------------------
gbm.brain.scabc.pval = as.data.frame(gbm.brain.scabc$pvalue)
rownames(gbm.brain.scabc.pval) <- rownames(ncounts)
colnames(gbm.brain.scabc.pval) <- c("Brain","GBM")
gbm.brain.scabc.pval$id <- rownames(gbm.brain.scabc.pval)
gbm.brain.scabc.pval <- gbm.brain.scabc.pval[order(gbm.brain.scabc.pval$id),]

#--------------------------------
## Combine dataset
#--------------------------------
opname="GBM.BRAIN"
write.table(gbm.brain.scabc.pval, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

#--------------------------------
## Plot images
#--------------------------------
brainsp <- subset(gbm.brain.scabc.pval$id, gbm.brain.scabc.pval$Brain < 0.05)
gbmsp <- subset(gbm.brain.scabc.pval$id, gbm.brain.scabc.pval$GBM < 0.05)

pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
image(t(apply(as.matrix(ncounts[c(brainsp,gbmsp),]),2,rev)), col=scalered)
dev.off()
