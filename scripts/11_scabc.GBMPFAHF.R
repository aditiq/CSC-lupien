#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify cluster specific peaks using scABC's function on GBM, PFA and HF
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
binarymatds="data/ConsensusSet/PCSC1/GBM.PFA.HF.Consensus.Catalogue.Binarymat.txt"
binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

#--------------------------------
## Annotations
#--------------------------------

bg <- c(rep(1,33))
gbmpfa.hf.cl <- c(rep(1,23),rep(2,3), rep(1,7))

names(bg) <- colnames(ncounts)
names(gbmpfa.hf.cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
gbmpfa.hf.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =gbmpfa.hf.cl, background_medians = bg)
save(gbmpfa.hf.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.GBMPFAvsHF.peakselection.Rdata"))

#--------------------------------
## Extract P values
#--------------------------------

gbmpfa.hf.scabc.pval = as.data.frame(gbmpfa.hf.scabc$pvalue)
rownames(gbmpfa.hf.scabc.pval) <- rownames(ncounts)
colnames(gbmpfa.hf.scabc.pval) <- c("GBMPFA","HF")
gbmpfa.hf.scabc.pval$id <- rownames(gbmpfa.hf.scabc.pval)
gbmpfa.hf.scabc.pval <- gbmpfa.hf.scabc.pval[order(gbmpfa.hf.scabc.pval$id),]

#--------------------------------
## Combine dataset
#--------------------------------
opname="GBMPFAvsHF"
write.table(gbmpfa.hf.scabc.pval, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

#--------------------------------
## Plot images
#--------------------------------
gbmsp <- subset(gbmpfa.hf.scabc.pval$id, gbmpfa.hf.scabc.pval$GBMPFA < 0.05)
hfsp <- subset(gbmpfa.hf.scabc.pval$id, gbmpfa.hf.scabc.pval$HF < 0.05)

pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
image(t(apply(as.matrix(ncounts[c(gbmsp,hfsp),]),2,rev)), col=scalered)
dev.off()
