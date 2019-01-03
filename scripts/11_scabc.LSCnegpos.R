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
binarymatds="data/ConsensusSet/PCSC1/LSCposandneg.Consensus.Catalogue.Binarymat.txt"
binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

anno <- read.table("data/Alignment.Stats.txt", header=T, sep="\t", stringsAsFactors = F)
anno$grp2 <- paste0(anno$Cancer,".", anno$Group)
rownames(anno) <- paste0(anno$Name, "_peaks")
anno2 <- anno[colnames(ncounts),]

#--------------------------------
## Annotations
#--------------------------------

bg <- c(rep(1,52))
lscp.lscn.cl <- as.numeric(gsub("LSC.negative",1,gsub("LSC.positive",2,anno2$grp2)))

names(bg) <- colnames(ncounts)
names(lscp.lscn.cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
lscp.lscn.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =lscp.lscn.cl, background_medians = bg)
save(lscp.lscn.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSCposvsLSCneg.peakselection.Rdata"))

#--------------------------------
## Extract P values
#--------------------------------
lscp.lscn.scabc.pval = as.data.frame(lscp.lscn.scabc$pvalue)
rownames(lscp.lscn.scabc.pval) <- rownames(ncounts)
colnames(lscp.lscn.scabc.pval) <- c("LSC.negative","LSC.positive")
lscp.lscn.scabc.pval$id <- rownames(lscp.lscn.scabc.pval)
lscp.lscn.scabc.pval <- lscp.lscn.scabc.pval[order(lscp.lscn.scabc.pval$id),]

#--------------------------------
## Combine dataset
#--------------------------------
opname="LSCnegvsLSCpos"
write.table(lscp.lscn.scabc.pval, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

#--------------------------------
## Plot images
#--------------------------------
lscnegsp <- subset(lscp.lscn.scabc.pval$id, lscp.lscn.scabc.pval$LSC.negative < 0.05)
lscpossp <- subset(lscp.lscn.scabc.pval$id, lscp.lscn.scabc.pval$LSC.positive < 0.05)

anno2 <- anno2[order(anno2$grp2),]
ncounts <- ncounts[,paste0(anno2$Name, "_peaks")]
pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
image(t(apply(as.matrix(ncounts[c(lscnegsp,lscpossp),]),2,rev)), col=scalered)
dev.off()