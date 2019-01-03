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

binarymatds="data/ConsensusSet/PCSC1/LSCp.LSCb.Consensus.Catalogue.Binarymat.txt"
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

bg <- c(rep(1,32))
lscp.lscb.cl <- as.numeric(gsub("LSC.Bulk",1,gsub("LSC.positive",2,anno2$grp2)))

names(bg) <- colnames(ncounts)
names(lscp.lscb.cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
lscp.lscb.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =lscp.lscb.cl, background_medians = bg)
save(lscp.lscb.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSCpLSCbulk.peakselection.Rdata"))
