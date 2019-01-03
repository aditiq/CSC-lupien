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

binarymatds="data/ConsensusSet/PCSC1/LSCp.HEMATHSC.Consensus.Catalogue.Binarymat.txt"
binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

#--------------------------------
## Annotations
#--------------------------------

bg <- c(rep(1,25))
cl <- c(rep(1,18), rep(2,7))

names(bg) <- colnames(ncounts)
names(cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
LSCp.HEMATHSC.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =cl, background_medians = bg)
save(LSCp.HEMATHSC.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSCp.HEMATHSC.peakselection.Rdata"))
