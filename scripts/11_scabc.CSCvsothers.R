#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify cluster specific peaks using scABC's function on CSC vs others
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

binarymatds="data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt"
binarymat <- fread(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)
colnames(ncounts) <- gsub("_peaks","", colnames(ncounts))

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[colnames(ncounts),]
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2)


#--------------------------------
## Annotations
#--------------------------------

bg <- c(rep(1,208))
csc.nonstem.cl <- mapping$stem

names(bg) <- colnames(ncounts)
names(csc.nonstem.cl) <- colnames(ncounts)

#--------------------------------
## Run scabc 
#--------------------------------
csc.nonstem.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =csc.nonstem.cl, background_medians = bg)
save(csc.nonstem.scabc, file=paste0("results/KitchenSink2/Cluster/scABC/scABC.CSCvsOthers.peakselection.Rdata"))

#--------------------------------
## Extract P values
#--------------------------------
csc.nonstem.scabc.pval = as.data.frame(csc.nonstem.scabc$pvalue)
rownames(csc.nonstem.scabc.pval) <- rownames(ncounts)
colnames(csc.nonstem.scabc.pval) <- c("CSC","Others")
csc.nonstem.scabc.pval$id <- rownames(csc.nonstem.scabc.pval)
csc.nonstem.scabc.pval <- csc.nonstem.scabc.pval[order(csc.nonstem.scabc.pval$id),]

#--------------------------------
## Combine dataset
#--------------------------------
opname="CSCvsOthers"
write.table(csc.nonstem.scabc.pval, file=paste0("results/KitchenSink2/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

#--------------------------------
## Plot images
#--------------------------------
csc.sp <- subset(csc.nonstem.scabc.pval$id, csc.nonstem.scabc.pval$CSC < 0.00001)
non.csc.sp <- subset(csc.nonstem.scabc.pval$id, csc.nonstem.scabc.pval$Others < 0.00001)

mapping2 <- mapping[order(mapping$stem),]
ncounts2 <- ncounts[,(mapping2$sample)]

pdf(paste0("results/KitchenSink2/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
image(t(apply(as.matrix(ncounts2[c(csc.sp ,non.csc.sp),]),2,rev)), col=scalered)
dev.off()

