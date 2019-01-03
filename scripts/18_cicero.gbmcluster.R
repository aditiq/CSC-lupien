#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
 #R-3.5.0
 #mordor
 #Objective : Run Cicero to identify co-accessible peaks in each CSC
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------
#load dependencies
#--------------------------------------------------------------------------------------------------------------------------------
(library(data.table))
suppressMessages(library(cicero))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))

#--------------------------------------------------------------------------------------------------------------------------------
#Run cicero
#--------------------------------------------------------------------------------------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
ds <- "data/ConsensusSet/PCSC1/GBM.Consensus.Mapbed.Qnorm.txt"
opname <- "GBM.clusters"

a <- read.table(ds, sep="\t", stringsAsFactors=F, header=T)
a$id <- paste(a[,1], a[,2],a[,3], sep="_")

## Group by PCA12.PCSC1.Consensus.wdlabel.v2
grp1 <- c("G566","G577","G729","G498","G797","G551")
grp3 <- c("G705","G719","G571")
gsub("GBM.", "",colnames(a)[grep("GBM",colnames(a))])

'%!in%' <- function(x,y)!('%in%'(x,y))

colnames(a) <- gsub("GBM.", "",colnames(a))

colnames(a) %!in% c(grp1,grp3)
g1.summ <- apply(a[,grp1],1,mean)
g3.summ <- apply(a[,grp3],1,mean)
g2.summ <- apply(a[,colnames(a)[grep("G",colnames(a)[colnames(a) %!in% c(grp1,grp3)])]],1,mean)

gbm.summ <- cbind(g1.summ , g2.summ ,g3.summ )
colnames(gbm.summ) <- c("grp1","grp2","grp3")
gbm.summ <- as.data.frame(gbm.summ)
gbm.summ$id  <- a$id

cicero_data <- melt(gbm.summ, id.vars=c("id"))
cicero_cds <- make_atac_cds(cicero_data, binarize = FALSE)
exprs(cicero_cds) <- as.matrix(exprs(cicero_cds))

set.seed(2017)
saveRDS(cicero_cds, file=paste0("results/PCSC1/cicero/", opname, ".cicero_cds.rds") )
cicero_cds <- readRDS(paste0("results/PCSC1/cicero/", opname, ".cicero_cds.rds"))
human.hg38.genome <- read.table("../common.data/hg38_ChromInfo.txt")
connsfull <- run_cicero(cicero_cds, human.hg38.genome)
saveRDS(connsfull, file=paste0("results/PCSC1/cicero/", opname, ".connsfull.rds") )

CCAN_assigns <- generate_ccans(connsfull)
saveRDS(CCAN_assigns, file=paste0("results/PCSC1/cicero/", opname, ".CCAN_assigns.rds") )

write.table(connsfull, file=paste0("results/PCSC1/cicero/", opname, ".connsfull.bed"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(CCAN_assigns, file=paste0("results/PCSC1/cicero/", opname, ".CCAN_assigns.bed"), sep="\t", row.names=F, col.names=T, quote=F)

