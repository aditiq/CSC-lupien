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

args <- commandArgs(trailingOnly = TRUE)
ds <- as.character(args[1])
opname <- as.character(args[2])

print(ds)
print(opname)

new_exprs <- read.table(ds, sep="\t",  header=T, stringsAsFactors=F)
rownames(new_exprs) <- paste(new_exprs[,1], new_exprs[,2],new_exprs[,3], sep="_")
new_exprs <- as.matrix(new_exprs[,4:ncol(new_exprs)])

cicero_cds2 <-  suppressWarnings(newCellDataSet(new_exprs,
                              #phenoData = pd,
                              #featureData = fd,
                              expressionFamily=negbinomial.size(),
                              lowerDetectionLimit=0))

cicero_cds2 <- BiocGenerics::estimateSizeFactors(cicero_cds2)
cicero_cds2 <- suppressWarnings(BiocGenerics::estimateDispersions(cicero_cds2))
Biobase::exprs(cicero_cds2) <- t(t(Biobase::exprs(cicero_cds2))/Biobase::pData(cicero_cds2)$Size_Factor)

set.seed(2017)
saveRDS(cicero_cds2, file=paste0("results/PCSC1/cicero/", opname, ".cicero_cds.readcount.rds") )
cicero_cds2 <- readRDS(paste0("results/PCSC1/cicero/", opname, ".cicero_cds.readcount.rds"))
human.hg38.genome <- read.table("../common.data/hg38_ChromInfo.txt")

tmp <- readRDS(paste0("results/PCSC1/cicero/", opname, ".cicero_cds.rds"))
fData(cicero_cds2) <- fData(tmp)

connsfull2 <- run_cicero(cicero_cds2, human.hg38.genome)
saveRDS(connsfull2, file=paste0("results/PCSC1/cicero/", opname, ".connsfull.readcount.rds") )
write.table(connsfull2, file=paste0("results/PCSC1/cicero/", f, ".connsfull.readcount.bed"), sep="\t", row.names=F, col.names=T, quote=F)

#CCAN_assigns2 <- generate_ccans(connsfull2)
#saveRDS(CCAN_assigns2, file=paste0("results/PCSC1/cicero/", opname, ".CCAN_assigns.readcount.rds") )
