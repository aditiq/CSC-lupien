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

a <- read.table(ds, sep="\t", stringsAsFactors=F, header=T)
a$id <- paste(a[,1], a[,2],a[,3], sep="_")
cicero_data <- melt(a[,4:ncol(a)], id.vars=c("id"))
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
