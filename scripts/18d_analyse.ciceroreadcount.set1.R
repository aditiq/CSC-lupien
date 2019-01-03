#============================================================================================================================================================
# Objective: First set of analysis from Cicero results
#============================================================================================================================================================

    # 1. What is the distribution of the co-accessible regions ?
        # A. E-E , E-P and P-P

    # 2. Co-accessibility regions of LSC17 genes

#============================================================================================================================================================

### R-3.4.1
### mordor

# ==============================================================================
# load dependencies
# ==============================================================================
library(data.table)
library(stringr)

# ==============================================================================
# Run analysis
# ==============================================================================

analysfunc=function(name, opname, ...){

    ds <- fread(paste0("results/PCSC1/cicero/run2.readcount/", name,".connsfull.bed"), header=T, sep="\t", stringsAsFactors=F, data.table=F)

    ds.coaccess <- subset(ds, ds$coaccess >0) ## keep only co-accessible regions
    ds.coaccess$pairID <- paste(ds.coaccess$Peak1, ds.coaccess$Peak2, sep=":")

    ds.dedup <- fread(paste0("results/PCSC1/cicero/run2.readcount/", name,".deduped.connsfull.pairs.bed"), header=F, sep="\t", stringsAsFactors=F, data.table=F)
    ds.dedup$pairID <- paste(ds.dedup$V1, ds.dedup$V2, sep=":")

    # Keep only dedup pairs
    ds.dedup.coaccess <- subset(ds.coaccess,ds.coaccess$pairID %in% ds.dedup$pairID)

    # Identify promoter regions
    ds.prom <- fread(paste0("results/PCSC1/cicero/run2.readcount/", name,".promoters.bed"), header=F, sep="\t", stringsAsFactors=F, data.table=F) ## read in regions mapping to promoters
    colnames(ds.prom) <- c("chr","start","end","geneid","txtid","gene","strand","distance","id")

    ds.dedup.coaccess.anno <- merge(ds.dedup.coaccess, ds.prom[,c("id","geneid","txtid")], by.x="Peak1", by.y="id", all.x=T)
    colnames(ds.dedup.coaccess.anno) <- c("Peak1", "Peak2", "coaccess", "pairID","geneid.peak1", "txtid.peak1")

    ds.dedup.coaccess.anno <- merge(ds.dedup.coaccess.anno, ds.prom[,c("id","geneid","txtid")], by.x="Peak2", by.y="id", all.x=T)
    colnames(ds.dedup.coaccess.anno) <- c("Peak1", "Peak2", "coaccess", "pairID","geneid.peak1", "txtid.peak1", "geneid.peak2", "txtid.peak2")

    write.table(ds.dedup.coaccess.anno, file=paste0("results/PCSC1/cicero/run2.readcount/", name, ".Deduplicated.coaccess.Annotated.txt"), sep="\t", quote=F,row.names=F, col.names=T)

    ## Summary 1 -- How many E-E, E-P and P-P interactiobs
    ds.dedup.coaccess.anno[is.na(ds.dedup.coaccess.anno) ] <- "NA"

    summ1 <- ds.dedup.coaccess.anno[,c("pairID","geneid.peak1","geneid.peak2")]
    summ1$peak1.anno <- ifelse(summ1$geneid.peak1 =="NA", "E","P")
    summ1$peak2.anno <- ifelse(summ1$geneid.peak2 =="NA", "E","P")
    summ1 <- summ1[,c("pairID", "peak1.anno","peak2.anno")]
    summ1 <- summ1[!duplicated(summ1$pairID),]

    write.table(summ1, file=paste0("results/PCSC1/cicero/run2.readcount/", name, ".E-P.Anno.txt"), sep="\t", quote=F,row.names=F, col.names=T)
    summ1tab <- as.data.frame(table(summ1$peak1.anno, summ1$peak2.anno))
    colnames(summ1tab)[1:2] <- c("Peak1", "Peak2")
    write.table(summ1tab, file=paste0("results/PCSC1/cicero/run2.readcount/", name, ".E-P.Summary.txt"), sep="\t", quote=F,row.names=F, col.names=T)

    ## Co-accessibility hubs for LSC17 genes

    lscgene17 <- read.table("data/lsc17genes.txt", header=F, sep="\t", stringsAsFactors=F)
    ds.dedup.coaccess.anno.lsc17 <- subset(ds.dedup.coaccess.anno, ds.dedup.coaccess.anno$geneid.peak1 %in% lscgene17$V2 | ds.dedup.coaccess.anno$geneid.peak2 %in% lscgene17$V2 )
    write.table(ds.dedup.coaccess.anno.lsc17, file=paste0("results/PCSC1/cicero/run2.readcount/", name, ".LSCp17geneshubs.txt"), sep="\t", quote=F,row.names=F, col.names=T)

}

analysfunc("GBM")
analysfunc("PFA")
analysfunc("LSCp")
