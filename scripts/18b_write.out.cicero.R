#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
 #R-3.5.0
 #mordor
 #Objective : Analyse Cicero to identify chromatin accessibility hubs in PCSC1
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------
#load dependencies
#--------------------------------------------------------------------------------------------------------------------------------
(library(data.table))
suppressMessages(library(cicero))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))

opname=c("PCSC1","GBM","PFA","LSCp")

for ( f in opname){

  connsfull <- readRDS(paste0("results/PCSC1/cicero/", f, ".connsfull.rds"))
  CCAN_assigns <- readRDS(paste0("results/PCSC1/cicero/", f, ".CCAN_assigns.rds") )
  write.table(connsfull, file=paste0("results/PCSC1/cicero/", f, ".connsfull.bed"), sep="\t", row.names=F, col.names=T, quote=F)
  write.table(CCAN_assigns, file=paste0("results/PCSC1/cicero/", f, ".CCAN_assigns.bed"), sep="\t", row.names=F, col.names=T, quote=F)

}

# ## CICERO GENE ACTIVITY
# gene_annotation_sample <- fread("../common.data/gencode.v24.for.cicero.bed", header=F,data.table=F, sep="\t", stringsAsFactors=F)
# colnames(gene_annotation_sample) <- c("chromosome","start","end","strand","feature","gene","transcript","symbol")
# gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]
# names(gene_annotation_sub)[4] <- "gene"
#
# opname="GBM"
# cicero_cds <- readRDS(paste0("results/PCSC1/cicero/", opname, ".cicero_cds.rds"))
# cicero_cds <- annotate_cds_by_site(cicero_cds, gene_annotation_sub)
# connsfull <-  readRDS(paste0("results/PCSC1/cicero/", opname, ".connsfull.rds") )
#
# unnorm_ga <- build_gene_activity_matrix(cicero_cds, connsfull)
# num_genes <- pData(cicero_cds)$num_genes_expressed
# names(num_genes) <- row.names(pData(input_cds))
#
# # normalize
# cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
