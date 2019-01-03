#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Extract clusters from KCCA
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(data.table)
library(prabclus)
library("RColorBrewer")
library(gplots)
library(pheatmap)
library(ggplot2)
library(cluster)
library(GGally)
library(Rtsne)
library(rGREAT)
library(grid)
library(gridExtra)
library(ggthemes)
library(scales)
library(ggrepel)
library(stringr)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026')
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract clusters for kitchen sink
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------
# Load binary matrix
#-----------------------------------
binarymat <- fread("data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")

ksenh <- fread("data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Enhancers.bed", header=F, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
rownames(ksenh) <- paste(ksenh[,1],ksenh[,2],ksenh[,3],sep="_")
ksenhmat <- binarymat[rownames(ksenh),]
colnames(ksenhmat) <- gsub("_peaks","", colnames(ksenhmat))

ksprom <- fread("data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Promoters.bed", header=F, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
rownames(ksprom) <- paste(ksprom[,1],ksprom[,2],ksprom[,3],sep="_")
kspromat <- binarymat[rownames(ksprom),]
colnames(kspromat) <- gsub("_peaks","", colnames(kspromat))

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[colnames(ksenhmat)[4:ncol(ksenhmat)],]
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2)
mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2))

visfunc=function(dsname, rownamesds, mat, opname,opname2,opname3,... ){

  load(dsname) ## load dataset
  load(rownamesds)

  flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
  rownames(flexclust.clus) <- rownames.obj

  flexclust.clus2 <- cbind(as.data.frame(str_split_fixed(rownames(flexclust.clus), "_", 3)), flexclust.clus[,1])
  write.table(flexclust.clus2, row.names=F, col.names=F, quote=F, sep="\t", file=opname3 )

  cc <- as.matrix(kcca.cl$centers)
  colnames(cc) <- mapping$sample
  mapping2 <- mapping[order(mapping$stem),]
  rownames(cc) <- seq(1, max(kcca.cl$cluster),1)

  annocol <- data.frame(stringsAsFactors=F, stem=mapping$stem, tissue=mapping$tissue)
  rownames(annocol) <- colnames(cc)

  ## Plot heatmap unsuper
  pdf(opname)
  pheatmap(as.matrix(cc[,mapping2$sample]),
          col=scalered, scale="none",
          fontsize_row=1.5,fontsize_col=1.5,
          cluster_rows=TRUE, cluster_cols=FALSE,
          clustering_method="ward.D2",
          clustering_distance_rows = "correlation",
          annotation_col=annocol,
          trace="none")
  dev.off()

  ## Plot heatmap unsuper

  pcaPRComp <- prcomp(t(as.matrix(cc)))
  df_out <- as.data.frame(pcaPRComp$x)
  df_out$group <- mapping$name2

  p <- ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
  geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
  theme(legend.position="none") +
  geom_point()

  pdf(opname2)
  print(p)
  dev.off()

}


datadir="results/KitchenSink1/Cluster/KCCA.Flexclust/kccadata/"
figdir="results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/"

for (f in c(100,75,50,25)) {
  visfunc(
  dsname=paste0(datadir,"KCCA.KitchenSink1.Enhancer.",f, ".Rdata"),
  rownamesds=paste0(datadir,"KCCA.KitchenSink1.Enhancer.rownames.",f,".Rdata"),
  mat=ksenhmat,
  opname=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Enhancer.unsuper.pdf"),
  opname2=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Enhancer.unsuper.PCA.pdf"),
  opname3=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Enhancer.txt")
  )

  visfunc(
  dsname=paste0(datadir,"KCCA.KitchenSink1.Promoter.",f, ".Rdata"),
  rownamesds=paste0(datadir,"KCCA.KitchenSink1.Promoter.rownames.",f,".Rdata"),
  mat=kspromat,
  opname=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Promoter.unsuper.pdf"),
  opname2=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Promoter.unsuper.PCA.pdf"),
  opname3=paste0(figdir,"K", f, ".KitchenSink1.Flexclust.Promoter.txt")

  )
}
