### R-3.4.1
### mordor
### Objective : Exploratory analysis of binary matrix

### to do --  --- REALLY NEED TO MAKE A FUNCTION !!!
### too many of the same plots repeated over and over again !!!!!!!!!!

######################################
### load dependencies
######################################
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

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
#library(PharmacoGx)
library(data.table)
brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026')

source("scripts/binarymatanalysis.funcs.R")


######################################
### Define function
######################################

runfunc=function(ds, opname,cellTypes,...){

  pdf(paste0("results/PCSC1/BinaryMat/MDS.",opname,".pdf"), useDingbats = F)
  MDSClust <- plotMDS(ds, k=3, groups=cellTypes, ret.val=TRUE, text.label=FALSE, title="")
  dev.off()

  pamx <- pam(as.dist(getJaccardDist(ds[,4:dim(ds)[2]])), k=3, diss=TRUE)
  clusterInfo <- data.frame(CellName = names(pamx$clustering), Cluster=pamx$clustering)

  CellClust <- paste0('Cluster-',as.character(pamx$clustering))
  pdf(paste0("results/PCSC1/BinaryMat/Kmedoids.MDS",opname,".pdf"), useDingbats = F)
  MDSClust <- plotMDSClust(ds, 4, Clusters=CellClust, ret.val=TRUE, text.label=FALSE, title="")
  dev.off()

  ## PCA

  pdf(paste0("results/PCSC1/BinaryMat/PCA.",opname,".pdf"), useDingbats = F)
  plotMultiplePCAJaccard <- plotMultiplePCAJaccard(ds,nPCAToDisplay=3, groups=cellTypes ,ret.val=FALSE, text.label=FALSE, title="")
  dev.off()

  #pdf(paste0("results/PCSC1/BinaryMat/VarPCA.",opname,".pdf"))
  #plotVarExplained <- plotVarExplained(ds,nPCAToDisplay=3, groups=cellTypes ,ret.val=FALSE, text.label=FALSE, title="")
  #dev.off()


  ## t-sne
  pdf(paste0("results/PCSC1/BinaryMat/Tsne.",opname,".pdf"), useDingbats = F)
  plot_tSNE <- plot_tSNE(ds,nPCAToUSE=3, groups=cellTypes ,perplexity_division=8,ret.val=FALSE, text.label=FALSE, title="")
  dev.off()

  ## Heatmap
  annocol <- data.frame(cellTypes)
  rownames(annocol) <- colnames(ds[,4:ncol(ds)])
  pdf(paste0("results/PCSC1/BinaryMat/Heatmap.",opname,".pdf"), useDingbats = F)
  pheatmap(getJaccardDist(ds),scale="none",
           annotation_row = annocol,
           annotation_col = annocol,
           col=rev(colors),
           cluster_rows=T,
           cluster_cols=T,cex=0.5)
  dev.off()

}

############################################################################################
######################################### Run analysis######################################
############################################################################################


##----------------------------------------
## PCSC1
##----------------------------------------
### R-3.4.1
### mordor
### Objective : Exploratory analysis of binary matrix

### to do --  --- REALLY NEED TO MAKE A FUNCTION !!!
### too many of the same plots repeated over and over again !!!!!!!!!!

######################################
### load dependencies
######################################
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

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
#library(PharmacoGx)
library(data.table)
brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026')

source("scripts/binarymatanalysis.funcs.R")


######################################
### Define function
######################################

runfunc=function(ds, opname,cellTypes,...){

  pdf(paste0("results/PCSC1/BinaryMat/MDS.",opname,".pdf"), useDingbats = F)
  MDSClust <- plotMDS(ds, k=3, groups=cellTypes, ret.val=TRUE, text.label=FALSE, title="")
  dev.off()

  pamx <- pam(as.dist(getJaccardDist(ds[,4:dim(ds)[2]])), k=3, diss=TRUE)
  clusterInfo <- data.frame(CellName = names(pamx$clustering), Cluster=pamx$clustering)

  CellClust <- paste0('Cluster-',as.character(pamx$clustering))
  pdf(paste0("results/PCSC1/BinaryMat/Kmedoids.MDS",opname,".pdf"), useDingbats = F)
  MDSClust <- plotMDSClust(ds, 4, Clusters=CellClust, ret.val=TRUE, text.label=FALSE, title="")
  dev.off()

  ## PCA

  pdf(paste0("results/PCSC1/BinaryMat/PCA.",opname,".pdf"), useDingbats = F)
  plotMultiplePCAJaccard <- plotMultiplePCAJaccard(ds,nPCAToDisplay=3, groups=cellTypes ,ret.val=FALSE, text.label=FALSE, title="")
  dev.off()

  #pdf(paste0("results/PCSC1/BinaryMat/VarPCA.",opname,".pdf"))
  #plotVarExplained <- plotVarExplained(ds,nPCAToDisplay=3, groups=cellTypes ,ret.val=FALSE, text.label=FALSE, title="")
  #dev.off()


  ## t-sne
  pdf(paste0("results/PCSC1/BinaryMat/Tsne.",opname,".pdf"), useDingbats = F)
  plot_tSNE <- plot_tSNE(ds,nPCAToUSE=3, groups=cellTypes ,perplexity_division=8,ret.val=FALSE, text.label=FALSE, title="")
  dev.off()

  ## Heatmap
  annocol <- data.frame(cellTypes)
  rownames(annocol) <- colnames(ds[,4:ncol(ds)])
  pdf(paste0("results/PCSC1/BinaryMat/Heatmap.",opname,".pdf"), useDingbats = F)
  pheatmap(getJaccardDist(ds),scale="none",
           annotation_row = annocol,
           annotation_col = annocol,
           col=rev(colors),
           cluster_rows=T,
           cluster_cols=T,cex=0.5)
  dev.off()

}

############################################################################################
######################################### Run analysis######################################
############################################################################################


##----------------------------------------
## PCSC1
##----------------------------------------

## Load binary mat
binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)

cellNames <- colnames(binarymat[,4:ncol(binarymat)])
cellTypes <- factor(c( rep("LSC",18), rep("GBM", 23), rep("PFA",7)    ))
cellCol <- (c( rep("darkgoldenrod1",18), rep("lightblue", 23), rep("green",7)    ))

rownames(binarymat) <- paste(binarymat$seqname, binarymat$start, binarymat$end, sep="_")

## Entire Consensus
runfunc(binarymat, "PCSC1.Consensus",cellTypes=cellTypes)

## Add other annotation to heatmap
annocol <- data.frame(cellTypes)
rownames(annocol) <- gsub("_peaks","",colnames(binarymat[,4:ncol(binarymat)]))
annocol$id <- rownames(annocol)
stats <- read.table("data/Alignment.Stats.txt", header=T, sep="\t", stringsAsFactors=F)
annocol <- merge(annocol,stats,by.x="id", by.y="Name")
rownames(annocol) <- annocol$id
annocol$Readsrank <- (rank(annocol$Final.filtered.and.deduped.reads))
annocol$Peaksrank <- (rank(annocol$Peaks))
annocol <- annocol[,c(2,12,13)]

opname="PCSC1.Consensus"
pdf(paste0("results/PCSC1/BinaryMat/Heatmap.",opname,".V2.pdf"), useDingbats = F)
pheatmap(annocol[,2:3],scale="none",
         #annotation_row = annocol,
         #annotation_col = annocol,
         col=rev(colors),
         cluster_rows=T,
         cluster_cols=T,cex=0.5)
dev.off()

binarymat <- fread("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
colnames(binarymat) <- gsub("_peaks", "", colnames(binarymat))
cellNames <- colnames(binarymat[,4:ncol(binarymat)])

mapping <- read.delim("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[cellNames,]
rownames(binarymat) <- paste(binarymat$seqname, binarymat$start, binarymat$end, sep="_")
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,0)

opname="PCSC1.Consensus"
jacc <- (getJaccardDist(binarymat[,4:dim(binarymat)[2]]))
pcaPRComp <- prcomp(t(jacc))
df_out <- as.data.frame(pcaPRComp$x)
df_out$group <- mapping$name2


pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$tissue) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".v3.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$study) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".stemlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/PCSC1/BinaryMat/PCA12.",opname,".wdlabel.v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$name2)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

annocol <- data.frame(name=mapping$group1, stringsAsFactors=F)
rownames(annocol) <- colnames(jacc)
jacc2 <- jacc
colnames(jacc2) <- mapping$name2
rownames(jacc2) <- mapping$name2

pdf(paste0("results/PCSC1/BinaryMat/Heatmap.",opname,".pdf"), useDingbats = F)
pheatmap(jacc2,
         scale="none",
         annotation_row = annocol,
         annotation_col = annocol,
         col=rev(bluered299),
         fontsize_row=5,
         fontsize_col=5,
         cluster_rows=T,
         cluster_cols=T,cex=0.5)
dev.off()




#--------------------------------------------------------------------------
## KitchenSink --- REALLY NEED TO MAKE A FUNCTION !!!
#--------------------------------------------------------------------------

## Load binary mat
binarymat <- fread("data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
colnames(binarymat) <- gsub("_peaks", "", colnames(binarymat))
cellNames <- colnames(binarymat[,4:ncol(binarymat)])

mapping <- read.delim("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[cellNames,]
rownames(binarymat) <- paste(binarymat$seqname, binarymat$start, binarymat$end, sep="_")
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,0)

opname="KitchenSink1.Consensus"
jacc <- (getJaccardDist(binarymat[,4:dim(binarymat)[2]]))
pcaPRComp <- prcomp(t(jacc))
df_out <- as.data.frame(pcaPRComp$x)
df_out$group <- mapping$name2


pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$tissue) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v3.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$study) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".stemlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$name2)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

annocol <- data.frame(name=mapping$group1, stringsAsFactors=F)
rownames(annocol) <- colnames(jacc)
jacc2 <- jacc
colnames(jacc2) <- mapping$name2
rownames(jacc2) <- mapping$name2

pdf(paste0("results/KitchenSink1/BinaryMat/Heatmap.",opname,".pdf"), useDingbats = F)
pheatmap(jacc2,
         scale="none",
         annotation_row = annocol,
         annotation_col = annocol,
         col=rev(bluered299),
         fontsize_row=5,
         fontsize_col=5,
         cluster_rows=T,
         cluster_cols=T,cex=0.5)
dev.off()

#--------------------------------------------------------------------------
### Run this only on common promoter and common enhancer regions
#--------------------------------------------------------------------------

## Subset binarymat for common promoter regions

#intersectBed -a data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak -b results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Promotergrp.txt -u > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.Consensus.Catalogue.narrowPeak.intersect.u.Common.Promotergrp.txt
# intersectBed -a data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak -b results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Enhancergrp.txt -u > results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.Consensus.Catalogue.narrowPeak.intersect.u.Common.Enhancergrp.txt


common.prom <- read.table("results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.Consensus.Catalogue.narrowPeak.intersect.u.Common.Promotergrp.txt", header=F, sep="\t", stringsAsFactors=F)
common.enh <- read.table("results/KitchenSink1/Cluster/KCCA.Flexclust/ExtractedClusters/KitchenSink1.Consensus.Catalogue.narrowPeak.intersect.u.Common.Enhancergrp.txt", header=F, sep="\t", stringsAsFactors=F)
rownames(common.prom ) <- paste(common.prom$V1, common.prom$V2, common.prom$V3, sep="_")
rownames(common.enh) <- paste(common.enh$V1, common.enh$V2, common.enh$V3, sep="_")

common.prom.binarymat <- binarymat[rownames(common.prom),]
common.enh.binarymat <- binarymat[rownames(common.enh),]

jacc.prom <- (getJaccardDist(common.prom.binarymat[,4:dim(common.prom.binarymat)[2]]))
jacc.enh <- (getJaccardDist(common.enh.binarymat[,4:dim(common.enh.binarymat)[2]]))

pcaPRComp.prom <- prcomp(t(jacc.prom))
df_out.prom <- as.data.frame(pcaPRComp.prom$x)
df_out.prom$group <- mapping$name2

pcaPRComp.enh <- prcomp(t(jacc.enh))
df_out.enh <- as.data.frame(pcaPRComp.enh$x)
df_out.enh$group <- mapping$name2

opname="common.prom"
pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()


pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
#geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$stem)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
#theme(legend.position="none") +
geom_point()
dev.off()

opname="common.enh"
pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
#geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$stem)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
#theme(legend.position="none") +
geom_point()
dev.off()



#--------------------------------------------------------------------------
### Run this only on common promoter and common enhancer regions --- V2
#--------------------------------------------------------------------------

#  Rscript $SCRIPTDIR/createbinarymat.R results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Enhancergrp.txt data/ConsensusSet/KitchenSink1/temp.peaks ".narrowPeak" results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/ Common.Enhancergrp.KitchenSink1.Binarymat.txt

common.prom <- read.table("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Promotergrp.KitchenSink1.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F)
common.enh <- read.table("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/Common.Enhancergrp.KitchenSink1.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(common.prom ) <- paste(common.prom[,1],common.prom[,2],common.prom[,3],sep="_")
rownames(common.enh) <- paste(common.enh[,1],common.enh[,2],common.enh[,3],sep="_")

jacc.prom <- (getJaccardDist(common.prom[,4:dim(common.prom)[2]]))
jacc.enh <- (getJaccardDist(common.enh[,4:dim(common.enh)[2]]))

pcaPRComp.prom <- prcomp(t(jacc.prom))
df_out.prom <- as.data.frame(pcaPRComp.prom$x)
df_out.prom$group <- mapping$name2

pcaPRComp.enh <- prcomp(t(jacc.enh))
df_out.enh <- as.data.frame(pcaPRComp.enh$x)
df_out.enh$group <- mapping$name2

opname="Common.Promotergrp.KitchenSink1"
pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out.prom,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
#geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$stem)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
#theme(legend.position="none") +
geom_point()
dev.off()

opname="Common.Enhancergrp.KitchenSink1"
pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink1/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out.enh,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
#geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$stem)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
#theme(legend.position="none") +
geom_point()
dev.off()



#--------------------------------------------------------------------------
## KitchenSink2 --- REALLY NEED TO MAKE A FUNCTION !!!
#--------------------------------------------------------------------------

## Load binary mat
binarymat <- fread("data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
colnames(binarymat) <- gsub("_peaks", "", colnames(binarymat))
cellNames <- colnames(binarymat[,4:ncol(binarymat)])

mapping <- read.delim("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[cellNames,]
rownames(binarymat) <- paste(binarymat$seqname, binarymat$start, binarymat$end, sep="_")
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,0)
mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,0))

opname="KitchenSink2.Consensus"
jacc <- (getJaccardDist(binarymat[,4:dim(binarymat)[2]]))
pcaPRComp <- prcomp(t(jacc))
df_out <- as.data.frame(pcaPRComp$x)
df_out$group <- mapping$name2


pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$tissue) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".v3.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$study) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".stemlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$stem) )) +
# geom_text_repel(segment.size  = 0.2,aes( label = mapping$V3)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".wdlabel.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$group1)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".wdlabel.v2.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$group1) )) +
geom_text_repel(size =1,segment.size  = 0.2,aes( label = mapping$name2)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()


pdf(paste0("results/KitchenSink2/BinaryMat/PCA12.",opname,".wdlabel.v3.pdf"), useDingbats = F)
ggplot(df_out,aes(x=PC1,y=PC2,color=factor(mapping$name2) )) +
geom_text_repel(size=2,segment.size  = 0.2,aes( label = mapping$name2)) +
scale_fill_manual(values = colorRampPalette(solarized_pal()(30))(colourCount),guide = guide_legend(nrow=2)) +
theme(legend.position="none") +
geom_point()
dev.off()

annocol <- data.frame(name=mapping$tissue, stem=mapping$stem,stringsAsFactors=T)
rownames(annocol) <- colnames(jacc)
jacc2 <- jacc
#colnames(jacc2) <- mapping$name2
#rownames(jacc2) <- mapping$name2

pdf(paste0("results/KitchenSink2/BinaryMat/Heatmap.",opname,".pdf"), useDingbats = F)
pheatmap(jacc2,
         scale="none",
         clustering_method="ward.D2",
         annotation_row = annocol,
         annotation_col = annocol,
         col=rev(brewer_celsius),
         #fontsize_row=5,
         #fontsize_col=5,
         cluster_rows=T,
         cluster_cols=T,
         cex=0.5)
dev.off()
