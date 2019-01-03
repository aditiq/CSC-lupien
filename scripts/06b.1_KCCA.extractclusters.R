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

bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract clusters for Enhancers
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------
# Load binary matrix
#-----------------------------------
enhmat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt",check.names=F, stringsAsFactors = F, header=T, sep="\t")
rownames(enhmat) <- paste(enhmat[,1], enhmat[,2], enhmat[,3], sep="_")

#-----------------------------------
# Load KCCA dataset
#-----------------------------------
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Enhancer.100.Rdata")
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Enhancer.rownames.100.Rdata")

flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
rownames(flexclust.clus) <- rownames.obj

cc <- as.matrix(kcca.cl$centers)
colnames(cc) <- colnames(enhmat)[4:ncol(enhmat)]
rownames(cc) <- seq(1, max(kcca.cl$cluster),1)

#-----------------------------------
# Assign groups
#-----------------------------------
common <- c(92)
shared <- c(10,74,95,20,6,36,79,38,77,72,8,40,60,94,91,30,28,39,53,57,61,67,49,21,3,69,78,7,71,87,90)
lsc.sp <- c(12,83,68,27,32,31,19,85,18,14,55,46,48,4,22,59,93,100,88,56,70)
gbm.sp <- c(9,29,41,42,45,50,52,73,75,82,63,15,89,81,34,2,11,16,47,24,23,54,17,51,13,86,66,62,64,76,80,96,84,98,99)
pfa.sp <- c(25,97,44,43,58,65,37,33,35, 26, 1,5)

#-----------------------------------
# Plot image
#-----------------------------------
pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Enhancer.Flexclust.pdf"))
heatmap.2(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),
          col=hmcols, scale="none",
          Rowv=NULL, Colv=NULL,
          trace="none",
          RowSideColors=c(rep("#F1D302",length(common)),rep("pink",length(shared)),rep("#C1292E",length(lsc.sp)),
                          rep("#8CAE68",length(gbm.sp)),rep("#235789",length(pfa.sp))),
          cexRow=0.5,
          add.expr = abline(v=c(18, 41))
)
dev.off()

pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Enhancer.Flexclust.Image.pdf"))
image(t(apply(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),2,rev)) ,col=scalered)
dev.off()

pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Enhancer.Flexclust.Image.v2.pdf"))
image(t(apply(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),2,rev)) ,col=c("white","red"))
dev.off()


#-----------------------------------
# Write the groups out
#-----------------------------------
flexclust.clus$id <- rownames(flexclust.clus)
flexclust.clus$FlexClust.group <- ifelse(flexclust.clus[,1] %in% common, "Common",
                                         ifelse(flexclust.clus[,1] %in% shared, "Shared",
                                                ifelse(flexclust.clus[,1] %in% lsc.sp, "LSC",
                                                       ifelse(flexclust.clus[,1] %in% gbm.sp, "GBM",
                                                              ifelse(flexclust.clus[,1] %in% pfa.sp, "PFA",NA)))))


colnames(flexclust.clus) <- c("Flexclust.ClusterNo","id","FlexClust.group")
write.table(flexclust.clus[,c("id","Flexclust.ClusterNo","FlexClust.group")],
            file="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/FlexClust.Enhancer.Groups.txt", row.names=F, col.names=T, sep="\t", quote=F)


for ( f in names(table(flexclust.clus$FlexClust.group))) {

  print(f)
  dat <- as.matrix(enhmat[subset(flexclust.clus$id, flexclust.clus$FlexClust.group==f),4:ncol(enhmat)])
  x <- (1:nrow(dat))
  y <- (1:ncol(dat))
  pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/",f,".Enhancer.Flexclust.Image.pdf"))
  image(y, x, t(dat), col=c("white","red"), axes=FALSE,xlab="",ylab="",srt=45)
  axis(3, at = 1:ncol(dat), labels=colnames(dat),srt=45,tick=FALSE)
  axis(2, at = 1:nrow(dat), labels=rownames(dat),srt=45,tick=FALSE)
  abline(v=c(18,41))
  dev.off()

}

for ( f in names(table(flexclust.clus$FlexClust.group))[2:5]) {

  print(f)
  pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/",f,".Enhancer.Flexclust.heatmap.pdf"))
  heatmap.2(as.matrix(cc[unique(subset(flexclust.clus$Flexclust.ClusterNo, flexclust.clus$FlexClust.group==f)),]),
            col=hmcols, scale="none",
            trace="none",cexRow=0.5,Colv=NULL,
            add.expr = abline(v=c(18, 41)),useRaster=TRUE,
            hclustfun=function(x) hclust(x,method="ward.D2"),
            distfun=function(x) as.dist((1 - cor(  t(x), method="spearman"  ))))
  dev.off()
}


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract clusters for Promoters
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------
# Load binary matrix
#-----------------------------------
promat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoter.Binarymat.txt",check.names=F, stringsAsFactors = F, header=T, sep="\t")
rownames(promat) <- paste(promat[,1], promat[,2], promat[,3], sep="_")

#-----------------------------------
# Load KCCA dataset
#-----------------------------------
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Promoter.100.Rdata")
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Promoter.rownames.100.Rdata")

flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
rownames(flexclust.clus) <- rownames.obj

cc <- as.matrix(kcca.cl$centers)
colnames(cc) <- colnames(promat)[4:ncol(promat)]
rownames(cc) <- seq(1, max(kcca.cl$cluster),1)

## Plot image unsuper
pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Promoter.Flexclust.unsuper.pdf"))
pheatmap(as.matrix(cc),
          col=hmcols, scale="none",fontsize_row=3.5,
          cluster_rows=TRUE, cluster_cols=FALSE,
          clustering_method="average",
          clustering_distance_rows = "correlation",
          trace="none")
dev.off()

#-----------------------------------
# Assign groups
#-----------------------------------
common<-c(87,86,81,79,58,17)
shared <-c(94,89,88,84,83,82,80,75,73,72,71,65,64,54,53,52,51,44, 41,30,29,27,26,25,20,19,14,5,3)
lsc.sp <-c(100,97,93,90,85,78,74,69,66,61, 57, 50, 46, 40, 38,28,23,21,16,11,9,8,6,4,2)
gbm.sp <-c(99,98,95,92,76,70,68,67,63,62,60,59, 56,55,49,48,47,45,43,42 , 36,35,34,33,32,24,22,18,15, 10,7,1)
pfa.sp <-c(96,91,77,39,37,31, 13,12)

## 22,35

#-----------------------------------
# Plot image
#-----------------------------------
pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Promoter.Flexclust.pdf"))
heatmap.2(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),
          col=hmcols, scale="none",
          Rowv=NULL, Colv=NULL,
          trace="none",
          RowSideColors=c(rep("#F1D302",length(common)),rep("pink",length(shared)),rep("#C1292E",length(lsc.sp)),
                          rep("#8CAE68",length(gbm.sp)),rep("#235789",length(pfa.sp))),
          cexRow=0.5,
          add.expr = abline(v=c(18, 41))
)
dev.off()

pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Promoter.Flexclust.Image.pdf"))
image(t(apply(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),2,rev)) ,col=scalered)
dev.off()

pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Promoter.Flexclust.Image.v2.pdf"))
image(t(apply(as.matrix(cc[c(common, shared, lsc.sp, gbm.sp, pfa.sp),]),2,rev)) ,col=c("white","red"))
dev.off()


#-----------------------------------
# Write the groups out
#-----------------------------------
flexclust.clus$id <- rownames(flexclust.clus)
flexclust.clus$FlexClust.group <- ifelse(flexclust.clus[,1] %in% common, "Common",
                                         ifelse(flexclust.clus[,1] %in% shared, "Shared",
                                                ifelse(flexclust.clus[,1] %in% lsc.sp, "LSC",
                                                       ifelse(flexclust.clus[,1] %in% gbm.sp, "GBM",
                                                              ifelse(flexclust.clus[,1] %in% pfa.sp, "PFA",NA)))))


colnames(flexclust.clus) <- c("Flexclust.ClusterNo","id","FlexClust.group")
write.table(flexclust.clus[,c("id","Flexclust.ClusterNo","FlexClust.group")],
            file="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/FlexClust.Promoter.Groups.txt", row.names=F, col.names=T, sep="\t", quote=F)


for ( f in names(table(flexclust.clus$FlexClust.group))) {

  print(f)
  dat <- as.matrix(promat[subset(flexclust.clus$id, flexclust.clus$FlexClust.group==f),4:ncol(promat)])
  x <- (1:nrow(dat))
  y <- (1:ncol(dat))
  pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/",f,".Promoter.Flexclust.Image.pdf"))
  image(y, x, t(dat), col=c("white","red"), axes=FALSE,xlab="",ylab="",srt=45)
  axis(3, at = 1:ncol(dat), labels=colnames(dat),srt=45,tick=FALSE)
  axis(2, at = 1:nrow(dat), labels=rownames(dat),srt=45,tick=FALSE)
  abline(v=c(18,41))
  dev.off()

}

for ( f in names(table(flexclust.clus$FlexClust.group))[1:5]) {

  print(f)
  pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/",f,".Promoter.Flexclust.heatmap.pdf"))
  heatmap.2(as.matrix(cc[unique(subset(flexclust.clus$Flexclust.ClusterNo, flexclust.clus$FlexClust.group==f)),]),
            col=hmcols, scale="none",
            trace="none",cexRow=0.5,Colv=NULL,
            add.expr = abline(v=c(18, 41)),useRaster=TRUE,
            hclustfun=function(x) hclust(x,method="average"),
            distfun=function(x) as.dist((1 - cor(  t(x), method="spearman"  ))))
  dev.off()
}


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract clusters for CSC.ESC
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------
# Load binary matrix
#-----------------------------------
cscescmat <- read.table("data/ConsensusSet/PCSC1/CSC.ESC.Consensus.Catalogue.Binarymat.txt",check.names=F, stringsAsFactors = F, header=T, sep="\t")
rownames(cscescmat) <- paste(cscescmat[,1], cscescmat[,2], cscescmat[,3], sep="_")

#-----------------------------------
# Load KCCA dataset
#-----------------------------------
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.CSC.ESC.100.Rdata")
load("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.CSC.ESC.rownames.100.Rdata")

flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
rownames(flexclust.clus) <- rownames.obj

cc <- as.matrix(kcca.cl$centers)
colnames(cc) <- colnames(cscescmat)[4:ncol(cscescmat)]
rownames(cc) <- seq(1, max(kcca.cl$cluster),1)

## Plot image unsuper
pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.CSC.ESC.Flexclust.unsuper.pdf"))
pheatmap(as.matrix(cc),
          col=hmcols, scale="none",fontsize_row=3.5,
          cluster_rows=TRUE, cluster_cols=TRUE,
          clustering_method="average",
          clustering_distance_rows = "correlation",
          trace="none")
dev.off()
