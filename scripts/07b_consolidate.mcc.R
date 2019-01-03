### R-3.4.1
### mordor
### Objective : Consolidate MCC
### Not a good way of clustering because it doens't distinguigh shared regions very well when clustered
### NOT RUN


############################################################################
### load dependencies
############################################################################
library(data.table)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(fpc)
hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)

## Followed https://manning-content.s3.amazonaws.com/download/e/dc31390-3cb7-49dd-ab02-937c1af1c2e1/PDSwR_CH08.pdf
############################################################################
## Define functions
############################################################################

sqr_edist <- function(x, y) {
  sum((x-y)^2)
}
wss.cluster <- function(clustermat) {
  c0 <- apply(clustermat, 2, FUN=mean)
  sum(apply(clustermat, 1, FUN=function(row){sqr_edist(row,c0)}))
}
wss.total <- function(dmatrix, labels) {
  wsstot <- 0
  k <- length(unique(labels))
  for(i in 1:k)
    wsstot <- wsstot + wss.cluster(subset(dmatrix, labels==i))
  wsstot
}

totss <- function(dmatrix) {
  grandmean <- apply(dmatrix, 2, FUN=mean)
  sum(apply(dmatrix, 1, FUN=function(row){sqr_edist(row, grandmean)}))
}

ch_criterion <- function(dmatrix, kmax, method="kmeans") {
  if(!(method %in% c("kmeans", "hclust"))) {
    stop("method must be one of c('kmeans', 'hclust')")
  }
  npts <- dim(dmatrix)[1] # number of rows.
  totss <- totss(dmatrix)
  wss <- numeric(kmax)
  crit <- numeric(kmax)
  wss[1] <- (npts-1)*sum(apply(dmatrix, 2, var))
  for(k in 2:kmax) {
    if(method=="kmeans") {
      clustering<-kmeans(dmatrix, k, nstart=10, iter.max=100)
      wss[k] <- clustering$tot.withinss
    }else { # hclust
      d <- dist(dmatrix, method="euclidean")
      pfit <- hclust(d, method="ward")
      labels <- cutree(pfit, k=k)
      wss[k] <- wss.total(dmatrix, labels)
    }
  }
  bss <- totss - wss
  crit.num <- bss/(0:(kmax-1))
  crit.denom <- wss/(npts - 1:kmax)
  list(crit = crit.num/crit.denom, wss = wss, totss = totss)
}


############################################################################
### Read data
############################################################################

lscfiles <- dir("results/PCSC1/MCC/mcc.split/", pattern="MCC.LSC")
gbmfiles <- dir("results/PCSC1/MCC/mcc.split/", pattern="MCC.GBM")
pfafiles <- dir("results/PCSC1/MCC/mcc.split/", pattern="MCC.PFA")

############################################################################
## Consolidate
############################################################################

lsc <- NULL
for ( f in 1:length(lscfiles)){
  
  a <- fread(paste0("results/PCSC1/MCC/",lscfiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  lsc <- rbind(lsc,a)
}


gbm <- NULL
for ( f in 1:length(gbmfiles)){
  
  b <- fread(paste0("results/PCSC1/MCC/",gbmfiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  gbm <- rbind(gbm,b)
}


pfa <- NULL
for ( f in 1:length(pfafiles)){
  
  c <- fread(paste0("results/PCSC1/MCC/",pfafiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  pfa <- rbind(pfa,c)
}


## Combine

c1 <- merge(lsc, gbm, by.x="id", by.y="id")
c2 <- merge(c1, pfa, by.x="id", by.y="id")
write.table(c2, file=paste0("results/PCSC1/MCC/Consolidated.results/Combined.MCC.txt"), sep="\t", row.names=F, col.names=T, quote=F)


############################################################################
## Determine no. of clusters
############################################################################

######################################
## Using within sum of squares
######################################

wss <- sapply(1:15, function(k){kmeans(c2[,2:4], k , n=30)$tot.withinss})

pdf(paste0("results/PCSC1/MCC/Consolidated.results/Kmeans.Elbow.pdf"))
par(mar=c(2,2,2,2))
plot(1:15, wss,
     type="b", pch = 20, frame = FALSE, 
     xlab="Number of clusters K",xlim=c(0,15),xaxt="n",las=2,
     ylab="Total within-clusters sum of squares")
axis(side = 1, at=1:15)
dev.off()

write.table(wss, file=paste0("results/PCSC1/MCC/Consolidated.results/Kmeans.WSS.txt"), row.names=T, col.names=T, sep="\t", quote=F)


######################################
## Using WSS and CH criterion
######################################

clustcrit <- ch_criterion(as.matrix(c2[,2:4]), 10, method="kmeans")
critframe <- data.frame(k=1:10, ch=(clustcrit$crit),wss=(clustcrit$wss))
critframe <- melt(critframe, id.vars=c("k"),variable.name="measure",value.name="score")

pdf(paste0("results/PCSC1/MCC/Consolidated.results/Kmeans.Criterion.pdf"))
ggplot(critframe, aes(x=k, y=score, color=measure)) +
  geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
  scale_x_continuous(breaks=1:10, labels=1:10)
dev.off()

## k=5 looks like the winner here

############################################################################
## automatically run k means over a range of K and choose best K
############################################################################
clustering.ch <- kmeansruns(as.matrix(c2[,2:4]), krange=1:10, criterion="ch") ## CH criterion. 

critframe <- data.frame(k=1:10, ch=scale(clustering.ch$crit))
critframe <- melt(critframe, id.vars=c("k"),variable.name="measure", value.name="score")

pdf("results/PCSC1/MCC/Consolidated.results/CH.criterion.pdf")
ggplot(critframe, aes(x=k, y=score, color=measure)) +
  geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
  scale_x_continuous(breaks=1:10, labels=1:10)
dev.off()

## highest CH is for k=10 but k=5 is also close

############################################################################
## Use bootstrap resampling to evaluate stability of cluster at k=5
############################################################################
kbest.p<-5
cboot<-clusterboot(as.matrix(c2[,2:4]), clustermethod=kmeansCBI,runs=100,iter.max=100,  krange=kbest.p, seed=15555)

## clusters are stable
groups <- cboot$result$partition

kbest.p <-10
cboot10 <-clusterboot(as.matrix(c2[,2:4]), clustermethod=kmeansCBI,runs=100,iter.max=100,  krange=kbest.p, seed=15555)

## clusters are stable
groups <- cboot$result$partition

############################################################################
## Now that we have evaluated that k=5 is stable and can be used, plot heatmap with k=5
############################################################################

## note this doesn't include regions that were common in all...i.e rowsum=48

pdf(paste0("results/PCSC1/MCC/Consolidated.results/PHeatmap.Kmeans5.pdf"))
ph <-pheatmap(as.matrix(cboot$result$result$centers), col=hmcols,scale="none",
              cluster_rows=F, cluster_cols=F); 
dev.off()


pdf(paste0("results/PCSC1/MCC/Consolidated.results/PHeatmap.Kmeans10.pdf"))
ph <-pheatmap(as.matrix(c2[,2:4]), col=hmcols,scale="none",kmeans_k=10,
              cluster_rows=T, cluster_cols=F); 
dev.off()

pdf(paste0("results/PCSC1/MCC/Consolidated.results/PHeatmap.Kmeans10.pdf"))
ph <-pheatmap(as.matrix(c2[,2:4]), col=hmcols,scale="none",
              kmeans_k=10, clustering_distance_rows="correlation",
              clustering_method="ward.D2",
              cluster_rows=T, cluster_cols=F); 
dev.off()


## 2d representation
pdf(paste0("results/PCSC1/MCC/Consolidated.results/PlotCluster.Kmeans5.pdf"))
plotcluster(c2[,2:4],groups)
dev.off()

## Image
c3 <- c2
c3$group <- cboot$result$result$cluster
c3 <- c3[order(c3$group),]
#c3 <- rbind(subset(c2, c2$group==1), subset(c2, c2$group==4), subset(c2, c2$group==2), subset(c2, c2$group==3), subset(c2, c2$group==5) )
pdf(paste0("results/PCSC1/MCC/Consolidated.results/Image.Kmeans5.pdf"))
image(t(apply(as.matrix(c3[,2:4]),2,rev)), col=hmcols)
dev.off()

## Write out regions with its cluster annotations
write.table(c3, file=paste0("results/PCSC1/MCC/Consolidated.results/Combined.MCC.WithK5ClustersAnno.txt"), sep="\t", row.names=F, col.names=T, quote=F)

## Also write back the regions with rowSums ==48 that were excluded in the MCC calculation
enh.mat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F)
write.table(subset(enh.mat[,1:3], rowSums(enh.mat[,4:ncol(enh.mat)])==48), file="results/PCSC1/MCC/Consolidated.results/Ubi.regions.rowsum48.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

