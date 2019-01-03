### R-3.4.1
### mordor
### Objective : Consolidate Jaccard distance
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
### Read data
############################################################################

lscfiles <- dir("results/PCSC1/Cluster/Jaccard/jaccard.split/", pattern="jaccard.LSC")
gbmfiles <- dir("results/PCSC1/Cluster/Jaccard/jaccard.split/", pattern="jaccard.GBM")
pfafiles <- dir("results/PCSC1/Cluster/Jaccard/jaccard.split/", pattern="jaccard.PFA")

############################################################################
## Consolidate
############################################################################

lsc <- NULL
for ( f in 1:length(lscfiles)){
  print(f)
  a <- fread(paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/",lscfiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  lsc <- rbind(lsc,a)
}


gbm <- NULL
for ( f in 1:length(gbmfiles)){
  print(f)
  b <- fread(paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/",gbmfiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  gbm <- rbind(gbm,b)
}


pfa <- NULL
for ( f in 1:length(pfafiles)){
  print(f)
  c <- fread(paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/",pfafiles[f]), header=T, sep="\t", stringsAsFactors = F, data.table=F)
  pfa <- rbind(pfa,c)
}


## Combine

c1 <- merge(lsc, gbm, by.x="id", by.y="id")
c2 <- merge(c1, pfa, by.x="id", by.y="id")
write.table(c2, file=paste0("results/PCSC1/Cluster/Jaccard/Consolidated.results/Combined.jaccard.txt"), sep="\t", row.names=F, col.names=T, quote=F)

#### REMEMBER THIS IS A DISSIMILARITY MATRIX. Higher values mean LOWER JACCARD COEFFICIENT 

########################################################################################################################################################
## Determine no. of clusters
########################################################################################################################################################

############################################################################
## Using within sum of squares -- Kmeans elbow
############################################################################

wss <- sapply(1:15, function(k){kmeans(c2[,2:4], k , n=30)$tot.withinss})

pdf(paste0("results/PCSC1/Cluster/Jaccard/Consolidated.results/Kmeans.Elbow.pdf"))
par(mar=c(2,2,2,2))
plot(1:15, wss,
     type="b", pch = 20, frame = FALSE, 
     xlab="Number of clusters K",xlim=c(0,15),xaxt="n",las=2,
     ylab="Total within-clusters sum of squares")
axis(side = 1, at=1:15)
dev.off()

write.table(wss, file=paste0("results/PCSC1/Cluster/Jaccard/Consolidated.results/Kmeans.WSS.txt"), row.names=T, col.names=T, sep="\t", quote=F)


############################################################################
## automatically run k means over a range of K and choose best K
############################################################################
clustering.ch <- kmeansruns(as.matrix(c2[,2:4]), krange=1:10, criterion="ch") ## CH criterion. 

critframe <- data.frame(k=1:10, ch=scale(clustering.ch$crit))
critframe <- melt(critframe, id.vars=c("k"),variable.name="measure", value.name="score")

pdf("results/PCSC1/Cluster/Jaccard/Consolidated.results/CH.criterion.pdf")
ggplot(critframe, aes(x=k, y=score, color=measure)) +
  geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
  scale_x_continuous(breaks=1:10, labels=1:10)
dev.off()

## highest CH is for k=5

############################################################################
## Use bootstrap resampling to evaluate stability of cluster at k=5
############################################################################
 
kbest.p<-5
cboot<-clusterboot(as.matrix(c2[,2:4]), clustermethod=kmeansCBI,runs=100,iter.max=100,  krange=kbest.p, seed=15555)
 
## clusters are stable
groups <- cboot$result$partition
