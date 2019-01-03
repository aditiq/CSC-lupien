### R-3.4.1
### mordor
### NOT RUN -- NOT CORRECT TO RUN K-MEDOIDS
### Objective : Run K-medoids clustering on sampling of the consensus set

args = commandArgs(trailingOnly=TRUE)
ds <- as.character(args[1])
opname <- as.character(args[2])

############################################################################
### load dependencies
############################################################################
library(prabclus)
library(factoextra)
library(fpc)
library(ggplot2)
library(reshape2)

binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")

## Run k-medoid clustering on each sampling and determine no. of clusters using gap statistic

enh <- read.table(ds, header=F, sep="\t", stringsAsFactors=F)
rownames(enh) <- paste(enh[,1],enh[,2],enh[,3],sep="_")

matenh <- as.matrix(binarymat[rownames(enh),4:ncol(binarymat)])

## compute distance matrix
dat1 <- prabclus::jaccard(t(as.matrix(matenh)))

## K-medoids to identify best K based on average silhouette width

pamx.asw <- pamk(dat1,krange=2:100,criterion="asw",seed=14355,scaling=F,diss=TRUE, critout=T,ns=10,usepam=TRUE)
save(pamx.asw, file=paste0("results/PCSC1/Cluster/Sampling/PamK.asw.", opname, ".Rdata"))

pamx.ch <- pamk(dat1,krange=2:100,criterion="ch",seed=14355,scaling=F,diss=TRUE, critout=T,ns=10,usepam=TRUE)
save(pamx.ch, file=paste0("results/PCSC1/Cluster/Sampling/PamK.ch.", opname, ".Rdata"))

## Visualise
library(ggplot2)

for ( f in 0:99){
  
  load(paste0("results/PCSC1/Cluster/Sampling/PamK.asw.Enh.sample.",f,".Rdata"))
  load(paste0("results/PCSC1/Cluster/Sampling/PamK.ch.Enh.sample.",f,".Rdata"))
  
  print(pamx.ch$nc)
  print(pamx.asw$nc)
  
  critframe <- data.frame(k=1:100, ch=scale(pamx.ch$crit),asw=scale(pamx.asw$crit))
  critframe <- melt(critframe, id.vars=c("k"),
                    variable.name="measure",
                    value.name="score")
  
  #pdf(paste0("results/PCSC1/Cluster/Sampling/PamK.Sample",f,".pdf"))
  p <- ggplot(critframe, aes(x=k, y=score, color=measure)) +
    geom_point(aes(shape=measure)) + geom_line(aes(linetype=measure)) +
    scale_x_continuous(breaks=1:100, labels=1:100)
  #print(p)
  #dev.off()
}