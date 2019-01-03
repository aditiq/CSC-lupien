### R-3.4.1
### mordor
### Objective : Combine coveragebed into readcount matrices

#----------------------------------------------------------------
### load dependencies
#----------------------------------------------------------------

library(data.table)
library("RColorBrewer")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
library(gplots)
library(pheatmap)
library(ggplot2)
library(matrixStats)
library(preprocessCore)
library(ggrepel)


combinfunc=function(dir1, pattern1, ...){

  #----------------------------------------------------------------
  ## Load data
  #----------------------------------------------------------------

  samples<-list.files(path = paste0("data/ConsensusSet/PCSC1/readcounts/", dir1, "/"),recursive=TRUE,pattern=pattern1)
  samples<-gsub(".Consensus.coveragebed.ReadCount.txt","",basename(samples))

  #create list
  map.l<-list()
  for (i in 1:length(samples)) {
    map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/readcounts/",dir1,"/",samples[i],".Consensus.coveragebed.ReadCount.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
    colnames(map.l[[i]])<-c("chr","start","end",as.character(samples[i]))
  }
  names(map.l)<-samples

  #----------------------------------------------------------------
  #put them together in a dataframe
  #----------------------------------------------------------------

  map.df<-map.l[[1]][,1:3]
  for (i in 1:length(samples)) {
    map.df<-cbind(map.df,map.l[[i]][,5])
  }
  colnames(map.df)<-c("chr","start","end",as.character(samples))

  write.table(map.df,paste0("data/ConsensusSet/PCSC1/", pattern1, "Consensus.ReadCount.txt"),col.names=T,sep="\t", row.names=F, quote=F)

}

combinfunc("lscp", "LSCp.")
combinfunc("gbm", "GBM.")
combinfunc("pfa", "PFA.")
