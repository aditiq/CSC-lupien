### R-3.4.1
### mordor
### Objective : Combine mapbed into signal matrices

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

#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------
## LSCP ##
#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------

#----------------------------------------------------------------
## Load data
#----------------------------------------------------------------

samples<-list.files(path = "data/ConsensusSet/PCSC1/mapbed/lscp/",recursive=TRUE,pattern="LSCp.")
samples<-gsub(".Consensus.mapbed.SignalMat.txt","",basename(samples))

#create list
map.l<-list()
for (i in 1:length(samples)) {
  map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/mapbed/lscp/",samples[i],".Consensus.mapbed.SignalMat.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
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

write.table(map.df,"data/ConsensusSet/PCSC1/LSCp.Consensus.Mapbed.txt",col.names=T,sep="\t", row.names=F, quote=F)

#----------------------------------------------------------------
## Quantile normalise
#----------------------------------------------------------------

map.df.qnorm <- cbind(map.df[,1:3],normalize.quantiles(as.matrix(map.df[,4:ncol(map.df)])))
colnames(map.df.qnorm) <- colnames(map.df)
write.table(map.df.qnorm,"data/ConsensusSet/PCSC1/LSCp.Consensus.Mapbed.Qnorm.txt",col.names=T,sep="\t", row.names=F, quote=F)


#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------
## PFA ##
#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------

#----------------------------------------------------------------
## Load data
#----------------------------------------------------------------

samples<-list.files(path = "data/ConsensusSet/PCSC1/mapbed/pfa/",recursive=TRUE,pattern="PFA.")
samples<-gsub(".Consensus.mapbed.SignalMat.txt","",basename(samples))

#create list
map.l<-list()
for (i in 1:length(samples)) {
  map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/mapbed/pfa/",samples[i],".Consensus.mapbed.SignalMat.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
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

write.table(map.df,"data/ConsensusSet/PCSC1/GBM.Consensus.Mapbed.txt",col.names=T,sep="\t", row.names=F, quote=F)

#----------------------------------------------------------------
## Quantile normalise
#----------------------------------------------------------------

map.df.qnorm <- cbind(map.df[,1:3],normalize.quantiles(as.matrix(map.df[,4:ncol(map.df)])))
colnames(map.df.qnorm) <- colnames(map.df)
write.table(map.df.qnorm,"data/ConsensusSet/PCSC1/PFA.Consensus.Mapbed.Qnorm.txt",col.names=T,sep="\t", row.names=F, quote=F)



#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------
## GBM ##
#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------

#----------------------------------------------------------------
## Load data
#----------------------------------------------------------------

samples<-list.files(path = "data/ConsensusSet/PCSC1/mapbed/gbm/",recursive=TRUE,pattern="GBM.")
samples<-gsub(".Consensus.mapbed.SignalMat.txt","",basename(samples))

#create list
map.l<-list()
for (i in 1:length(samples)) {
  map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/mapbed/gbm/",samples[i],".Consensus.mapbed.SignalMat.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
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

write.table(map.df,"data/ConsensusSet/PCSC1/GBM.Consensus.Mapbed.txt",col.names=T,sep="\t", row.names=F, quote=F)

#----------------------------------------------------------------
## Quantile normalise
#----------------------------------------------------------------

map.df.qnorm <- cbind(map.df[,1:3],normalize.quantiles(as.matrix(map.df[,4:ncol(map.df)])))
colnames(map.df.qnorm) <- colnames(map.df)
write.table(map.df.qnorm,"data/ConsensusSet/PCSC1/GBM.Consensus.Mapbed.Qnorm.txt",col.names=T,sep="\t", row.names=F, quote=F)




#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------
## PCSC1 ##
#----------------------------------------------------------------#----------------------------------------------------------------#----------------------------------------------------------------

#----------------------------------------------------------------
## Load data
#----------------------------------------------------------------

samples<-list.files(path = "data/ConsensusSet/PCSC1/mapbed/pcsc1/",recursive=TRUE,pattern="PCSC1.")
samples<-gsub(".Consensus.mapbed.SignalMat.txt","",basename(samples))

#create list
map.l<-list()
for (i in 1:length(samples)) {
  map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/mapbed/pcsc1/",samples[i],".Consensus.mapbed.SignalMat.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
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

write.table(map.df,"data/ConsensusSet/PCSC1/PCSC1.Consensus.Mapbed.txt",col.names=T,sep="\t", row.names=F, quote=F)

#----------------------------------------------------------------
## Quantile normalise
#----------------------------------------------------------------

map.df.qnorm <- cbind(map.df[,1:3],normalize.quantiles(as.matrix(map.df[,4:ncol(map.df)])))
colnames(map.df.qnorm) <- colnames(map.df)
write.table(map.df.qnorm,"data/ConsensusSet/PCSC1/PCSC1.Consensus.Mapbed.Qnorm.txt",col.names=T,sep="\t", row.names=F, quote=F)
