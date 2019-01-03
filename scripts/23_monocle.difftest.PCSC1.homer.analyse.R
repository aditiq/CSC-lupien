### R-3.4.1
### mordor
### Objective : parse homer results

# ==============================================================================
# load dependencies
# ==============================================================================
library(gplots)
library(RColorBrewer)
library(pheatmap)
require(ggplot2)
library(reshape)
library(stringr)


hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(100)

# ==============================================================================
# Make function
# ==============================================================================

ipdir="results/PCSC1/monocle.difftest/homer/"
searchstring="beta"


samples<-list.files(ipdir,  full.names=T,recursive=T, pattern="knownResults.txt" )
samples <- samples[grep(searchstring,samples)]
#create motif list, keeping only those with qval < 0.001
# and fold enrichment in target seqs vs bg seqs of >=1.5
motifs.l<-list()
number_motifs<-vector()
for (i in 1:length(samples)) {
  motifs.l[[i]]<-read.delim(paste(samples[i],sep=""),sep="\t",header=T,stringsAsFactors = F)
  motifs.l[[i]]<-motifs.l[[i]][which(p.adjust(motifs.l[[i]][,3],method="BH")<=0.001),]
  motifs.l[[i]][,1]<-gsub("/Homer","",motifs.l[[i]][,1])
  motifs.l[[i]][,1]<-gsub("/HOMER","",motifs.l[[i]][,1])
  #motifs.l[[i]][,1]<-gsub("/.*","",motifs.l[[i]][,1])
  #motifs.l[[i]][,1]<-gsub("\\(.*","",motifs.l[[i]][,1])
  motifs.l[[i]][,7]<-gsub("%.*","",motifs.l[[i]][,7])
  motifs.l[[i]][,9]<-gsub("%.*","",motifs.l[[i]][,9])
  motifs.l[[i]]$NewCol<-as.numeric(motifs.l[[i]][,7])/as.numeric(motifs.l[[i]][,9])
  names(motifs.l[[i]])<-c(names(motifs.l[[i]])[1:9],paste("FoldPercEnr_",samples[i],sep=""))
  motifs.l[[i]]<-motifs.l[[i]][which(motifs.l[[i]][,10]>=1.5),]
  #motifs.l[[i]]<-motifs.l[[i]][which(as.numeric(gsub("%","",motifs.l[[i]][,7]))>=10),]
  number_motifs<-c(number_motifs,nrow(motifs.l[[i]]))
}
names(motifs.l)<-samples


#simplify the list keeping only motif name, fold difference between target and background and pval

mots.l<-list()
number_mots<-vector()
for (i in 1:length(samples)) {
  mots.l[[i]]<-motifs.l[[i]][,c(1,3,10)]
  name <- gsub("/","",gsub("FoldPercEnr_results/motif//","",
                           gsub("knownResults.txt","",
                                gsub(gsub("/","", ipdir),"",colnames( mots.l[[i]][3])))))

  colnames( mots.l[[i]]) <- c( "MotifName", paste0(name, ".NegLog10Pval"), paste0(name, ".FE") )
}
names(mots.l) <- gsub("results/motif//","",dirname(samples))


#merge into a data.frame
mots.df<-Reduce(function(x,y) merge (x,y,all=TRUE),mots.l)
write.table(mots.df,file=paste(ipdir,"Motif.Summary",searchstring,".txt", sep=""),sep="\t",row.names=F)

#create heatmap
data<-mots.df[,2:ncol(mots.df)]
data <- data[, grep(".bg.FE", colnames(data))]
data[!is.na(data)]<-1
data[is.na(data)]<- 0
data<-as.matrix(data)
mode(data)<-"numeric"
row.names(data)<-mots.df[,1]

rownames(data) <- str_split_fixed(rownames(data), "/",2)[,1]
rownames(data) <- str_split_fixed(rownames(data), "\\(",2)[,1]

colnames(data) <- gsub("betaFE","",gsub(".betapatch.qval.0.05.bg.","",gsub("FoldPercEnr_","",  gsub(gsub("/","", ipdir),"",colnames(data)))))

## remove ESC
data2 <- data[,-grep("ESC", colnames(data))]
colnames(data2) <- gsub("negative", ".negative", gsub("positive",".positive",colnames(data2)))
data2 <- data2[,colSums(data2)>0]
data2 <- data2[rowSums(data2)>0,]

pdf(file=paste(ipdir,"Heatmap.",searchstring,".pdf",sep=""),height=12,width=12, useDingbats = F)
pheatmap((as.matrix(data2)), col=c("white","red"),
          clustering_method="ward.D2",
          clustering_distance_rows="euclidean",
          #cellwidth=12,
          #cellheight=12,
         cluster_rows=T, cluster_cols=T,  scale="none", fontsize_row=5)
dev.off()

## Remove LSCp vs LSC neg and LSCp vs LSCb
data3 <- data2[, -grep("LSCposandneg", colnames(data2))]
data3 <- data3[, -grep("LSCp.LSCb", colnames(data3))]

pdf(file=paste(ipdir,"Heatmap.",searchstring,".v2.pdf",sep=""),height=12,width=12, useDingbats = F)
pheatmap((as.matrix(data3)), col=c("white","red"),
          clustering_method="ward.D2",
          clustering_distance_rows="euclidean",
          #cellwidth=12,
          #cellheight=12,
         cluster_rows=T, cluster_cols=T,  scale="none", fontsize_row=5)
dev.off()
