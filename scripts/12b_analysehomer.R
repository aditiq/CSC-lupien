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

analysehomer=function(ipdir,searchstring,...){

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
    if ( nrow(motifs.l[[i]]) ==0 ) {
      namescol <- colnames(motifs.l[[i]])
      motifs.l[[i]] <-   (rbind(rep("nothing", dim(motifs.l[[i]])[2]), motifs.l[[i]]))
      colnames(motifs.l[[i]]) <- namescol
      for (f in 1:ncol( motifs.l[[i]]) ){

          motifs.l[[i]][,f] <- as.character(motifs.l[[i]][,f])

      }
    }
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
  mots.df<-Reduce(function(x,y) merge(x,y,all=TRUE),mots.l)
  write.table(mots.df,file=paste(ipdir,"Motif.Summary",searchstring,".txt", sep=""),sep="\t",row.names=F)

  #create heatmap
  data<-mots.df[,2:ncol(mots.df)]
  data <- data[, grep(".FE", colnames(data))]
  data[!is.na(data)]<-1
  data[is.na(data)]<- 0
  data<-as.matrix(data)
  mode(data)<-"numeric"
  row.names(data)<-mots.df[,1]

  rownames(data) <- str_split_fixed(rownames(data), "/",2)[,1]
  rownames(data) <- str_split_fixed(rownames(data), "\\(",2)[,1]

  colnames(data) <- gsub("FoldPercEnr_","",  gsub(gsub("/","", ipdir),"",colnames(data)))

  data2 <- data[,colSums(data)>0]
  write.table(data2,file=paste(ipdir,"Motif.HeatmapData",searchstring,".txt", sep=""),sep="\t",row.names=T, col.names=T, quote=F)

  pdf(file=paste(ipdir,"Heatmap.",searchstring,".pdf",sep=""),height=12,width=12, useDingbats = F)
  pheatmap((as.matrix(data2)), col=c("white","red"), clustering_method="ward.D2",
           cluster_rows=T, cluster_cols=T,  scale="none", fontsize_row=5)
  dev.off()
}

# ==============================================================================
# Run analysis
# ==============================================================================

#analysehomer(ipdir="results/PCSC1/Cluster/scABC/motif/", searchstring="Enhancer")
#analysehomer(ipdir="results/PCSC1/Cluster/scABC/motif/", searchstring="Promoter")
analysehomer(ipdir="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg", searchstring="Brain.bg")
analysehomer(ipdir="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/homer/Brain.bg", searchstring="Brain.bg")

### Cicero
for ( name1 in c("common", "\\.shared", "notshared"))
{
  for (name in c("PCSC1","KitchenSink1","KitchenSink2"))
  {
      print(paste0(name, "-",name1))
      analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.promoters.", name))
      #analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.enhancer.", name))
      analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.5.promoters.", name))
      #analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.5.enhancer.", name))
      analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.2.promoters.", name))
      #analysehomer(ipdir="results/PCSC1/cicero/run2.readcount/homer/", searchstring=paste0(name1, ".coaccessgt0.2.enhancer.", name))
    }
}


for ( name1 in c("common", "\\.shared", "notshared"))
{
  for (name in c("PCSC1","KitchenSink1","KitchenSink2"))
  {
      print(paste0(name, "-",name1))
      analysehomer(ipdir="results/PCSC1/C3D/finalresults/homer/",
                    searchstring=paste0(name1, ".", name,".Consensus.Catalogue.narrowPeak"))
    }
}

plothm2 = function(filename,opname,...){
  mots.df <- read.table(filename, header=T, sep="\t", stringsAsFactors=F)
  data<-mots.df[,2:ncol(mots.df)]
  data <- data[, grep(".bg.FE", colnames(data))]
  row.names(data)<-mots.df[,1]
  colnames(data) <- gsub("FoldPercEnr_","",  gsub(gsub("/","", "results/PCSC1/Cluster/scABC/motif/"),"",colnames(data)))
  data[,paste0("Else.", opname,".bg.FE")] <- NULL
  colnames(data) <- gsub(paste0(".",opname,".bg.FE"),"",colnames(data))
  data[is.na(data)] <- 0
  data<-as.matrix(data)
  mode(data)<-"numeric"
  data2 <- subset(data,rowSums(data)>0)
  rownames(data2) <- str_split_fixed(rownames(data2), "/",2)[,1]
  rownames(data2) <- str_split_fixed(rownames(data2), "\\(",2)[,1]
  colnames(data2) <- gsub(opname,"",colnames(data2))
  colnames(data2) <- gsub("\\..",".",colnames(data2))

  col2 <- brewer.pal(n = 8, name = "Reds")

  data2 <- data2[order(data2[,4],data2[,1],data2[,3],data2[,2] ),]
  pdf(paste0("results/PCSC1/Cluster/scABC/motif/",opname,".Heatmap.v2.pdf" ))
  pheatmap(as.matrix(data2[,c(1,3,2,4)]),
           cluster_cols = F,cluster_rows = T,
           clustering_method="ward.D2",
           clustering_distance_rows = "euclidean",
           color=scalered,
           trace="none",scale="none")
  dev.off()
}

plothm2(filename='results/PCSC1/Cluster/scABC/motif/Motif.SummaryPromoter.txt', opname="Promoter")
plothm2(filename='results/PCSC1/Cluster/scABC/motif/Motif.SummaryEnhancer.txt', opname="Enhancer")
