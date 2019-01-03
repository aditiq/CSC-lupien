### R-3.4.1
### mordor
### Objective : Exploratory analysis of lsc+ and LSC- to check which behave abnormally and exclude them from downstream analysis

#----------------------------------------------------------------
### load dependencies
#----------------------------------------------------------------

library(data.table)
library("RColorBrewer")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
library(gplots)
library(pheatmap)
library(vioplot)
library(ggplot2)
library(matrixStats)
library(preprocessCore)
library(ggrepel)


#----------------------------------------------------------------
## Load data
#----------------------------------------------------------------

samples<-list.files(path = "data/ConsensusSet/PCSC1/mapbed",recursive=TRUE,pattern="All.LSC.Consensus.mapbed")
samples<-gsub(".All.LSC.Consensus.mapbed.SignalMat.txt","",basename(samples))

#create list
map.l<-list()
for (i in 1:length(samples)) {
  map.l[[i]]<-fread(paste("data/ConsensusSet/PCSC1/mapbed/LSCn.LSCp/",samples[i],".All.LSC.Consensus.mapbed.SignalMat.txt",sep=""),sep="\t",header=F, data.table=F,skip=1)
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

write.table(map.df,"data/ConsensusSet/PCSC1/LSCp.LSCn.Consensus.Mapbed.txt",col.names=T,sep="\t", row.names=F, quote=F)

#----------------------------------------------------------------
## Quantile normalise
#----------------------------------------------------------------

map.df.qnorm <- cbind(map.df[,1:3],normalize.quantiles(as.matrix(map.df[,4:ncol(map.df)])))
colnames(map.df.qnorm) <- colnames(map.df)

#----------------------------------------------------------------
## Process
#----------------------------------------------------------------

## PCA Clustering
ds <- map.df.qnorm[,4:ncol(map.df.qnorm)]
PCA<-prcomp(t(as.matrix(ds)))
PCA_dim <- data.frame(x = PCA$x[,1], y = PCA$x[,2],  z = PCA$x[,3])
PCA_eig <- (PCA$sdev)[1:3]^2/sum(PCA$sdev^2)*100
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

df_out <- as.data.frame(PCA$x)
df_out$group <- substr(rownames(df_out),1,4)
df_out$name <- substr(rownames(df_out),6,nchar(rownames(df_out))-15 )


q <- ggplot(df_out,aes(x=PC1,y=PC2,color=group ))+
  geom_point(size=3,alpha=0.5,show.legend=T) +
  scale_color_manual(values = brewer.pal(7,"Dark2")) +
  geom_text_repel(aes(PC1, PC2, label = df_out$name), size = 2) +
  theme(legend.position = "none") +
  xlab( paste( "PC1 (" , round(PCA_eig[1],2) , "%", ")", sep="") )  +
  ylab( paste( "PC2 (", round(PCA_eig[2],2) ,  "%", ")", sep="") )

pdf(paste0("results/PCSC1/LSCp.LSCn.exploratory/PCA12.qnorm.pdf"),useDingbats=F) ; print(q) ; dev.off()

p <- ggplot(df_out,aes(x=PC2,y=PC3,color=group )) +
  geom_point(size=3,alpha=0.5,show.legend=T) +
  scale_color_manual(values = brewer.pal(7,"Dark2")) +
  geom_text_repel(aes(PC2, PC3, label = df_out$name), size = 2) +
  theme(legend.position = "none") +
  xlab( paste( "PC2 (" , round(PCA_eig[2],2) , "%", ")", sep="") )  +
  ylab( paste( "PC3 (", round(PCA_eig[3],2) ,  "%", ")", sep="") )

pdf(paste0("results/PCSC1/LSCp.LSCn.exploratory/PCA23.qnorm.pdf"),useDingbats=F) ; print(p) ; dev.off()

## Pearson Correlation

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
cor <- cor(as.matrix(ds))
annotation=data.frame(stringsAsFactors = T, name=substr(colnames(ds),1,4))
annotation$number <- seq(1,ncol(ds),1)
rownames(annotation) <- colnames(ds)

## Add NPM1c information
mut=c("8081_TAAGGCGA",
      "8082_CGTACTAG",
      "8083_TCCTGAGC",
      "90239",
      "8096_TCCTGAGC",
      "8097_TAAGGCGA",
      "8098_CGTACTAG",
      "8099_TCCTGAGC",
      "100091",
      "8090_TCCTGAGC",
      "8026_AGGCAGAA",
      "8027_AGGCAGAA",
      "8028_AGGCAGAA",
      "8029_TAAGGCGA",
      "8030_CGTACTAG",
      "8101_TAAGGCGA",
      "8102_CGTACTAG",
      "8103_TCCTGAGC",
      "8104_TAAGGCGA",
      "8105_CGTACTAG",
      "8066_TCCTGAGC",
      "8067_TAAGGCGA",
      "8068_CGTACTAG",
      "8069_TCCTGAGC",
      "8070_CGTACTAG",
      "8026_AGGCAGAA",
      "8027_AGGCAGAA",
      "8028_AGGCAGAA",
      "8029_TAAGGCGA",
      "8030_CGTACTAG")

annotation$NPM1c <- "FALSE"
for (f in 1:length(mut)) { annotation[grep(mut[f], rownames(annotation)),"NPM1c"] <- "TRUE" }
annotation$NPM1c <- ifelse(annotation$name=="LSCp", "FALSE",annotation$NPM1c)
annotation$NPM1c <- as.factor(annotation$NPM1c)

pdf(paste0("results/PCSC1/LSCp.LSCn.exploratory/Pearson.qnorm.pdf")) ;
pheatmap(cor,cluster_rows=T, cluster_cols=T, col=colors,cex=0.7,
         annotation_row=annotation[,c(1,3)],
         annotation_col=annotation[,c(1,3)])
dev.off()


### Remove LSC- that are NPM1c
ds.subset <- ds[,-c(subset(annotation$number, annotation$name=="LSCn" & annotation$NPM1c=="TRUE"))]
annotation.subset <- annotation[colnames(ds.subset),]
cor2 <- cor(as.matrix(ds.subset))
pdf(paste0("results/PCSC1/LSCp.LSCn.exploratory/Pearson.qnorm.wdtLSCn.NPM1c.pdf")) ;
pheatmap(cor2,cluster_rows=T, cluster_cols=T, col=colors,cex=0.7, annotation_row=annotation.subset[,c(1,3)],annotation_col=annotation.subset[,c(1,3)])
dev.off()


## Top variant rows
ds$var <- apply(ds,1,var)
ds2 <- subset(ds, ds$var >= quantile(ds$var, 0.95))
ds2$var <- NULL

pdf("results/PCSC1/LSCp.LSCn.exploratory/Heatmap.qnorm.pdf")
pheatmap(as.matrix(ds2),
         color = bluered299,
         cluster_rows =T,
         clustering_distance_rows="euclidean",
         cluster_cols =T,
         clustering_distance_cols="euclidean",
         scale = "row",clustering_method="ward.D2",
         show_rownames = F,
         annotation_col=annotation[,c(1,3)])
dev.off()
