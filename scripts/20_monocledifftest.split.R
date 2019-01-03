
library(stringr)
library(pheatmap)
library(data.table)


files <- dir("results/PCSC1/monocle.difftest/", pattern="difftest.txt", full.names=T)
for ( f in 1:length(files)){

  a <- read.table(files[f], header=T,sep="\t", stringsAsFactors=F)
  a2 <- subset(a$Peak, a$qval <0.05)

  a3 <- as.data.frame(str_split_fixed(a2, "_",3))
  write.table(a3, file=paste0("results/PCSC1/monocle.difftest/" , gsub(".difftest.txt","", basename(files[f]) ),".qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)
}

## beta patch -- PCSC1

files <- dir("results/PCSC1/monocle.difftest/", pattern=".difftest.betapatch.txt", full.names=T)
for ( f in 1:length(files)){

  name <- gsub(".difftest.betapatch.txt","",basename(files[f]))
  a <- read.table(files[f], header=T,sep="\t", stringsAsFactors=F)
  a2 <- a[order(a$beta),]
  a3 <- subset(a2$Peak, a2$qval <0.05)

  bmat <- fread(paste0("data/ConsensusSet/PCSC1/",name,".Consensus.Catalogue.Binarymat.txt"), data.table=F,sep="\t", stringsAsFactors=F, header=T, check.names=F)
  colnames(bmat) <- gsub("_peaks","", colnames(bmat))
  rownames(bmat) <- paste(bmat[,1], bmat[,2],bmat[,3], sep="_")
  mat <- as.matrix(bmat[a3,4:ncol(bmat)])

  write.table(as.data.frame(str_split_fixed(a3, "_",3)), file=paste0("results/PCSC1/monocle.difftest/" , name,".betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(as.data.frame(str_split_fixed(subset(a2$Peak, a2$qval <0.05 & a2$beta >0 ), "_",3)), file=paste0("results/PCSC1/monocle.difftest/" , name,"positivebeta.betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(as.data.frame(str_split_fixed(subset(a2$Peak, a2$qval <0.05 & a2$beta <0 ), "_",3)), file=paste0("results/PCSC1/monocle.difftest/" , name,"negativebeta.betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)

  groupinfo <- read.table(paste0("results/PCSC1/monocle.difftest/",name,".anno.txt"), sep="\t", header=F, stringsAsFactors=F)
  rownames(groupinfo) <- groupinfo[,1]
  groupinfo <- groupinfo[order(groupinfo$V2),]
  mat <- mat[,groupinfo$V1]

  ## plot image
  pdf(paste0("results/PCSC1/monocle.difftest/" , gsub(".difftest.betapatch.txt","", basename(files[f])),".Heatmap.pdf"), useDingbats=F)
  pheatmap(mat,
          cluster_rows=F, cluster_cols=F,
          scale="none")
  dev.off()

  negbeta <- subset(a2$Peak, a2$qval <0.05 & a2$beta <0 )
  posbeta <- subset(a2$Peak, a2$qval <0.05 & a2$beta >0 )

  annorow <- data.frame(stringsAsFactors=F, group=c(rep("neg",length(negbeta)), rep("pos",length(posbeta))))
  rownames(annorow) <- c(negbeta,posbeta)

  print(f)
  png(paste0("results/PCSC1/monocle.difftest/" , gsub(".difftest.betapatch.txt","", basename(files[f])),".Heatmap.png"))
  pheatmap(mat[  rownames(annorow) ,],
          cluster_rows=F, cluster_cols=F,
          annotation_row=annorow,
          scale="none")
  dev.off()

}


## beta patch -- KitchenSink2


files <- dir("results/KitchenSink2/monocle.difftest/", pattern=".difftest.betapatch.txt", full.names=T)
for ( f in 1:length(files)){

  name <- gsub(".difftest.betapatch.txt","",basename(files[f]))
  a <- read.table(files[f], header=T,sep="\t", stringsAsFactors=F)
  a2 <- a[order(a$beta),]
  a3 <- subset(a2$Peak, a2$qval <0.05)

  bmat <- fread(paste0("data/ConsensusSet/KitchenSink2/","KitchenSink2",".Consensus.Catalogue.Binarymat.txt"), data.table=F,sep="\t", stringsAsFactors=F, header=T, check.names=F)
  colnames(bmat) <- gsub("_peaks","", colnames(bmat))
  rownames(bmat) <- paste(bmat[,1], bmat[,2],bmat[,3], sep="_")
  mat <- as.matrix(bmat[a3,4:ncol(bmat)])

  write.table(as.data.frame(str_split_fixed(a3, "_",3)), file=paste0("results/KitchenSink2/monocle.difftest/" , name,".betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(as.data.frame(str_split_fixed(subset(a2$Peak, a2$qval <0.05 & a2$beta >0 ), "_",3)), file=paste0("results/KitchenSink2/monocle.difftest/" , name,"positivebeta.betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(as.data.frame(str_split_fixed(subset(a2$Peak, a2$qval <0.05 & a2$beta <0 ), "_",3)), file=paste0("results/KitchenSink2/monocle.difftest/" , name,"negativebeta.betapatch.qval.0.05.bed"), sep="\t", quote=F, row.names=F, col.names=F)

  groupinfo <- read.table(paste0("results/KitchenSink2/monocle.difftest/",name,".anno.txt"), sep="\t", header=F, stringsAsFactors=F)
  rownames(groupinfo) <- groupinfo[,1]
  groupinfo <- groupinfo[order(groupinfo$V2),]
  mat <- mat[,groupinfo$V1]

  ## plot image
  pdf(paste0("results/KitchenSink2/monocle.difftest/" , gsub(".difftest.betapatch.txt","", basename(files[f])),".Heatmap.pdf"), useDingbats=F)
  pheatmap(mat,
          cluster_rows=F, cluster_cols=F,
          scale="none")
  dev.off()

  negbeta <- subset(a2$Peak, a2$qval <0.05 & a2$beta <0 )
  posbeta <- subset(a2$Peak, a2$qval <0.05 & a2$beta >0 )

  annorow <- data.frame(stringsAsFactors=F, group=c(rep("neg",length(negbeta)), rep("pos",length(posbeta))))
  rownames(annorow) <- c(negbeta,posbeta)

  mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
  rownames(mapping) <- mapping$sample
  mapping <- mapping[colnames(mat),]
  mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2)
  mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2))

  annocol <- data.frame(stringsAsFactors=F, stem=mapping$stem, tissue=mapping$tissue)
  rownames(annocol) <- colnames(mat)

  print(f)
  png(paste0("results/KitchenSink2/monocle.difftest/" , gsub(".difftest.betapatch.txt","", basename(files[f])),".Heatmap.posbeta.png"))
  pheatmap(mat[posbeta,],
          cluster_rows=F, cluster_cols=F,
          annotation_col=annocol,
          annotation_names_row=FALSE,
          annotation_names_col=FALSE,
          scale="none")
  dev.off()

  png(paste0("results/KitchenSink2/monocle.difftest/" , gsub(".difftest.betapatch.txt","", basename(files[f])),".Heatmap.negbeta.png"))
  pheatmap(mat[negbeta,],
    cluster_rows=F, cluster_cols=F,
    annotation_col=annocol,
    annotation_names_row=FALSE,
    annotation_names_col=FALSE,
          scale="none")
  dev.off()

}
