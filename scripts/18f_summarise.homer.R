# Objective : Identify TFs that are shared, common and not shared between CSCs

library(data.table)
library(reshape2)
library(ggplot2)

scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)

motiffiles <- dir("results/PCSC1/cicero/run2.readcount/homer/", pattern="Motif.HeatmapData")
common <- motiffiles[grep("common", motiffiles)]
notshared <- motiffiles[grep("notshared", motiffiles)]
shared <- motiffiles[grep("\\.shared", motiffiles)]


for ( f in notshared){
  tryCatch({
    name <- gsub(".txt", "", gsub("Motif.HeatmapDatanotshared.", "", f))
    print(name)
    a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",f), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a) <- gsub(paste0(name, ".FE"), "", colnames(a))
    colnames(a)[1] <- "TF"
    a$sum <- apply(a[,-1],1,sum)
    print(range(a$sum))
    tf <- subset(a$TF, a$sum>=2)
    write.table(tf, file=paste0("results/PCSC1/cicero/run2.readcount/homer/", "TFs.notshared.", name, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
    }, error=function(e){})
}


for ( f in common){
  tryCatch({
    name <- gsub(".txt", "", gsub("Motif.HeatmapDatacommon.", "", f))
    print(name)
    a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",f), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a) <- gsub(paste0(name, ".FE"), "", colnames(a))
    colnames(a)[1] <- "TF"
    a$sum <- apply(a[,-1],1,sum)
    print(range(a$sum))
    tf <- subset(a$TF, a$sum ==3)
    write.table(tf, file=paste0("results/PCSC1/cicero/run2.readcount/homer/", "TFs.common.", name, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
    }, error=function(e){})
}

## common but not in notshared

for ( f in common){
  tryCatch({
    name <- gsub(".txt", "", gsub("Motif.HeatmapDatacommon.", "", f))
    print(name)
    a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",f), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a) <- gsub(paste0(name, ".FE"), "", colnames(a))
    colnames(a)[1] <- "TF"

    name2 <- gsub("Motif.HeatmapDatacommon." , "Motif.HeatmapDatanotshared.", f)
    a2 <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",name2), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a2)[1] <- "TF"

    gbm <- subset(a2$TF, a2$GBM.notshared.==1)
    pfa <- subset(a2$TF, a2$PFA.notshared.==1)
    lsc <- subset(a2$TF, a2$LSCp.notshared.==1)

    a$gbm2 <- ifelse(a$TF %in% gbm, 1, 0 )
    a$pfa2 <- ifelse(a$TF %in% pfa, 1, 0 )
    a$lsc2 <- ifelse(a$TF %in% lsc, 1, 0 )

    a$gbm3 <- ifelse(a$GBM.common. ==1 & a$gbm2 ==1,0, a$GBM.common. )
    a$pfa3 <- ifelse(a$PFA.common. ==1 & a$pfa2 ==1,0, a$PFA.common. )
    a$lsc3 <- ifelse(a$LSCp.common. ==1 & a$lsc2 ==1,0, a$LSCp.common. )

    a$sum <- apply(a[,c("gbm3", "pfa3","lsc3")],1,sum)
    print(range(a$sum))
    tf <- subset(a$TF, a$sum == 3)
    write.table(tf, file=paste0("results/PCSC1/cicero/run2.readcount/homer/", "TFs.common.butnotinnotshared.", name, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
    }, error=function(e){})
}


for ( f in shared){
  tryCatch({
    name <- gsub(".txt", "", gsub("Motif.HeatmapData\\.shared.", "", f))
    print(name)
    a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",f), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a) <- gsub(paste0(name, ".FE"), "", colnames(a))
    colnames(a)[1] <- "TF"
    a$sum <- apply(a[,-1],1,sum)
    tf <- subset(a$TF, a$sum >=2)
    write.table(tf, file=paste0("results/PCSC1/cicero/run2.readcount/homer/", "TFs.shared.", name, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  }, error=function(e){})
}

# shared not in notshared

for ( f in shared){
  tryCatch({
    name <- gsub(".txt", "", gsub("Motif.HeatmapData\\\\.shared", "", f))
    print(name)
    a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",f), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a) <- gsub(paste0(name, ".FE"), "", colnames(a))
    colnames(a)[1] <- "TF"

    name2 <- gsub("\\.shared", "notshared",gsub("\\\\", "", f))
    a2 <- read.table(paste0("results/PCSC1/cicero/run2.readcount/homer/",name2), sep="\t", header=T, stringsAsFactors=F,  row.names=NULL)
    colnames(a2)[1] <- "TF"

    colnames(a2) <- gsub(".FE","",gsub(gsub(".txt", "", gsub("Motif.HeatmapData\\\\.shared", "", f)),"", colnames(a2)))
    gbm <- subset(a2$TF, a2$GBM.notshared==1)
    pfa <- subset(a2$TF, a2$PFA.notshared==1)
    lsc <- subset(a2$TF, a2$LSCp.notshared==1)

    a$gbm2 <- ifelse(a$TF %in% gbm, 1, 0 )
    a$pfa2 <- ifelse(a$TF %in% pfa, 1, 0 )
    a$lsc2 <- ifelse(a$TF %in% lsc, 1, 0 )

    a$gbm3 <- ifelse(a$GBM.shared ==1 & a$gbm2 ==1,0, a$GBM.shared )
    a$pfa3 <- ifelse(a$PFA.shared ==1 & a$pfa2 ==1,0, a$PFA.shared )
    a$lsc3 <- ifelse(a$LSCp.shared ==1 & a$lsc2 ==1,0, a$LSCp.shared )

    a$sum <- apply(a[,c("gbm3", "pfa3","lsc3")],1,sum)
    print(range(a$sum))
    tf <- subset(a$TF, a$sum >=2)
    write.table(tf, file=paste0("results/PCSC1/cicero/run2.readcount/homer/", "TFs.shared.butnotinnotshared.", name, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
    }, error=function(e){})
}


## Plot identified TFs

tfs <- dir("results/PCSC1/cicero/run2.readcount/homer/", pattern="TFs.", full.names=T)
tfs <- tfs[-grep("butnotinnotshared", tfs)]

for ( threshold in c(0, 0.2, 0.5)){

  tflistall <- tfs[grep(paste0("gt", threshold,".enhancer.KitchenSink1.txt"), tfs)]
  tflist <- list()

  for (f in 1:length(tflistall)){
    tflist[[f]] <- read.table(tflistall[f], header=F, sep="\t", stringsAsFactors=F)
    tflist[[f]]$flag <- 1
    colnames(tflist[[f]]) <- c("tf", basename(tflistall[f]))
  }

  names(tflist) <- basename(tflistall[f])
  tf.df<-Reduce(function(x,y) merge(x,y,all=TRUE),tflist)
  tf.df[is.na(tf.df)] <- 0
  colnames(tf.df) <- gsub(paste0(".coaccessgt", threshold),"", gsub("TFs.","",gsub(".Motif.HeatmapData\\\\.shared", "", gsub(".enhancer.KitchenSink1.txt", "", colnames(tf.df)))))
  rownames(tf.df) <- tf.df[,1]
  tf.df.long <- melt(tf.df[,c(1,2,4,3)])
  p <- ggplot(data = tf.df.long, mapping = aes(x = variable, y = tf,fill = value)) +
  scale_fill_gradient(low = 'white', high = 'red') +
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.title.x= element_blank(), text = element_text(size=16), legend.position="none") +
  geom_tile()

  pdf(paste0("results/PCSC1/cicero/run2.readcount/homer/Heatmap.TF.", paste0("gt", threshold,".enhancer.KitchenSink1.pdf")))
  print(p)
  dev.off()
}

for ( threshold in c(0, 0.2, 0.5)){

  tflistall <- tfs[grep(paste0("gt", threshold,".promoters.KitchenSink1.txt"), tfs)]
  tflist <- list()

  for (f in 1:length(tflistall)){
    tflist[[f]] <- read.table(tflistall[f], header=F, sep="\t", stringsAsFactors=F)
    tflist[[f]]$flag <- 1
    colnames(tflist[[f]]) <- c("tf", basename(tflistall[f]))
  }

  names(tflist) <- basename(tflistall[f])
  tf.df<-Reduce(function(x,y) merge(x,y,all=TRUE),tflist)
  tf.df[is.na(tf.df)] <- 0
  colnames(tf.df) <- gsub(paste0(".coaccessgt", threshold),"", gsub("TFs.","",gsub(".Motif.HeatmapData\\\\.shared", "", gsub(".promoters.KitchenSink1.txt", "", colnames(tf.df)))))
  tf.df <- tf.df[!duplicated(tf.df),]
  rownames(tf.df) <- tf.df[,1]
  tf.df <- subset(tf.df, !(tf.df$tf=="nothing"))
  tf.df.long <- melt(tf.df[,c(1,2,4,3)])
  p <- ggplot(data = tf.df.long, mapping = aes(x = variable, y = tf,fill = value)) +
  scale_fill_gradient(low = 'white', high = 'red') +
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.title.x= element_blank(), text = element_text(size=16), legend.position="none",) +
  geom_tile() +
  coord_fixed()


  pdf(paste0("results/PCSC1/cicero/run2.readcount/homer/Heatmap.TF.", paste0("gt", threshold,".promoters.KitchenSink1.pdf")))
  print(p)
  dev.off()
}
