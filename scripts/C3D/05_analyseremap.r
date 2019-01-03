### R-3.4.1
### mordor
### Objective : combine results of Remap TFBS overlap for cicero


######################################
### load dependencies
######################################
library(gplots)
library(RColorBrewer)
library(pheatmap)
require(ggplot2)
library(reshape)
library(stringr)
library(data.table)
library(qvalue)

######################################
### Run analysis
######################################


## all tfs present in remap database
tflist <- fread("tflist.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)

files <- dir("results/PCSC1/C3D/finalresults/remap//", pattern="Support.")
files <- files[grep("nhancer", files)]

combined.enh <- tflist
for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/remap/",files[f]),
  header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.enh <- merge(combined.enh, a, by.x="V1", by.y="tf", all.x=T)
}
combined.enh  <- combined.enh[!duplicated(combined.enh),]

files2 <- dir("results/PCSC1/cicero/run2.readcount/remap//", pattern="Support.")
files2 <- files2[grep("romoter", files2)]

combined.prom <- tflist
for ( f in 1:length(files2)){
  a <- read.table(paste0("results/PCSC1/cicero/run2.readcount/remap/",files2[f]),
  header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files2[f])
  combined.prom <- merge(combined.prom, a, by.x="V1", by.y="tf", all.x=T)
}


enh.sampleuni <- data.frame(name=c("Support.PCSC1.Consensus.Catalogue.Enhancers.bed.txt",
                                    "Support.KitchenSink1.Consensus.Catalogue.Enhancers.bed.txt",
                                    "Support.KitchenSink2.Consensus.Catalogue.Enhancers.bed.txt"),
                            count=c(333018,939483,780152))

prom.sampleuni <- data.frame(name=c("Support.PCSC1.Consensus.Catalogue.Promoters.bed.txt",
                                    "Support.KitchenSink1.Consensus.Catalogue.Promoters.bed.txt",
                                    "Support.KitchenSink2.Consensus.Catalogue.Promoters.bed.txt"),
                                    count=c(70111,156918,133340))

samplesizelist <- read.delim("results/PCSC1/cicero/run2.readcount/remap/samplesize.txt",
                              header=F, sep="\t", stringsAsFactors=F )

## Run Fisher test for each TF for each group
functiondt=function(ds,colsofinterest,popcol,uni, name,...){

  namelist=tflist2=pval=or=NULL

    for ( k in colsofinterest){

      opname=colnames(ds)[k]
      print(opname)

      for (f in 1:nrow(ds)){

        tfname <- as.character(ds[f,"V1"])
        print(tfname)
        print(f)
        hitInSample <- ds[f,k] ##hitInSample
        hitInPop <- ds[f,popcol]
        sampleSize=(subset(samplesizelist$V1, samplesizelist$V2==gsub("Support.","",gsub(".txt",".bed",opname))))
        failInPop <- (subset(uni$count, uni$name==colnames(ds)[popcol]))-hitInPop ## Population Size - Hit in population

        if((hitInPop-hitInSample)<=0){
          ft$p.value=1
            ft$estimate=0
        } else {
        ft <-fisher.test(matrix(c(hitInSample, (hitInPop-hitInSample),
                                (sampleSize-hitInSample),
                                (failInPop-sampleSize +hitInSample)),
                                2, 2), alternative='greater')
        }
        namelist <- c(namelist, opname)
        tflist2 <- c(tflist2, tfname)
        pval <- c(pval,ft$p.value)
        or <- c(or, ft$estimate)
      }
      }

  ft.ds <- data.frame(name=namelist, tf=tflist2, pval=pval, or=or)

  # ## Histograms for p value dist for each group
  # for (f in names(table(ft.ds$name))){
  #
  #   dat <- subset(ft.ds, ft.ds$name==f)
  #   name=gsub(".txt","", f)
  #
  #   pdf(paste0("results/PCSC1/cicero/run2.readcount/remap/Pvaldist.Greater.Fishertest.",name,".pdf" ))
  #   plot(hist(dat$pval, nclass=20))
  #   dev.off()
  #
  # }
  ft.ds$enriched <- ifelse(ft.ds$pval <0.001  & ft.ds$or > 5,1,0)
  ft.ds2 <- ft.ds[,c("tf","name","enriched")]

  ft.ds.wide <- dcast(tf ~ name, data =   ft.ds2)
  rownames(ft.ds.wide) <- ft.ds.wide$tf
  ft.ds.wide$tf <- NULL
  write.table(ft.ds, file=paste0("results/PCSC1/cicero/run2.readcount/remap/",name,".FisherTest.txt" ), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(ft.ds.wide,
          file=paste0("results/PCSC1/cicero/run2.readcount/remap/Heatmap.",name,".FisherTest.txt" ),
            row.names=T, col.names=T, sep="\t", quote=F)

  # pdf(paste0("results/PCSC1/cicero/run2.readcount/remap/Heatmap.",name,".pdf" ))
  # par(mar=c(5.1, 4.1, 5, 5))
  # heatmap.2(as.matrix(ft.ds.wide), col=c("white","red"),
  #           trace="none",scale="none", key=F,cexRow=0.2, cexCol=0.2)
  # dev.off()
}

functiondt(combined.enh, c(2,3,4,5,6,7,8,9,10), 11, enh.sampleuni, "GBM.KS1.Enhancer")
functiondt(combined.enh, c(2,3,4,5,6,7,8,9,10), 12, enh.sampleuni, "GBM.KS2.Enhancer")
functiondt(combined.enh, c(2,3,4,5,6,7,8,9,10), 22, enh.sampleuni, "GBM.PCSC1.Enhancer")

functiondt(combined.enh, c(13,14,15,16,17,18,19,20,21), 11, enh.sampleuni, "LSCp.KS1.Enhancer")
functiondt(combined.enh, c(13,14,15,16,17,18,19,20,21), 12, enh.sampleuni, "LSCp.KS2.Enhancer")
functiondt(combined.enh, c(13,14,15,16,17,18,19,20,21), 22, enh.sampleuni, "LSCp.PCSC1.Enhancer")

functiondt(combined.enh, c(23,24,25,26,27,28,29,30,31), 11, enh.sampleuni, "PFA.KS1.Enhancer")
functiondt(combined.enh, c(23,24,25,26,27,28,29,30,31), 12, enh.sampleuni, "PFA.KS2.Enhancer")
functiondt(combined.enh, c(23,24,25,26,27,28,29,30,31), 22, enh.sampleuni, "PFA.PCSC1.Enhancer")


functiondt(combined.prom, c(2,3,4,5,6,7,8,9,10), 11, prom.sampleuni, "GBM.KS1.promoter")
functiondt(combined.prom, c(2,3,4,5,6,7,8,9,10), 12, prom.sampleuni, "GBM.KS2.promoter")
functiondt(combined.prom, c(2,3,4,5,6,7,8,9,10), 22, prom.sampleuni, "GBM.PCSC1.promoter")

functiondt(combined.prom, c(13,14,15,16,17,18,19,20,21), 11, prom.sampleuni, "LSCp.KS1.promoter")
functiondt(combined.prom, c(13,14,15,16,17,18,19,20,21), 12, prom.sampleuni, "LSCp.KS2.promoter")
functiondt(combined.prom, c(13,14,15,16,17,18,19,20,21), 22, prom.sampleuni, "LSCp.PCSC1.promoter")

functiondt(combined.prom, c(23,24,25,26,27,28,29,30,31), 11, prom.sampleuni, "PFA.KS1.promoter")
functiondt(combined.prom, c(23,24,25,26,27,28,29,30,31), 12, prom.sampleuni, "PFA.KS2.promoter")
functiondt(combined.prom, c(23,24,25,26,27,28,29,30,31), 22, prom.sampleuni, "PFA.PCSC1.promoter")


## Summarise
hmfiles <- dir("results/PCSC1/cicero/run2.readcount/remap/", pattern="Heatmap")

hmfiles.e.pcsc1 <- hmfiles[grep("PCSC1.Enhancer",hmfiles)]

tf.gt0.2 <- c()
tf.gt0.5 <- c()
tf.gt0 <- c()

for (f in 1:length(hmfiles.e.pcsc1)){
  a <- read.delim(paste0("results/PCSC1/cicero/run2.readcount/remap/",hmfiles.e.pcsc1[f]),sep="\t",header=T,stringsAsFactors = F)
  tf.gt0.2 <- c(tf.gt0.2,subset(rownames(a), a[,1]==1 & a[,4]==0 & a[,7]==1))
  tf.gt0.5 <- c(tf.gt0.5,subset(rownames(a), a[,2]==1 & a[,5]==0 & a[,8]==1))
  tf.gt0 <- c(tf.gt0,subset(rownames(a), a[,3]==1 & a[,6]==0 & a[,9]==1))
}

df.tf.gt0.2 <- as.data.frame(table(tf.gt0.2 ))
df.tf.gt0.5 <- as.data.frame(table(tf.gt0.5 ))
df.tf.gt0 <- as.data.frame(table(tf.gt0 ))

df.tf.gt0 <- subset(df.tf.gt0, df.tf.gt0$Freq==3)
df.tf.gt0.2 <- subset(df.tf.gt0.2, df.tf.gt0.2$Freq==3)
df.tf.gt0.5 <- subset(df.tf.gt0.5, df.tf.gt0.5$Freq==3)
