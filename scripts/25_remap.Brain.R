#!/bin/bash
## Intersect with Remap for brain clustering

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
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)

#module load bedtools/2.26.0
#files="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Brain.CSC.sp.bed
#data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak
#data/ConsensusSet/PCSC1/Combined.Brain.Consensus.Catalogue.narrowPeak" ;
#
# for f in  $files ;
# do
#     name=$(echo $f |awk -F"/" '{print $NF}' )
#     echo $name $f
#     for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/remap.tf.hg38/*bed ;
#     do
#         echo $k
#         count=$(intersectBed -a $f -b $k -u -f 0.1 | wc -l)
#         tf=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/.nr_macs2_hg38_v1_2.bed//g')
#         echo $tf $count >> results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/remap/Support.$name.txt
#     done
# done
#

######################################
### Run analysis
######################################


## all tfs present in remap database
tflist <- fread("tflist.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)

files <- dir("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/remap/", pattern="Support.")
files <- files[grep(".txt", files)]
combined.enh <- tflist
for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/remap/",files[f]), header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.enh <- merge(combined.enh, a, by.x="V1", by.y="tf", all.x=T)
}

85465 results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/K100.Brain.CSC.sp.bed
1096401 data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak"
541574 data/ConsensusSet/PCSC1/Combined.Brain.Consensus.Catalogue.narrowPeak

enh.sampleuni <- data.frame(name=colnames(combined.enh)[2:ncol(combined.enh)],
                            count=c(541574,85465,1096401))



## Run Fisher test for each TF for each group
functiondt=function(ds, uni,opdir,...){

  namelist=tflist2=pval=or=NULL

  for ( k in c(3)){

    opname=colnames(ds)[k]
    print(opname)

    for (f in 1:nrow(ds)){

      tfname <- as.character(ds[f,"V1"])
      print(tfname)
      hitInSample <- ds[f,k] ##hitInSample
      hitInPop <- ds[f,4]
      sampleSize=(subset(uni$count, uni$name==opname))
      failInPop <- (subset(uni$count, uni$name==colnames(ds)[4]))-hitInPop ## Population Size - Hit in population

      ft <-fisher.test(matrix(c(hitInSample,
                                (hitInPop-hitInSample),
                                (sampleSize-hitInSample),
                                (failInPop-sampleSize +hitInSample)),
                              2, 2), alternative='greater')
      namelist <- c(namelist, opname)
      tflist2 <- c(tflist2, tfname)
      pval <- c(pval,ft$p.value)
      or <- c(or, ft$estimate)

    }
  }

  ft.ds <- data.frame(name=namelist, tf=tflist2, pval=pval, or=or)

  ## Histograms for p value dist for each group
  for (f in names(table(ft.ds$name))){

    dat <- subset(ft.ds, ft.ds$name==f)
    name=gsub(".txt","", f)

    pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/remap/",opdir,".Pvaldist.Greater.Fishertest.",name,".pdf" ))
    plot(hist(dat$pval, nclass=20))
    dev.off()

  }
  ft.ds <- ft.ds[order(-ft.ds$or),]
  ft.ds$qval <- p.adjust(ft.ds$pval, method="BH")
  write.table(ft.ds, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/remap/",opdir,".FisherTest.txt" ), row.names=F, col.names=T, sep="\t", quote=F)

}

functiondt(combined.enh, enh.sampleuni, "Brain")
functiondt(combined.enh, enh.sampleuni, "Brain.KS1")


## Downloaded from remapbluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)

files <- dir("./", pattern="remap.hg38.tsv")

tflist <- read.table("tflist.bed", header=F, sep="\t", stringsAsFactors=F)
tflist$name <- tolower(tflist[,1])
combined <- tflist

for (f  in 1:length(files)){
  a <- read.table(files[f], header=T, sep="\t", stringsAsFactors=F)
  a$fc <- log2(a[,2]/a[,3])
  colnames(a)[4:5] <- c(paste0("Log10Eval.", files[f]), paste0("FC.", files[f]))
  a[,5] <- ifelse(a[,4]>=5, a[,5], 0)
  combined <- merge(combined, a[,c(1,4,5)], by.y="Transcription.Factor", by.x="name", all.x=T)
}

colnames(combined) <- gsub(".remap.hg38.tsv","", colnames(combined))

combined.fc <- combined[,grep("FC", colnames(combined))]
rownames(combined.fc) <- combined$V1

pdf("chk.pdf") ; pheatmap((as.matrix(combined.fc)),
                        clustering_method="complete",
                        clustering_distance_rows="euclidean",
                scale="none",fontsize_row=2) ; dev.off()

combined.fc <- combined.fc[order(combined.fc[,3]),]
write.table(combined.fc, file="Combined.remap.Brain.bed", sep="\t", row.names=T, col.names=T, quote=F)
