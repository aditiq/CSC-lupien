### R-3.4.1
### mordor
### Objective : combine results of Remap TFBS overlap for scABC peaks


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

files <- dir("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Enhancer/", pattern="Support.")
combined.enh <- tflist
for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Enhancer/",files[f]), header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.enh <- merge(combined.enh, a, by.x="V1", by.y="tf", all.x=T)
}

files <- dir("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Promoter/", pattern="Support.")
combined.prom <- tflist

for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/Promoter/",files[f]), header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.prom <- merge(combined.prom, a, by.x="V1", by.y="tf", all.x=T)
}

enh.sampleuni <- data.frame(name=c("Support.Common.txt","Support.GBM.txt","Support.LSC.txt" ,
                               "Support.PCSC1.Consensus.Catalogue.Enhancers.bed.txt",
                               "Support.PFA.txt","Support.Shared.txt"),
                        count=c(18089,115592,91480,333018,57030,50827))

prom.sampleuni <- data.frame(name=c("Support.Common.txt","Support.GBM.txt","Support.LSC.txt" ,
                                   "Support.PCSC1.Consensus.Catalogue.Promoters.bed.txt",
                                   "Support.PFA.txt","Support.Shared.txt"),
                            count=c(19521,26304,24755,99113,10694,17839))



## Run Fisher test for each TF for each group
functiondt=function(ds, uni,opdir,...){
  
  namelist=tflist2=pval=or=NULL
  
  for ( k in c(2,3,4,6,7)){ 
    
    opname=colnames(ds)[k]
    print(opname)
    
    for (f in 1:nrow(ds)){
      
      tfname <- as.character(ds[f,"V1"])
      print(tfname)
      hitInSample <- ds[f,k] ##hitInSample
      hitInPop <- ds[f,5] 
      sampleSize=(subset(uni$count, uni$name==opname))
      failInPop <- (subset(uni$count, uni$name==colnames(ds)[5]))-hitInPop ## Population Size - Hit in population
     
      ft <-fisher.test(matrix(c(hitInSample, (hitInPop-hitInSample), (sampleSize-hitInSample), (failInPop-sampleSize +hitInSample)),
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
    
    pdf(paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/",opdir,"/Pvaldist.Greater.Fishertest.",name,".pdf" ))
    plot(hist(dat$pval, nclass=20))
    dev.off()
  
  }
  
  ft.ds$enriched <- ifelse(ft.ds$pval <0.01 ,1,0)   
  ft.ds.wide <- dcast(tf ~ name, data = ft.ds, value.var = "enriched")
  rownames(ft.ds.wide) <- ft.ds.wide$tf
  ft.ds.wide$tf <- NULL
  pdf(paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/",opdir,"/Heatmap.pdf" ))
  par(mar=c(5.1, 4.1, 0.5, 0.5))
  heatmap.2(as.matrix(ft.ds.wide[,2:4]), col=c("white","red"), trace="none",scale="none", key=F,cexRow=0.2, cexCol=0.5)
  dev.off()
  
  write.table(ft.ds, file=paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/remap/",opdir,"/FisherTest.txt" ), row.names=F, col.names=T, sep="\t", quote=F)

}

functiondt(combined.enh, enh.sampleuni, "Enhancer")
functiondt(combined.prom, prom.sampleuni, "Promoter")

