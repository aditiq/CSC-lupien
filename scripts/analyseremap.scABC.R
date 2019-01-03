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
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)

######################################
### Run analysis
######################################


## all tfs present in remap database
tflist <- fread("tflist.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)

files <- dir("results/PCSC1/Cluster/scABC/remap/Enhancer/", pattern="Support.")
files <- files[grep(".txt", files)]
combined.enh <- tflist
for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/Cluster/scABC/remap/Enhancer/",files[f]), header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.enh <- merge(combined.enh, a, by.x="V1", by.y="tf", all.x=T)
}

files <- dir("results/PCSC1/Cluster/scABC/remap/Promoter/", pattern="Support.")
files <- files[grep(".txt", files)]

combined.prom <- tflist

for ( f in 1:length(files)){
  a <- read.table(paste0("results/PCSC1/Cluster/scABC/remap/Promoter/",files[f]), header=F, sep=" ", stringsAsFactors = F)
  colnames(a) <- c("tf",files[f])
  combined.prom <- merge(combined.prom, a, by.x="V1", by.y="tf", all.x=T)
}


# 215165 Else.Enhancer.p0.05.bed
# 43351 Else.Promoter.p0.05.bed
# 9248 GBM.Enhancer.p0.05.bed
# 1708 GBM.Promoter.p0.05.bed
# 49804 LSC.Enhancer.p0.05.bed
# 9075 LSC.Promoter.p0.05.bed
# 32548 PFA.Enhancer.p0.05.bed
# 4254 PFA.Promoter.p0.05.bed
# 55255 Shared.Enhancer.p0.05.bed
# 11723 Shared.Promoter.p0.05.bed

enh.sampleuni <- data.frame(name=colnames(combined.enh)[2:ncol(combined.enh)],
                            count=c(215165,9248,49804,362021,32548,55255))

prom.sampleuni <- data.frame(name=colnames(combined.prom)[2:ncol(combined.prom)],
                             count=c(43351,1708,9075,70112,4254,11723))



## Run Fisher test for each TF for each group
functiondt=function(ds, uni,opdir,...){
  
  namelist=tflist2=pval=or=NULL
  
  for ( k in c(2:4,6:7)){ 
    
    opname=colnames(ds)[k]
    print(opname)
    
    for (f in 1:nrow(ds)){
      
      tfname <- as.character(ds[f,"V1"])
      print(tfname)
      hitInSample <- ds[f,k] ##hitInSample
      hitInPop <- ds[f,5] 
      sampleSize=(subset(uni$count, uni$name==opname))
      failInPop <- (subset(uni$count, uni$name==colnames(ds)[5]))-hitInPop ## Population Size - Hit in population
      
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
    
    pdf(paste0("results/PCSC1/Cluster/scABC/remap/",opdir,"/Pvaldist.Greater.Fishertest.",name,".pdf" ))
    plot(hist(dat$pval, nclass=20))
    dev.off()
    
  }
  
  ft.ds$enriched <- ifelse(ft.ds$pval <0.00001 & ft.ds$or > 2,1,0)   
  ft.ds.wide <- dcast(tf ~ name, data = ft.ds, value.var = "enriched")
  rownames(ft.ds.wide) <- ft.ds.wide$tf
  ft.ds.wide$tf <- NULL
  ft.ds.wide <- subset(ft.ds.wide,rowSums(ft.ds.wide)>0)
  pdf(paste0("results/PCSC1/Cluster/scABC/remap/",opdir,"/Heatmap.pdf" ))
  pheatmap(as.matrix(ft.ds.wide),clustering_method="ward.D2",
           clustering_distance_rows = "correlation",
           color=hmcols, fontsize_row=0.8,
            trace="none",scale="none")
  dev.off()
  
  write.table(ft.ds, file=paste0("results/PCSC1/Cluster/scABC/remap/",opdir,"/FisherTest.txt" ), row.names=F, col.names=T, sep="\t", quote=F)
  
}

functiondt(combined.enh, enh.sampleuni, "Enhancer")
functiondt(combined.prom, prom.sampleuni, "Promoter")

## Custom functions
library(tidyr)

plot2=function(filename, opdir,fc,...){
  
  promfish <- read.table(filename, header=T, sep="\t", stringsAsFactors=F)
  promfish <- promfish[order(-promfish$or),]
  promfish <- subset(promfish, promfish$name!='Support.Else.txt')
  promfish$or2 <- ifelse(promfish$pval <0.01 & promfish$or >fc,promfish$or,0)
  promfish.wide <- spread(promfish[,c("name" , "tf" ,"or2")], name, or2)
  promfish.wide <- subset(promfish.wide, rowSums(promfish.wide[,2:5])>0)
  rownames(promfish.wide)<- promfish.wide$tf
  pdf(paste0("results/PCSC1/Cluster/scABC/remap/",opdir,"/Heatmap.v2.pdf" ))
  pheatmap(as.matrix(promfish.wide[,c(2,4,3,5)]),
           cluster_cols = F,
           clustering_method="ward.D2",
           clustering_distance_rows = "euclidean",
           color=scalered, 
           trace="none",scale="none")
  dev.off()
  
}

plot2(filename="results/PCSC1/Cluster/scABC/remap/Promoter/FisherTest.txt", opdir="Promoter",fc=2)
plot2(filename="results/PCSC1/Cluster/scABC/remap/Enhancer/FisherTest.txt", opdir="Enhancer", fc=4)