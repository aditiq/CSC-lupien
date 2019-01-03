#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Identify cluster specific peaks using scABC's function...
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
source("scripts/scABC.cluster.R") #source("https://raw.githubusercontent.com/timydaley/scABC/master/R/cluster.R")
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)


#----------------------------------------------------------------
# Define function to run scABC on different pairwise comparisons
#----------------------------------------------------------------


runscabcfunc=function(binarymatds, opname,...){
  
  #--------------------------------
  ## Read in binary matrix
  #--------------------------------
  binarymat <- read.table(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
  rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
  ncounts <- binarymat[,4:ncol(binarymat)]
  rownames(ncounts) <- rownames(binarymat)
  
  #--------------------------------
  ## Annotations
  #--------------------------------
  
  bg <- c(rep(1,48))
  all.three.cl <- c(rep(1,18),rep(2,23), rep(3,7))
  lsc.vs.gbmpfa.cl <- c(rep(1,18),rep(2,23), rep(2,7))
  lscgbm.vs.pfa.cl <- c(rep(1,18),rep(1,23), rep(2,7))
  lscpfa.vs.gbm.cl <- c(rep(1,18),rep(2,23), rep(1,7))

  names(bg) <- colnames(ncounts)
  
  names(all.three.cl) <- colnames(ncounts)
  names(lsc.vs.gbmpfa.cl) <- colnames(ncounts)
  names(lscgbm.vs.pfa.cl) <- colnames(ncounts)
  names(lscpfa.vs.gbm.cl) <- colnames(ncounts)
  
  
  #--------------------------------
  ## Run scabc 
  #--------------------------------
  all.three.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =all.three.cl, background_medians = bg)
  save(all.three.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSC.vs.GBM.vs.PFA.peakselection.", opname,".Rdata"))
  
  lsc.vs.gbmpfa.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =lsc.vs.gbmpfa.cl, background_medians = bg)
  save(lsc.vs.gbmpfa.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSC.vs.GBMPFA.peakselection.", opname,".Rdata"))
  
  lscgbm.vs.pfa.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =lscgbm.vs.pfa.cl, background_medians = bg)
  save(lscgbm.vs.pfa.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSCGBM.vs.PFA.peakselection.", opname,".Rdata"))
  
  lscpfa.vs.gbm.scabc = getClusterSpecificPvalue(ForeGround=as.matrix(ncounts), cluster_assignments =lscpfa.vs.gbm.cl, background_medians = bg)
  save(lscpfa.vs.gbm.scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.LSCPFA.vs.GBM.peakselection.", opname,".Rdata"))
  
  #--------------------------------
  ## Extract P values
  #--------------------------------
  
  all.three.scabc.pval = as.data.frame(all.three.scabc$pvalue)
  rownames(all.three.scabc.pval) <- rownames(ncounts)
  colnames(all.three.scabc.pval) <- c("LSC","GBM","PFA")
  all.three.scabc.pval$id <- rownames(all.three.scabc.pval)
  all.three.scabc.pval <- all.three.scabc.pval[order(all.three.scabc.pval$id),]
  
  lsc.vs.gbmpfa.scabc.pval = as.data.frame(lsc.vs.gbmpfa.scabc$pvalue)
  rownames(lsc.vs.gbmpfa.scabc.pval) <- rownames(ncounts)
  colnames(lsc.vs.gbmpfa.scabc.pval) <- c("LSCr","GBM.PFA") ## LSCr would be LSC specific and LSC- shared with GBM or PFA. Likewise for GBMr and PFAr
  lsc.vs.gbmpfa.scabc.pval$id <- rownames(lsc.vs.gbmpfa.scabc.pval)
  lsc.vs.gbmpfa.scabc.pval <- lsc.vs.gbmpfa.scabc.pval[order(lsc.vs.gbmpfa.scabc.pval$id),]
  
  lscgbm.vs.pfa.scabc.pval = as.data.frame(lscgbm.vs.pfa.scabc$pvalue)
  rownames(lscgbm.vs.pfa.scabc.pval) <- rownames(ncounts)
  colnames(lscgbm.vs.pfa.scabc.pval) <- c("LSC.GBM","PFAr")
  lscgbm.vs.pfa.scabc.pval$id <- rownames(lscgbm.vs.pfa.scabc.pval)
  lscgbm.vs.pfa.scabc.pval <- lscgbm.vs.pfa.scabc.pval[order(lscgbm.vs.pfa.scabc.pval$id),]
  
  lscpfa.vs.gbm.scabc.pval = as.data.frame(lscpfa.vs.gbm.scabc$pvalue)
  rownames(lscpfa.vs.gbm.scabc.pval) <- rownames(ncounts)
  colnames(lscpfa.vs.gbm.scabc.pval) <- c("LSC.PFA","GBMr")
  lscpfa.vs.gbm.scabc.pval$id <- rownames(lscpfa.vs.gbm.scabc.pval)
  lscpfa.vs.gbm.scabc.pval <- lscpfa.vs.gbm.scabc.pval[order(lscpfa.vs.gbm.scabc.pval$id),]
  
  #--------------------------------
  ## Combine dataset
  #--------------------------------

  combined <- merge(all.three.scabc.pval,lsc.vs.gbmpfa.scabc.pval[,c("id","GBM.PFA")], by.x="id", by.y="id",all.x=T)
  combined <- merge(combined,lscgbm.vs.pfa.scabc.pval[,c("id","LSC.GBM")], by.x="id", by.y="id",all.x=T)
  combined <- merge(combined,lscpfa.vs.gbm.scabc.pval[,c("id","LSC.PFA")], by.x="id", by.y="id",all.x=T)

  rownames(combined) <- combined$id
  write.table(combined, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".txt"), row.names=F, col.names=T, sep="\t", quote=F)

  #--------------------------------
  ## Plot images
  #--------------------------------

  dat <- c(subset(combined$id, combined$LSC<0.05),
           subset(combined$id, combined$GBM<0.05),
           subset(combined$id, combined$PFA<0.05),
           subset(combined$id, combined$GBM.PFA<0.05),
           subset(combined$id, combined$LSC.GBM<0.05),
           subset(combined$id, combined$LSC.PFA<0.05))
  
  remain1 <- combined[!(combined$id %in% dat),]
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[c(unique(dat), rownames(remain1)),]),2,rev)), col=scalered)
  dev.off()

  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.LSC.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$LSC<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.GBM.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$GBM<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.PFA.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$PFA<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.LSC.GBM.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$LSC.GBM<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.GBM.PFA.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$GBM.PFA<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.LSC.PFA.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[subset(combined$id, combined$LSC.PFA<0.05), ]),2,rev)), col=scalered)
  dev.off()
  
  pdf(paste0("results/PCSC1/Cluster/scABC/scABC.P0.05.Image.Remaining.sp.", opname,".pdf"))
  image(t(apply(as.matrix(ncounts[remain1$id, ]),2,rev)), col=scalered)
  dev.off()
  
}

#----------------------------------------------------------------
# Run function
#----------------------------------------------------------------

runscabcfunc("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt","Enhancer")
runscabcfunc("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoter.Binarymat.txt","Promoter")
