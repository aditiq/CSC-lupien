### R-3.4.1
### mordor
### Objective : Calculate jaccard distance for all regions with LSC, GBM and PFA to convert the data to euclidean distance
### NOT RUN

######################################
### load dependencies
######################################
library(prabclus)

args <- commandArgs(TRUE)
ds <- as.character(args[1])
pattern <- as.character(args[2])


######################################
## load data
######################################

matfull <- read.table(ds, sep="\t", header=F, stringsAsFactors = F, check.names=F)
rownames(matfull) <- paste0(matfull[,1],"_", matfull[,2],"_", matfull[,3])
mat <- matfull[,4:ncol(matfull)]

lsc <- c(rep(1,18))
gbm <- c(rep(1,23))
pfa <- c(rep(1,7))

jaccard.lsc.df <- NULL
jaccard.gbm.df <- NULL ; 
jaccard.pfa.df <- NULL ; 

for ( f in 1:nrow(mat))  { 
  print(f); 
  
  vec1=(as.numeric(mat[f,1:18]))
  jaccardlsc <- prabclus::jaccard(t(rbind(lsc, vec1)))[1,2]
  jaccard.lsc.df <- c(jaccard.lsc.df, jaccardlsc)
                      
  vec2=(as.numeric(mat[f,19:41]))
  jaccardgbm <- prabclus::jaccard(t(rbind(gbm, vec2)))[1,2]
  jaccard.gbm.df <- c(jaccard.gbm.df, jaccardgbm)
  
  vec3=(as.numeric(mat[f,42:48]))
  jaccardpfa <- prabclus::jaccard(t(rbind(pfa, vec3)))[1,2]
  jaccard.pfa.df <- c(jaccard.pfa.df, jaccardpfa)
}

jaccard.lsc.df2 <- as.data.frame(jaccard.lsc.df)
jaccard.gbm.df2 <- as.data.frame(jaccard.gbm.df)
jaccard.pfa.df2 <- as.data.frame(jaccard.pfa.df)

jaccard.lsc.df2$id <- rownames(mat)
jaccard.gbm.df2$id <- rownames(mat)
jaccard.pfa.df2$id <- rownames(mat)

  
write.table(jaccard.lsc.df2, file=paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/jaccard.LSC.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(jaccard.gbm.df2, file=paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/jaccard.GBM.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(jaccard.pfa.df2, file=paste0("results/PCSC1/Cluster/Jaccard/jaccard.split/jaccard.PFA.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)

