### R-3.4.1
### mordor
### Objective : MCC of binary matrix
### NOT RUN

######################################
### load dependencies
######################################
library(PharmacoGx)
args <- commandArgs(TRUE)
ds <- as.character(args[1])
pattern <- as.character(args[2])


######################################
## load data
######################################

matfull <- read.table(ds, sep="\t", header=F, stringsAsFactors = F, check.names=F)
rownames(matfull) <- paste0(matfull[,1],"_", matfull[,2],"_", matfull[,3])
mat <- subset(matfull[,4:ncol(matfull)], rowSums(matfull[,4:ncol(matfull)])<48)

lsc <- as.factor(c(rep(1,18),rep(0,23),rep(0,7)))
gbm <- as.factor(c(rep(0,18),rep(1,23),rep(0,7)))
pfa <- as.factor(c(rep(0,18),rep(0,23),rep(1,7)))

mcc.lsc.df <- NULL
mcc.gbm.df <- NULL ; 
mcc.pfa.df <- NULL ; 

for ( f in 1:nrow(mat))  { 
  print(f); 
  vec1=as.factor(as.numeric(mat[f,1:48]))
  
  mcclsc <- mcc(vec1, lsc, nperm=1,nthread=1)$estimate ; 
    mcc.lsc.df <- c(mcc.lsc.df, mcclsc)
                      
  mccgbm <- mcc(vec1, gbm, nperm=1, nthread=1) ; 
  mcc.gbm.df <- c(mcc.gbm.df, mccgbm$estimate)
                      
  mccpfa <- mcc(vec1, pfa,  nperm=1,nthread=1) ; 
  mcc.pfa.df <- c(mcc.pfa.df, mccpfa$estimate)
}

mcc.lsc.df2 <- as.data.frame(mcc.lsc.df)
mcc.gbm.df2 <- as.data.frame(mcc.gbm.df)
mcc.pfa.df2 <- as.data.frame(mcc.pfa.df)

mcc.lsc.df2$id <- rownames(mat)
mcc.gbm.df2$id <- rownames(mat)
mcc.pfa.df2$id <- rownames(mat)

  
write.table(mcc.lsc.df2, file=paste0("results/PCSC1/MCC/mcc.split/MCC.LSC.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(mcc.gbm.df2, file=paste0("results/PCSC1/MCC/mcc.split/MCC.GBM.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
write.table(mcc.pfa.df2, file=paste0("results/PCSC1/MCC/mcc.split/MCC.PFA.",pattern,".txt"), sep="\t", row.names=F, col.names=T, quote=F)

