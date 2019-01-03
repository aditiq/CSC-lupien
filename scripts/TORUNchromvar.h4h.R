#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# h4h cluster
# Objective : Run chromVAR on scABC peaks
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

#---------------------------------
# load dependencies
#---------------------------------
library(data.table)

#---------------------------------
# Load function and countdata
#---------------------------------
source("scripts/chromvarfunc.R")
enhmat <-  fread("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt",check.names=F,stringsAsFactors = F,data.table=F, sep="\t", header=T)
enhmat$id <- paste(enhmat[,1], enhmat[,2], enhmat[,3], sep="_")
promat <-  fread("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoter.Binarymat.txt",check.names=F,stringsAsFactors = F,data.table=F, sep="\t", header=T)
promat$id <- paste(promat[,1], promat[,2], promat[,3], sep="_")

#---------------------------------
# Run chromVAR on scABC peaks
#---------------------------------

files <- dir("results/PCSC1/Cluster/scABC/", pattern="bed", full.names = T)
scabcdf <- list() 
for ( f in 1:length(files)){
  a <- read.table(files[f], header=F, sep="\t", stringsAsFactors = F)
  scabcdf[[f]] <- paste(a[,1],a[,2],a[,3],sep="_")
}

names(scabcdf) <- basename(files)

## Enhancer scABC cancer specific peaks
lsc.scabc.e <- scabcdf["LSC.Enhancer.p0.05.bed"][[1]]
gbm.scabc.e <- scabcdf["GBM.Enhancer.p0.05.bed"][[1]]
pfa.scabc.e <- scabcdf["PFA.Enhancer.p0.05.bed"][[1]]
shared.scabc.e <- scabcdf["Shared.Enhancer.p0.05.bed"][[1]]

lsc.scabc.p <- scabcdf["LSC.Promoter.p0.05.bed"][[1]]
gbm.scabc.p <- scabcdf["GBM.Promoter.p0.05.bed"][[1]]
pfa.scabc.p <- scabcdf["PFA.Promoter.p0.05.bed"][[1]]
shared.scabc.p <- scabcdf["Shared.Promoter.p0.05.bed"][[1]]

enhmat.scabc.cs <- subset(enhmat, enhmat$id %in% c(lsc.scabc.e, gbm.scabc.e, pfa.scabc.e))
enhmat.scabc.cs$id <- NULL
promat.scabc.cs <- subset(promat, promat$id %in% c(lsc.scabc.p, gbm.scabc.p, pfa.scabc.p))
promat.scabc.cs$id <- NULL

enhmat.scabc.shared <- subset(enhmat, enhmat$id %in% c(shared.scabc.e))
enhmat.scabc.shared$id <- NULL
promat.scabc.shared <- subset(promat, promat$id %in% c(shared.scabc.p))
promat.scabc.shared$id <- NULL

enhmat.scabc.lsc <- subset(enhmat, enhmat$id %in% c(lsc.scabc.e))
enhmat.scabc.lsc$id <- NULL
enhmat.scabc.pfa <- subset(enhmat, enhmat$id %in% c(pfa.scabc.e))
enhmat.scabc.pfa$id <- NULL
enhmat.scabc.gbm <- subset(enhmat, enhmat$id %in% c(gbm.scabc.e))
enhmat.scabc.gbm$id <- NULL

promat.scabc.lsc <- subset(promat, promat$id %in% c(lsc.scabc.p))
promat.scabc.lsc$id <- NULL
promat.scabc.pfa <- subset(promat, promat$id %in% c(pfa.scabc.p))
promat.scabc.pfa$id <- NULL
promat.scabc.gbm <- subset(promat, promat$id %in% c(gbm.scabc.p))
promat.scabc.gbm$id <- NULL

runchromvar(ds=enhmat.scabc.cs,dsflag=T,name="scABC.Enhancer.CancerSpecific", opdir="results/PCSC1/chromvar/scabc/enhancer/" )
runchromvar(ds=promat.scabc.cs,dsflag=T,name="scABC.Promoter.CancerSpecific", opdir="results/PCSC1/chromvar/scabc/promoter/" )

runchromvar(ds=enhmat.scabc.shared,dsflag=T,name="scABC.Enhancer.Shared", opdir="results/PCSC1/chromvar/scabc/enhancer/" )
runchromvar(ds=promat.scabc.shared,dsflag=T,name="scABC.Promoter.Shared", opdir="results/PCSC1/chromvar/scabc/promoter/" )

runchromvar(ds=enhmat.scabc.lsc,dsflag=T,name="scABC.Enhancer.LSC", opdir="results/PCSC1/chromvar/scabc/enhancer/" )
runchromvar(ds=promat.scabc.lsc,dsflag=T,name="scABC.Promoter.LSC", opdir="results/PCSC1/chromvar/scabc/promoter/" )

#runchromvar(ds=enhmat.scabc.gbm,dsflag=T,name="scABC.Enhancer.GBM", opdir="results/PCSC1/chromvar/scabc/enhancer/" )
#runchromvar(ds=promat.scabc.gbm,dsflag=T,name="scABC.Promoter.GBM", opdir="results/PCSC1/chromvar/scabc/promoter/" )

runchromvar(ds=enhmat.scabc.pfa,dsflag=T,name="scABC.Enhancer.PFA", opdir="results/PCSC1/chromvar/scabc/enhancer/" )
#runchromvar(ds=promat.scabc.pfa,dsflag=T,name="scABC.Promoter.PFA", opdir="results/PCSC1/chromvar/scabc/promoter/" )
