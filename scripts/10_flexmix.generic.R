### R-3.4.1
### mordor
### Objective : Clustering using finite mixture model for binary data
### NOT RUN -- Needs to be optimised

args <- commandArgs(trailingOnly = TRUE)
ds <- as.character(args[1])
opname <- as.character(args[2])
num1 <- as.numeric(args[3])
num2 <- as.numeric(args[4])

######################################
### load dependencies
######################################
library(flexmix)          # For model based clustering

######################################
### Run analysis
######################################

## load data
binarymat <- read.table(ds, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
mat <- (binarymat[,4:ncol(binarymat)])

## Creates the input for EM
test.df <- data.frame(Freq=as.integer(rep(1,dim(as.matrix(mat))[1])))
test.df$Incidence <- as.matrix(mat)

## Perform model based clustering
## The probability of the data fitting the model for each of the k parameter is stored 
## and will be used for inferring the number of clones.

set.seed(0)
enhmat_fmm <- flexmix::stepFlexmix(Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = test.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list(minprior  = 0.005,iter=100000), 
                                   k       = num1:num2, 
                                   nrep    = 3)  
save(enhmat_fmm, file=paste0("results/PCSC1/Cluster/Flexmix/Binomial/", opname,".Binarymat.Fmm.K",num1,".",num2,".Rdata"))
