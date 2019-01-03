### R-3.4.1
### mordor
### Objective : KCCA Cluster analysis of binary matrix of different matrices following Wouter Muelman's script *Roadmap epigenomics DNAse clustering*

#-----------------------------------------
### load dependencies
#-----------------------------------------
library(flexclust)
library(methods)
library(Matrix)

source("scripts/wouter.kcca.R")


args <- commandArgs(trailingOnly = TRUE)
binarydsfile <- as.character(args[1])
opname <- as.character(args[2])
numk <- as.numeric(args[3])

#-----------------------------------------
## load data
#-----------------------------------------


### Define flexclust control
### flexClust control object, with custom distance function, which splits up computations.
ctrl <- list(verbose=1, iter=10000000)
ctrl <- as(ctrl ,"flexclustControl");
fam <- kccaFamily("ejaccard");

fam@dist <- function (x, centers) {
  if (ncol(x) != ncol(centers))
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")

  iter_size <- 5e5;
  num_iter <- ceiling(nrow(x)/iter_size);
  idxs <- rep(1:num_iter, length=nrow(x));
  z <- matrix(NA, nrow=nrow(x), ncol=nrow(centers));
  for (i in 1:num_iter) {
    sel <- which(idxs == i);
    x_sel <- x[sel,];

    xc <- x_sel %*% t(centers)
    nenner <- Matrix(rowSums(x_sel), nrow = nrow(x_sel), ncol = nrow(centers)) +
      Matrix(rowSums(centers), nrow = nrow(x_sel), ncol = nrow(centers), byrow = TRUE) - xc
    z_sel <- 1 - xc/nenner
    z_sel[nenner < sqrt(.Machine$double.eps)] <- 0
    z[sel,] <- as.matrix(z_sel);
  }
  z
}

#-----------------------------------------
### Perform clustering
#-----------------------------------------
 binarymat <- read.table(binarydsfile, header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")
mat <- as.matrix(binarymat[,4:ncol(binarymat)])

set.seed(1234)
kcca.cl <- kcca_Wouter(x=mat, k=numk, family=fam, control=ctrl,stop_percentage=0.0001);
save(kcca.cl, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.", opname,".",numk,".Rdata"))
rownames.obj <- rownames(mat)
save(rownames.obj, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.", opname,".rownames.",numk,".Rdata"))
