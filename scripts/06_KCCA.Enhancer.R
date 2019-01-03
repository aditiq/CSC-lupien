### R-3.4.1
### mordor
### Objective : KCCA Cluster analysis of binary matrix following Wouter Muelman's script *Roadmap epigenomics DNAse clustering*

args <- commandArgs(trailingOnly = TRUE)
numk <- as.numeric(args[1])


######################################
### load dependencies
######################################
library(flexclust)
library(methods)
library(Matrix)

source("scripts/wouter.kcca.R")

##############################################################################################
## load data
##############################################################################################

binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")
enh <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancers.bed", header=F, sep="\t", stringsAsFactors = F,check.names = F)
rownames(enh) <- paste(enh[,1],enh[,2],enh[,3],sep="_")

mat <- as.matrix(binarymat[rownames(enh),4:ncol(binarymat)])

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

### Perform clustering
set.seed(1234)
kcca.cl <- kcca_Wouter(x=mat, k=numk, family=fam, control=ctrl,stop_percentage=0.0001);
save(kcca.cl, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Enhancer.", numk,".Rdata"))
rownames.obj <- rownames(mat)
save(rownames.obj, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.Enhancer.rownames.", numk,".Rdata"))

## Perform clustering after removing singletons
mat2 <- as.data.frame(mat)
mat2$LSC <- rowSums(mat2[,1:18])
mat2$GBM <- rowSums(mat2[,19:41])
mat2$PFA <- rowSums(mat2[,42:48])

mat3 <- subset(mat2[,1:48], mat2$LSC>1 | mat2$GBM>1 | mat2$PFA >1)
set.seed(1234)
numk=100
kcca.cl <- kcca_Wouter(x=as.matrix(mat3), k=numk, family=fam, control=ctrl,stop_percentage=0.0001);
save(kcca.cl, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.EnhancerwdtSingletons.", numk,".Rdata"))
rownames.obj2 <- rownames(mat3)
save(rownames.obj2, file=paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.EnhancerwdtSingletons.rownames.", numk,".Rdata"))
