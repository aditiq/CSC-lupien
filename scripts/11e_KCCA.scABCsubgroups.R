### R-3.4.1
### mordor
### Objective : KCCA Cluster analysis of regions grouped by scABC following Wouter Muelman's script *Roadmap epigenomics DNAse clustering*

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

run_scabc_kcca <- function(matds, opname,...){
  
  binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)
  rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")
  ds <- read.table(matds, header=F, sep="\t", stringsAsFactors = F,check.names = F)
  rownames(ds) <- paste(ds[,1],ds[,2],ds[,3],sep="_")
  
  mat <- as.matrix(binarymat[rownames(ds),4:ncol(binarymat)])
  
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
  
  set.seed(1234)
  kcca.cl <- kcca_Wouter(x=mat, k=100, family=fam, control=ctrl,stop_percentage=0.01);
  save(kcca.cl, file=paste0(opname, ".100.Rdata"))
  
}

### Perform clustering

files <- dir("results/PCSC1/Cluster/scABC/", ".p0.05.bed", full.names=T)
files2 <- files[-grep("hg19", files)]

for (f  in 1:length(files2)){
  run_scabc_kcca( matds=files2[f],
                 opname=gsub(".bed","", basename(files2[f])) )
}
