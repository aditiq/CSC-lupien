### R-3.4.1
### mordor
### Objective : Clustering using -
###  TF-IDF followed by LSI for dimension reduction
### Maybe graph based clustering
### NOT RUN

#####################################
#### USEFUL LINKS
#####################################
## http://bwlewis.github.io/1000_genomes_examples/PCA.html
## https://bwlewis.github.io/irlba/comparison.html
## http://atlas.gs.washington.edu/fly-atac/docs/#use-case-1-id-clades-with-lsi
## 
######################################
### load dependencies
######################################
library(Matrix)
library(proxy)
library(gplots)
library(Rtsne)
library(densityClust)
library(irlba)
library(pheatmap)
color_pal = c("#1F78B4","#FFD700","#60CC52","#E31A1C")
hmcols = colorpanel(100, "steelblue", "white", "tomato")


######################################
## Run analysis
######################################

## Load and prepare data
binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)

num_cells_ncounted = rowSums(binarymat[,4:ncol(binarymat)])
options(repr.plot.width=4, repr.plot.height=4)
pdf("results/PCSC1/Cluster/LSI/AccessibilityPlot.pdf")
hist((num_cells_ncounted),main="No. of Samples Each Site is Observed In",breaks=50, xlab="No. of samples")
dev.off()

rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
ncounts <- binarymat[,4:ncol(binarymat)]
rownames(ncounts) <- rownames(binarymat)

## Normalisation and scaling of binary matrix -- 
nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
idf = as(log(1 + ncol(ncounts) / Matrix::rowSums(ncounts)), "sparseVector")
tf_idf_counts = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs

## Dimensionality reduction using partial SVD
## We have to identify the V matrix (http://www1.se.cuhk.edu.hk/~seem5680/lecture/LSI-Eg.pdf)

## How many dimensions should be considered ?
## We first need to assess number of dimensions to keep. This can be estimated from the S matrix 
## (https://stackoverflow.com/questions/9582291/how-do-we-decide-the-number-of-dimensions-for-latent-semantic-analysis)

## Running first with 20 dimensions to check how much of the variance is being captured
## (http://genomicsclass.github.io/book/pages/svd.html)
set.seed(0) #For reproducibility
SVD = irlba(tf_idf_counts, 47,47,verbose=T) 

## Plot Sum of Squares of UD

pdf("results/PCSC1/Cluster/LSI/SVD.Variability.pdf")
plot(cumsum(SVD$d^2)/sum(SVD$d^2)*100,ylab="Percent variability explained",ylim=c(0,100),type="l")
plot(SVD$d^2/sum(SVD$d^2)*100,ylab="Percent variability explained")
dev.off()

## 40 dimensions -- 95% variability
## 30 dimensions -- 86% variability

set.seed(0) #For reproducibility
k=30
SVD30 = irlba(tf_idf_counts, k,k,verbose=T) 
sk_diag = matrix(0, nrow=k, ncol=k)
diag(sk_diag) = SVD30$d
sk_diag[1,1] = 0 ## Since component 1 is associated with read depth

LSI_out = t(t(sk_diag %*% t(SVD30$v)) %*% t(SVD30$u))
colnames(LSI_out) <- colnames(ncounts)
rownames(LSI_out) <- rownames(binarymat)

LSI_out = t(scale(t(LSI_out)))
LSI_out[LSI_out > 3] = 3
LSI_out[LSI_out < -3] = -3

write.table(LSI_out, file="results/PCSC1/Cluster/LSI/LSI.SVD30.matrix.txt", sep="\t", row.names=T, col.names=T,quote=F )
saveRDS(LSI_out, file="LSI_out.SVD30.rds")

### Now that we have a reduced dimensionality matrix in the euclidean space, we can think about clustering it

## Hierarchial clustering would have been ideal but 
## https://stackoverflow.com/questions/42479854/merge-error-negative-length-vectors-are-not-allowed
## hclust_cells = hclust(proxy::dist(t(sk_diag %*% t(SVD30$v)), method="cosine"), method="ward.D2")
## hclust_genes = hclust(proxy::dist(t(sk_diag %*% t(SVD30$u)), method="cosine"), method="ward.D2")

## Ways to cluster it --
## Construct a KNN graph in cosine and use Louvain clustering .. ie. graph based clustering

## Find neighbours using cosine similarity -- Can use FALCONN
## Compute Jaccard coeff between neighbours
## Build graph



#### Alternative --- Finite mixture models 
pdf("results/PCSC1/Cluster/LSI/TF_IDF_dist.pdf")
