#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Run chromVAR
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------
## load dependencies
#----------------------------------------------------------------
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2016)
#register(MulticoreParam(8, progressbar = TRUE))
register(SerialParam())
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gplots)

## colors from devtools::install_github("caleblareau/BuenColors")
#library(BuenColors)
brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')
solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D')
Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

differentialDeviations2 <- function(object,
                                    groups,
                                    alternative = c("two.sided", "less",
                                                    "greater"),
                                    parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))) {
    groups <- colData(object)[[groups]]
  } else if (length(groups) != ncol(object)) {
    stop("invalid groups input, must be vector of lench ncol(object) or column",
         " name from colData(object)")
  }

  groups <- as.factor(groups)

  alternative <- match.arg(alternative)
  inputs <- deviations(object)
  inputs <- na.omit(inputs)
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }

  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}


t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]],splitx[[2]],
                alternative = alternative,
                paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- oneway.test(devs ~ groups, tmpdf, var.equal = FALSE)
  return(res$p.value)
}

kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
                     alternative = alternative,
                     paired = FALSE)$p.value)
}

differentialVariability2 <- function (object, groups, parametric = TRUE)
{
    stopifnot(is(object, "chromVARDeviations"))
    if (length(groups) == 1 && groups %in% colnames(colData(object))) {
        groups <- colData(object)[[groups]]
    }
    else if (length(groups) != ncol(object)) {
        stop("invalid groups input, must be vector of lench ncol(object) or column",
            " name from colData(object)")
    }
    groups <- as.factor(groups)
    inputs <- deviationScores(object)
    inputs <- na.omit(inputs)

    if (parametric) {
        p_val <- apply(inputs, 1, bf_var_test, groups)
    }
    else {
        p_val <- apply(inputs, 1, bf_kw_var_test, groups)
    }
    p_adj <- p.adjust(p_val, method = "BH")
    return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}



bf_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(anova(lm(median_diff ~ groups))[1, 5])
}

bf_kw_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(kruskal.test(median_diff ~ groups)$p.value)
}

#----------------------------------------------------------------
## Load Countdata
#----------------------------------------------------------------

#----------------------------------------------------------------
## Load Countdata
#----------------------------------------------------------------
set.seed(2017)
my_counts_matrix <- fread("data/ConsensusSet/PCSC1/Combined.Blood.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
## Keep only for LSC positive and differentiated hemat.

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping$stem <- ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2)
mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2))
mapping <- mapping[gsub("_peaks","",colnames(my_counts_matrix)[4:ncol(my_counts_matrix)]),]

mapping <- subset(mapping, mapping$tissue %in% c("pos.LSC","diff.Hemat"))
colnames(my_counts_matrix) <- gsub("_peaks","",colnames(my_counts_matrix))
my_counts_matrix <- my_counts_matrix[,c("seqnames","start","end",rownames(mapping))]
my_counts_matrix <- subset(my_counts_matrix, rowSums(my_counts_matrix[,4:ncol(my_counts_matrix)]) >0 )
rownames(my_counts_matrix) <- paste(my_counts_matrix[,1],my_counts_matrix[,2],my_counts_matrix[,3], sep="_")
fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix)
assayNames(fragment_counts) <- "counts"
colData(fragment_counts)$Stem <- as.factor(mapping$stem)
colData(fragment_counts)$Group <- as.factor(mapping$group1)

#----------------------------------------------------------------
## add gc content
#----------------------------------------------------------------

fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
#counts_filtered <- filterPeaks(fragment_counts,min_fragments_per_peak = 10, non_overlapping = TRUE)
counts_filtered <- fragment_counts
rm(fragment_counts)
save(counts_filtered, file="results/PCSC1/chromVAR/remap/blood/counts_filtered.Rdata")

#----------------------------------------------------------------
## Get motifs and what peaks contain motifs
#----------------------------------------------------------------
my_annotation_df <- fread("data/ConsensusSet/PCSC1/Combined.Blood.Consensus.Catalogue.Binarymat.with.remap.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
rownames(my_annotation_df) <- paste(my_annotation_df[,1],my_annotation_df[,2],my_annotation_df[,3], sep="_")
my_annotation_df <- my_annotation_df[rownames(my_counts_matrix),]

anno_ix <- getAnnotations(as.matrix(my_annotation_df[,4:ncol(my_annotation_df)]), rowRanges = rowRanges(counts_filtered))
save(anno_ix, file="results/PCSC1/chromVAR/remap/blood/anno_ix.Rdata")

#kmer_ix <- matchKmers(6, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)
#rm(motifs)
#save("motif_ix",file="results/All/chromVAR/repeats/All.extended.motif.ixe.Rdata")
#dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
#bg <- getBackgroundPeaks(object = counts_filtered)
#expected <- computeExpectations(counts_filtered, norm=F)

#----------------------------------------------------------------
## compute deviation
#----------------------------------------------------------------
dev <- computeDeviations(object = counts_filtered, annotations = anno_ix)
save(dev,file="results/PCSC1/chromVAR/remap/blood/Blood.extended.dev.Rdata")
z.scores = deviationScores(dev) ## deviation Z-score
dev.scores = deviations(dev) ## bias corrected deviations

write.table(z.scores, file="results/PCSC1/chromVAR/remap/blood/Dev.Zscore.txt", col.names=T, row.names=T, sep="\t", quote=F)
write.table(dev.scores, file="results/PCSC1/chromVAR/remap/blood/Deviations.txt", col.names=T, row.names=T, sep="\t", quote=F)

#----------------------------------------------------------------
## compute variablity
#----------------------------------------------------------------
variability <- computeVariability(dev)
pdf("results/PCSC1/chromVAR/remap/blood/Variability.pdf") ;
plotVariability(variability, use_plotly = FALSE) ;
dev.off()

write.table(variability, file="results/PCSC1/chromVAR/remap/blood/Variability.txt", col.names=T, row.names=F, sep="\t", quote=F)

#----------------------------------------------------------------
## Plotting with inlfection point calculated using the SUPER ENHANCER ROSE ALGO
#----------------------------------------------------------------
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

inputVector <- sort(variability$variability)
slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.


signalOrder = order(inputVector,decreasing=TRUE)
res_df <- inputVector[signalOrder]
topvartfs <- which(res_df> y_cutoff)

pdf("results/PCSC1/chromVAR/remap/blood/VariabilityRanked.pdf", useDingbats = F) ;
plot(length(res_df):1,res_df, col='dodgerblue2',xlab="Ranked TFs",ylab="Variability",pch=16)
abline(h=y_cutoff,col='grey',lty=2)
abline(v=length(res_df)-length(topvartfs),col='grey',lty=2)
#lines(length(res_df):1,res_df,lwd=4, col='dodgerblue2')
dev.off()

#----------------------------------------------------------------
## Heatmap of high variability TFs
#----------------------------------------------------------------

variability.df <- as.data.frame(variability)
variability.df <-  variability.df[order(-variability.df$variability),]

top_motifs = subset(variability.df$name, variability.df$p_value_adj < 0.01 & variability.df$variability > y_cutoff)
top_devs = z.scores[which(rownames(dev.scores) %in% (top_motifs)), ]

annocol1 <- data.frame(group=c(rep("LSCp",18), rep("Hemat.diff",34)), stringsAsFactors=F)
rownames(annocol1) <- colnames(top_devs)
rownames(top_devs) <- gsub(".nr_macs2_hg38_v1_2","", rownames(top_devs) )
colnames(top_devs) <- (data.frame(str_split_fixed((data.frame(str_split_fixed(colnames(top_devs),"-",2), stringsAsFactors=F))$X1,"_",2), stringsAsFactors=F))$X1

pdf("results/PCSC1/chromVAR/remap/blood/Heatmap.BiasCorrectedDev.Top.pdf")
pheatmap((top_devs),
          clustering_method="complete",
          clustering_distance_rows="euclidean",
          clustering_distance_cols="euclidean",
          col = brewer_yes,
          fontsize_col=5,fontsize_row=5,
          annotation_col=annocol1,
          scale="none")
dev.off()

#----------------------------------------------------------------
## T-SNE  plots
#----------------------------------------------------------------
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 30)
tsne_plots <- plotDeviationsTsne(dev, tsne_results,   annotation_name = NULL,   sample_column = "Stem", shiny = FALSE)
pdf("results/PCSC1/chromVAR/remap/blood/Tsne.stem.pdf"); print(tsne_plots) ; dev.off()

tsne_plots2 <- plotDeviationsTsne(dev, tsne_results,   annotation_name = NULL,   sample_column = "Group", shiny = FALSE)
pdf("results/PCSC1/chromVAR/remap/blood/Tsne.group.pdf"); print(tsne_plots2) ; dev.off()

for( f in 1:nrow(variability) ){
  tsne_plots <- plotDeviationsTsne(dev, tsne_results,   annotation_name = as.character(variability[f,1]),   sample_column = c("Group","Stem"), shiny = FALSE)
  pdf(paste0("results/PCSC1/chromVAR/remap/blood/tsne.tfplots/Tsne.",as.character(variability[f,1]),".pdf")); print(tsne_plots) ; dev.off()
}

#----------------------------------------------------------------
## Sample correlation
#----------------------------------------------------------------

sample_cor <- getSampleCorrelation(dev)
rownames(sample_cor) <- gsub("_peaks","", rownames(sample_cor))
colnames(sample_cor) <- gsub("_peaks","", colnames(sample_cor))
rownames(colData(dev)) <- gsub("_peaks","", rownames(colData(dev)))

pdf("results/PCSC1/chromVAR/remap/blood/SampleCorrelation.pdf")
pheatmap((sample_cor),
         annotation_row = colData(dev)[,1:2],  annotation_col=colData(dev)[,1:2],cex=0.5,
         clustering_distance_rows = as.dist(1-sample_cor),
         clustering_distance_cols = as.dist(1-sample_cor), labels_row=mapping$group1)
dev.off()

#----------------------------------------------------------------
## Differential accessibilty and variations
#----------------------------------------------------------------

## Function to see whether deviations differ between groups-- CSC and differentiated

diff_acc <- differentialDeviations2(dev, 'Stem')
diff_acc <- diff_acc[order(diff_acc$p_value_adjusted),]
diff_acc$repname <- rownames(diff_acc)
write.table(diff_acc, file="results/PCSC1/chromVAR/remap/blood/DiffAccessibility.txt", col.names=T, row.names=F, sep="\t", quote=F)

## Difference in variability
diff_var <- differentialVariability2(dev, "Stem", parametric = FALSE)
diff_var <- diff_var[order(diff_var$p_value_adjusted),]
diff_var$repname <- rownames(diff_var)
write.table(diff_var, file="results/PCSC1/chromVAR/remap/blood/DiffVariability.txt", col.names=T, row.names=F, sep="\t", quote=F)

colnames(z.scores) <- gsub("_peaks","", colnames(z.scores))
colnames(dev.scores) <- gsub("_peaks","", colnames(dev.scores))
annocol <- as.data.frame(colData(dev)@listData)
rownames(annocol) <- colnames(dev.scores)
#annocol$name <- rownames(annocol)

pdf("results/PCSC1/chromVAR/remap/blood/Heatmap.BiasCorrectedDev.pdf")
pheatmap(na.omit(dev.scores),
          annotation_col=annocol,clustering_method="ward.D2",
          clustering_distance_rows="euclidean",
          clustering_distance_cols="euclidean",
         cex=0.5,color=brewer_yes,scale="none",
         cluster_rows=T, cluster_cols=T)
dev.off()

pdf("results/PCSC1/chromVAR/remap/blood/Heatmap.ZscoreDev.pdf")
pheatmap(na.omit(z.scores),
         annotation_col=annocol,clustering_method="ward.D2",
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         cex=0.5,color=brewer_yes,scale="none",
         cluster_rows=T, cluster_cols=T)
dev.off()
