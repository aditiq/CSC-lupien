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

#----------------------------------------------------------------
## Load Countdata
#----------------------------------------------------------------
set.seed(2017)
my_counts_matrix <- fread("results/PCSC1/CREAM/PCSC1.CREAM.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
rownames(my_counts_matrix) <- paste(my_counts_matrix[,1],my_counts_matrix[,2],my_counts_matrix[,3], sep="_")

## remove random chromosomes
toMatch <- c("random", "alt", "chrUn", "chrM")
my_counts_matrix2 <- subset(my_counts_matrix, !(grepl(paste(toMatch, collapse="|"), my_counts_matrix$seqnames)))

fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix2)
assayNames(fragment_counts) <- "counts"

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping$stem <- ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2)
mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2))
mapping <- mapping[gsub("_peaks","",colnames(fragment_counts)),]
colData(fragment_counts)$Stem <- as.factor(mapping$stem)
colData(fragment_counts)$Group <- as.factor(mapping$group1)

#----------------------------------------------------------------
## add gc content
#----------------------------------------------------------------

fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
#counts_filtered <- filterPeaks(fragment_counts,min_fragments_per_peak = 10, non_overlapping = TRUE)
counts_filtered <- fragment_counts

#----------------------------------------------------------------
## Get motifs and what peaks contain motifs
#----------------------------------------------------------------
#----------------------------------------------------------------
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered,   genome = BSgenome.Hsapiens.UCSC.hg38)
save(motif_ix,file="results/PCSC1/chromVAR/CREAM/PCSC1/motif_ix.Rdata")

#----------------------------------------------------------------
## compute deviation
#----------------------------------------------------------------
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
save(dev,file="results/PCSC1/chromVAR/CREAM/PCSC1/Combined.PCSC1.extended.dev.Rdata")
z.scores = deviationScores(dev) ## deviation Z-score
dev.scores = deviations(dev) ## bias corrected deviations

write.table(z.scores, file="results/PCSC1/chromVAR/CREAM/PCSC1/Dev.Zscore.txt", col.names=T, row.names=T, sep="\t", quote=F)
write.table(dev.scores, file="results/PCSC1/chromVAR/CREAM/PCSC1/Deviations.txt", col.names=T, row.names=T, sep="\t", quote=F)

rownames(dev.scores) <- substr(rownames(dev.scores), 10, length(rownames(dev.scores)))
rownames(z.scores) <- substr(rownames(z.scores), 10, length(rownames(z.scores)))

#----------------------------------------------------------------
## compute variablity
#----------------------------------------------------------------
variability <- computeVariability(dev)
pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Variability.pdf") ;
plotVariability(variability, use_plotly = FALSE) ;
dev.off()

write.table(variability, file="results/PCSC1/chromVAR/CREAM/PCSC1/Variability.txt", col.names=T, row.names=F, sep="\t", quote=F)

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

pdf("results/PCSC1/chromVAR/CREAM/PCSC1/VariabilityRanked.pdf", useDingbats = F) ;
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

rownames(top_devs) <- gsub(".nr_macs2_hg38_v1_2","", rownames(top_devs) )
colnames(top_devs) <- gsub("-includes-blacklisted-regions_peaks.hg38","", colnames(top_devs))

annocol1 <- data.frame(group=c(rep("LSC",18),
                                rep("GBM",23),
                                rep("PFA",7)), stringsAsFactors=F)

rownames(annocol1) <- colnames(top_devs)

pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Heatmap.BiasCorrectedDev.Top.pdf")
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
pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Tsne.stem.pdf"); print(tsne_plots) ; dev.off()

tsne_plots2 <- plotDeviationsTsne(dev, tsne_results,   annotation_name = NULL,   sample_column = "Group", shiny = FALSE)
pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Tsne.group.pdf"); print(tsne_plots2) ; dev.off()

for( f in 1:nrow(variability) ){
  tsne_plots <- plotDeviationsTsne(dev, tsne_results,   annotation_name = as.character(variability[f,1]),   sample_column = c("Group","Stem"), shiny = FALSE)
  pdf(paste0("results/PCSC1/chromVAR/CREAM/PCSC1/tsne.tfplots/Tsne.",as.character(variability[f,1]),".pdf")); print(tsne_plots) ; dev.off()
}

#----------------------------------------------------------------
## Sample correlation
#----------------------------------------------------------------

sample_cor <- getSampleCorrelation(dev)
rownames(sample_cor) <- gsub("_peaks","", rownames(sample_cor))
colnames(sample_cor) <- gsub("_peaks","", colnames(sample_cor))
rownames(colData(dev)) <- gsub("_peaks","", rownames(colData(dev)))

pdf("results/PCSC1/chromVAR/CREAM/PCSC1/SampleCorrelation.pdf")
pheatmap((sample_cor),
         annotation_row = colData(dev)[,1:2],  annotation_col=colData(dev)[,1:2],cex=0.5,
         clustering_distance_rows = as.dist(1-sample_cor),
         clustering_distance_cols = as.dist(1-sample_cor), labels_row=mapping$group1)
dev.off()

#----------------------------------------------------------------
## Differential accessibilty and variations
#----------------------------------------------------------------

## Function to see whether deviations differ between groups
diff_acc <- differentialDeviations(dev, "Stem")
diff_acc <- diff_acc[order(diff_acc$p_value_adjusted),]
rownames(diff_acc) <- substr(rownames(diff_acc),10, nchar(rownames(diff_acc)))

## Difference in variability
diff_var <- differentialVariability(dev, "Stem", parametric = FALSE)

colnames(z.scores) <- gsub("_peaks","", colnames(z.scores))
colnames(dev.scores) <- gsub("_peaks","", colnames(dev.scores))
annocol <- as.data.frame(colData(dev)@listData)
rownames(annocol) <- colnames(dev.scores)
#annocol$name <- rownames(annocol)

pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Heatmap.BiasCorrectedDev.pdf")
pheatmap(dev.scores,
         annotation_col=annocol, clustering_method="complete",
         clustering_distance_rows="correlation",clustering_distance_cols="euclidean",
         cex=0.5,color=brewer_yes,scale="none",
         cluster_rows=T, cluster_cols=T)
dev.off()

pdf("results/PCSC1/chromVAR/CREAM/PCSC1/Heatmap.ZscoreDev.pdf")
pheatmap(z.scores,
         annotation_col=annocol,clustering_method="complete",
         clustering_distance_rows="correlation",clustering_distance_cols="euclidean",
         cex=0.5,color=brewer_yes,scale="none",
         cluster_rows=T, cluster_cols=T)
dev.off()
