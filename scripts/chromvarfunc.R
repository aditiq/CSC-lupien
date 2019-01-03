#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# h4h cluster
# Objective : generic function for chromvar 
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

#---------------------------------
# load dependencies
#---------------------------------
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2016)
register(SerialParam())
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gplots)
#register(MulticoreParam(8, progressbar = TRUE))

## colors from devtools::install_github("caleblareau/BuenColors")
#library(BuenColors)
brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')
solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D')
Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

#---------------------------------
# Define function
#---------------------------------

runchromvar=function(ds,dsflag, name,opdir,...){
  set.seed(2017)
  
  if(dsflag==T){
    my_counts_matrix <- ds
  } else {
    my_counts_matrix <- fread(ds,check.names=F,stringsAsFactors = F,data.table=F, sep="\t", header=T)
  }
  rownames(my_counts_matrix) <- paste(my_counts_matrix[,1],my_counts_matrix[,2],my_counts_matrix[,3], sep="_")
  fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix)
  assayNames(fragment_counts) <- "counts"
  
  
  colData(fragment_counts)$Cell_Type <- as.factor(c(rep("LSC",18), rep("GBM",23), rep("PFA",7)))
  ## add gc content
  fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  counts_filtered <- fragment_counts
  
  ## Get motifs and what peaks contain motifs
  motifs <- getJasparMotifs()
  motif_ix <- matchMotifs(motifs, counts_filtered,   genome = BSgenome.Hsapiens.UCSC.hg38)
  dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
  
  save("dev",file=paste0(opdir, name, ".extended.dev.Rdata"))
  
  ## compute deviation
  z.scores = deviationScores(dev) ## deviation Z-score
  dev.scores = deviations(dev) ## bias corrected deviations
  
  write.table(z.scores, file=paste0(opdir, name, ".dev.Zscore.txt"), col.names=T, row.names=T, sep="\t", quote=F)
  write.table(dev.scores, file=paste0(opdir, name, ".Deviations.txt"), col.names=T, row.names=T, sep="\t", quote=F)
  
  rownames(dev.scores) <- substr(rownames(dev.scores), 10, length(rownames(dev.scores)))
  rownames(z.scores) <- substr(rownames(z.scores), 10, length(rownames(z.scores)))
  
  
  ## compute variablity
  variability <- computeVariability(dev)
  pdf(paste0(opdir, name,".variability.pdf")) ;
  plotVariability(variability, use_plotly = FALSE) ; 
  dev.off()
  
  write.table(variability, file=paste0(opdir, name, ".variability.txt"), col.names=T, row.names=F, sep="\t", quote=F)
  
  ## Plotting with inlfection point calculated using the SUPER ENHANCER ROSE ALGO
  
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
  
  pdf(paste0(opdir, name, ".variability2.pdf"), useDingbats = F) ;
  plot(length(res_df):1,res_df, col='dodgerblue2',xlab="Ranked TFs",ylab="Variability",pch=16)	
  abline(h=y_cutoff,col='grey',lty=2)
  abline(v=length(res_df)-length(topvartfs),col='grey',lty=2)
  #lines(length(res_df):1,res_df,lwd=4, col='dodgerblue2')
  dev.off()
  
  
  ## Heatmap of high variability TFs
  variability.df <- as.data.frame(variability)
  variability.df <-  variability.df[order(-variability.df$variability),]
  
  top_motifs = subset(variability.df$name, variability.df$variability>y_cutoff)
  top_devs = dev.scores[which(rownames(dev.scores) %in% (top_motifs)), ]
  
  pdf(paste0(opdir, name, ".Heatmap.BiasCorrectedDev.Top.pdf"))
  heatmap.2((top_devs),
            Rowv=T, Colv=T,
            trace='none', col = brewer_yes,
            density.info = "none", 
            scale="none")
  dev.off()
  
  
  ## tsne
  tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)
  
  ## T-SNE
  tsne_plots <- plotDeviationsTsne(dev, tsne_results,   annotation_name = NULL,   sample_column = "Cell_Type", shiny = FALSE)
  pdf(paste0(opdir, name,".Tsne.pdf")); print(tsne_plots) ; dev.off()
  
  # for( f in 1:nrow(variability) ){
  #   tsne_plots <- plotDeviationsTsne(dev, tsne_results,   annotation_name = as.character(variability[f,1]),   sample_column = "Cell_Type2", shiny = FALSE)
  #   pdf(paste0("results/.Tsne.",as.character(variability[f,1]),".pdf")); print(tsne_plots) ; dev.off()
  # }
  
  ## Sample correlation
  sample_cor <- getSampleCorrelation(dev)
  
  pdf(paste0(opdir, name,".SampleCorrelation.pdf")); 
  pheatmap((sample_cor), 
           annotation_row = colData(dev),  annotation_col=colData(dev),cex=0.5,
           clustering_distance_rows = as.dist(1-sample_cor), 
           clustering_distance_cols = as.dist(1-sample_cor))
  dev.off()
  
  
  ## Differential accessibilty and variations
  
  ## Difference in bias corrected deviations for motifs
  #diff_acc <- differentialDeviations(dev, "Cell_Type")
  #diff_acc <- diff_acc[order(diff_acc$p_value_adjusted),]
  #rownames(diff_acc) <- substr(rownames(diff_acc),10, nchar(rownames(diff_acc)))
  
  ## Difference in variability of chromatin accessibility associated with motifs 
  #diff_var <- differentialVariability(dev, "Cell_Type", parametric = FALSE)  
  
  pdf(paste0(opdir, name,".Heatmap.BiasCorrectedDev.pdf")); 
  pheatmap(dev.scores,
           #annotation_col=colData(dev)[,1:2],
           cex=0.5,color=brewer_yes,scale="none",
           cluster_rows=T, cluster_cols=T)
  dev.off()
  
  pdf(paste0(opdir, name,".Heatmap.ZscoreDev.pdf")); 
  pheatmap(z.scores,
           #annotation_col=colData(dev)[,1:2],
           cex=0.5,color=brewer_yes,scale="none",
           cluster_rows=T, cluster_cols=T)
  dev.off()
  
  
}