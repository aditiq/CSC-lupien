#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Run cistopic to identify CSC specific peaks
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
library(cisTopic)
library(pheatmap)
library("RColorBrewer")
library(gplots)

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026')
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")

#----------------------------------------------------------------
#Load binary data
#----------------------------------------------------------------

#runcistopic=function(ds, ...){
    binarymatds="data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt"

    binarymat <- fread(binarymatds, header=T, sep="\t", stringsAsFactors = F,check.names = F, data.table=F)
    rownames(binarymat) <- paste(binarymat[,1], binarymat[,2], binarymat[,3], sep="_")
    ncounts <- binarymat[,4:ncol(binarymat)]
    rownames(ncounts) <- paste0(binarymat[,1], ":", binarymat[,2], "-",binarymat[,3] )
    colnames(ncounts) <- gsub("_peaks","", colnames(ncounts))
    cisTopicObject <- createcisTopicObject(ncounts, project.name='cisTopic_CSC')

    mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
    rownames(mapping) <- mapping$sample
    mapping <- mapping[colnames(ncounts),c("group1","tissue")]

    ## Add metadata
    cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = mapping) ## add metadata

    ## Run the model
    cisTopicObject <- runModels(cisTopicObject, topic=c(25, 50, 75, 100), seed=987, nCores=1, burnin = 250, iterations = 500, returnType="allModels") ## build the model
    saveRDS(cisTopicObject, file="results/KitchenSink2/Cluster/cisTopic/cisTopicObject.all.rds")

    cto <- readRDS("results/KitchenSink2/Cluster/cisTopic/cisTopicObject.all.rds")

    pdf("results/KitchenSink2/Cluster/cisTopic/logLikelihoodByIter.All.pdf"); logLikelihoodByIter(cto) ; dev.off()
    #cto <- selectModel(cto)



    ## Interpreting the models

    # A. Identification of cell states using the cell-cisTopic distributions

    ## 1. Add more metadata
    mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
    rownames(mapping) <- mapping$sample
    mapping$stem <- ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2)
    mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2))
    mapping2 <- mapping[colnames(ncounts),c("tissue","stem","group1")]
    mapping2$group2 <- ifelse(grepl("LSC|Hemat", mapping2$group1),"blood",ifelse(grepl("esc", mapping2$group1),"esc","brain"))

    for (f in c(25, 50,75,100)){

      cto <- selectModel(cto, select=f)
      cto@cell.data <- mapping2
      cto <- runtSNE(cto, perplexity=10, seed=987)
      cto <- runDM(cto)
      cto <- runPCA(cto)

      pdf(paste0("results/KitchenSink2/Cluster/cisTopic/", f,".PCA.pdf")); plotCellStates(cto,dim=2,col.low = "dodgerblue4", col.mid = "floralwhite", col.high = "brown2",  method='PCA', topic_contr='Zscore', topics='all', colorBy=c('group1', 'tissue','stem')) ; dev.off()
      #pdf("chk2.pdf"); plotCellStates(cto, method='Biplot', topic_contr='Zscore', topics = c(4,7,8), colorBy=c('group1', 'tissue'))  ; dev.off()
      pdf(paste0("results/KitchenSink2/Cluster/cisTopic/", f,".Heatmap.pdf")); cellTopicHeatmap(cto, col.low = "dodgerblue4", col.mid = "floralwhite", col.high = "brown2", colorBy=c('group1', 'tissue','stem'))  ; dev.off()

      mapping3 <- mapping2[order(mapping2$stem),]
      df <- cto@selected.model$document_expects
      rownames(df) <- paste("Topic", seq(1, nrow(df)))
      colnames(df) <- cto@cell.names
      df <- df[,rownames(mapping3)]
      topic.mat <- scale(df,   center = TRUE, scale = TRUE)
      colorPal <- grDevices::colorRampPalette(c("dodgerblue4", "floralwhite", "brown2"))
      pdf(paste0("results/KitchenSink2/Cluster/cisTopic/", f,".pHeatmap.pdf"));
      pheatmap(  topic.mat,
                cluster_rows=T, cluster_cols=F,
                clustering_method="ward.D2", clustering_distance_rows="euclidean",
                col=colorPal(20), scale="none",
                annotation_col=mapping3)
      dev.off()

    }

    ## Extract common regions from each cluster 25,50,75,100
    cto <- selectModel(cto, select=50)
    cto <- getRegionsScores(cto, method='Zscore', scale=TRUE)
    cto <- binarizecisTopics(cto, thrP=0.99, plot=FALSE)
    getBedFiles(cto, path='results/KitchenSink2/Cluster/cisTopic/bedfiles50/')
