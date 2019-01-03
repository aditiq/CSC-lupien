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
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)

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
    cisTopicObject <- runModels(cisTopicObject, topic=c(2, 5, 10, 15, 20, 25), seed=987, nCores=1, burnin = 250, iterations = 500, returnType="selectedModel") ## build the model
    saveRDS(cisTopicObject, file="results/KitchenSink2/Cluster/cisTopic/cisTopicObject.rds")
    cto <- readRDS("results/KitchenSink2/Cluster/cisTopic/cisTopicObject.rds")

    #cisTopicObject <- selectModel(cto)
    #logLikelihoodByIter(cto)

    ## Interpreting the models

    # A. Identification of cell states using the cell-cisTopic distributions

    ## 1. Add more metadata
    mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
    rownames(mapping) <- mapping$sample
    mapping$stem <- ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2)
    mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2))

    mapping2 <- mapping[colnames(ncounts),c("tissue","stem","group1")]

    cto@cell.data <- mapping2
    cto <- runtSNE(cto, perplexity=10, seed=987)
    cto <- runDM(cto)
    cto <- runPCA(cto)

    pdf("results/KitchenSink2/Cluster/cisTopic/PCA.pdf"); plotCellStates(cto,dim=2,col.low = "dodgerblue4", col.mid = "floralwhite", col.high = "brown2",  method='PCA', topic_contr='Zscore', topics='all', colorBy=c('group1', 'tissue','stem')) ; dev.off()
    #pdf("chk2.pdf"); plotCellStates(cto, method='Biplot', topic_contr='Zscore', topics = c(4,7,8), colorBy=c('group1', 'tissue'))  ; dev.off()
    pdf("results/KitchenSink2/Cluster/cisTopic/Heatmap.pdf"); cellTopicHeatmap(cto, col.low = "dodgerblue4", col.mid = "floralwhite", col.high = "brown2", colorBy=c('group1', 'tissue','stem'))  ; dev.off()



}
