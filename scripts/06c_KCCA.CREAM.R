#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Plot clusters from KCCA
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(stringr)

bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
hmcols = colorpanel(100, "steelblue", "white", "tomato")
#-----------------------------------
# Load binary matrix
#-----------------------------------

plothm=function(name, ...){

        mat <- read.table(paste0("results/PCSC1/CREAM/", name, ".Binarymat.txt"),check.names=F, stringsAsFactors = F, header=T, sep="\t")
        rownames(mat) <- paste(mat[,1], mat[,2], mat[,3], sep="_")

        #-----------------------------------
        # Load KCCA dataset
        #-----------------------------------
        load(paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.", name, ".100.Rdata"))
        load(paste0("results/PCSC1/Cluster/KCCA.Flexclust/kccadata/KCCA.", name, ".rownames.100.Rdata"))

        flexclust.clus <- as.data.frame(kcca.cl$cluster) ## cluster assignment
        rownames(flexclust.clus) <- rownames.obj

        opname=paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/", name, ".Flexclust.txt")
        flexclust.clus2 <- cbind(as.data.frame(str_split_fixed(rownames(flexclust.clus), "_", 3)), flexclust.clus[,1])
        write.table(flexclust.clus2, row.names=F, col.names=F, quote=F, sep="\t", file=opname )


        cc <- as.matrix(kcca.cl$centers)
        colnames(cc) <- gsub("_peaks","",colnames(mat)[4:ncol(mat)])
        rownames(cc) <- seq(1, max(kcca.cl$cluster),1)

        mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
        rownames(mapping) <- mapping$sample
        mapping <- mapping[colnames(cc),]
        mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,2)
        mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,2))
        colnames(cc) <- mapping$sample
        annocol <- data.frame(stringsAsFactors=F, stem=mapping$stem, tissue=mapping$tissue)
        rownames(annocol) <- colnames(cc)

        pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/", name, ".Flexclust.pdf"))
        pheatmap(as.matrix(cc),
                col=scalered, scale="none",
                cluster_rows=TRUE, cluster_cols=TRUE,
                clustering_method="ward.D2",
                fontsize_row=3,
                fontsize_col=3,
                annotation_col=annocol,
                clustering_distance_rows = "correlation",
                trace="none")
        dev.off()

        annocol <- annocol[order(annocol$stem),]

        pdf(paste0("results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters/", name, ".Flexclust.v2.pdf"))
        pheatmap(as.matrix(cc[,rownames(annocol)]),
                col=scalered, scale="none",
                cluster_rows=TRUE, cluster_cols=FALSE,
                clustering_method="ward.D2",
                fontsize_row=3,
                fontsize_col=3,
                annotation_col=annocol,
                clustering_distance_rows = "correlation",
                trace="none")
        dev.off()

}

plothm("KitchenSink2.CREAM")
plothm("KitchenSink1.CREAM")
plothm("PCSC1.CREAM")
