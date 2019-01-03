
## Step1 -- make binarymat for shared CSC regions
#DIR="results/PCSC1/Cluster/KCCA.Flexclust/ExtractedClusters" ;
#Rscript scripts/createbinarymat.R $DIR/Shared.PCSC1.notHK.bed data/ConsensusSet/KitchenSink2/temp.peaks/ ".narrowPeak" ./ Shared.PCSC1.notHK.KS2.Binarymat.txt


library(data.table)
library(pheatmap)
library(RColorBrewer)


colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
brewer_celsius = c('#313695', '#5083BB', '#8FC3DD', '#D2ECF4', '#FFFFBF', '#FDD384', '#F88D51', '#DE3F2E', '#A50026')
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)

binarymat <- fread("Shared.PCSC1.notHK.KS2.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F, data.table=F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")
colnames(binarymat) <- gsub("_peaks","", colnames(binarymat))
mat <- binarymat[,4:ncol(binarymat)]

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[colnames(mat),]
mapping$stem <- ifelse(grepl("\\.pos", (mapping$group1))==TRUE, 1,0)
mapping$stem <- ifelse( mapping$stem==1, 1, ifelse(grepl("pos\\.", (mapping$group1))==TRUE, 1,0))

mapping <- mapping[order(mapping$group1),]
mat2 <- mat[,rownames(mapping)]
colnames(mat2) <- mapping$group1

for ( f in 1:length( unique(mapping$group1))) {

  mat2$temp <- apply(mat2[,grepl(unique(mapping$group1)[f], colnames(mat2))] , 1, sum)
  colnames(mat2)[ncol(mat2)] <- paste0("Sum.",unique(mapping$group1)[f] )
}

mat3 <- mat2[,grepl("Sum.", colnames(mat2))]
mat3[mat3>=1] <- 1

pdf("Shared.PCSC1.notHK.KS2.Heatmap.pdf", useDingbats=F)
pheatmap(as.matrix(mat3),
        col=scalered, scale="none",
        #fontsize_row=1.5,fontsize_col=1.5,
        cluster_rows=F, cluster_cols=F,
        clustering_method="ward.D2",
        clustering_distance_rows = "correlation",
        #annotation_col=annocol,
        trace="none")
dev.off()
