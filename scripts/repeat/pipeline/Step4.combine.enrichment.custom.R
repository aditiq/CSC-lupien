#--------------------------------------------------------------
## Script to combine enrichment results from all analysis
#--------------------------------------------------------------

#--------------------------------------------------------------
# Load dependencies
#--------------------------------------------------------------

library(stringr)
library(gplots)
library(RColorBrewer)
library(pheatmap)
require(ggplot2)
library(stringr)

hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(100)
col2 <- brewer.pal(n = 8, name = "Reds")


files <- dir("results/PCSC1/Repeatenrich/permutation/results/", pattern="Enrichment", full.names=T, recursive=T)
#files <- files[-9]

## removing LSCn for now

repeats.l <-list()
number_reps <-vector()

for (i in 1:length(files)) {
    name <- gsub(".Consensus.Catalogue", "", str_split_fixed(files[i],"/", 9)[7])
    a <- read.table(files[i], header=T, sep="\t", stringsAsFactors=F)
    a2 <- a[,c("name", "qvalue")]
    colnames(a2)[2] <- paste0(name,".", colnames(a2)[2])
    repeats.l[[i]] <- a2
}

names(repeats.l) <- gsub("results/PCSC1/Repeatenrich/permutation/results//","",dirname(files))

repeats.df<-Reduce(function(x,y) merge (x,y,all=TRUE),repeats.l)
repeats.df$name <- gsub(".sorted", "", repeats.df$name)

write.table(repeats.df, file="results/PCSC1/Repeatenrich/permutation/results/Combined.enrichment.results.bed", sep="\t", quote=F, row.names=F, col.names=T)

repeats.df2 <- repeats.df[,2:ncol(repeats.df)]
rownames(repeats.df2) <- repeats.df[,1]
repeats.df2$PCSC1.qvalue <- NULL

colnames(repeats.df2) <- gsub(".qvalue","",colnames(repeats.df2) )

repeats.df3 <- subset(repeats.df2, repeats.df2$Brain >=0.01 & repeats.df2$Hemat.differentiated >=0.01 & 
                                      repeats.df2$PFA <0.01 & repeats.df2$GBM <0.01 & repeats.df2$LSCp <0.01 )

repeats.df3[repeats.df3 >= 0.01] <- 100
repeats.df3[repeats.df3 < 0.01] <- 1
repeats.df3[repeats.df3==100] <- 0

pdf("results/PCSC1/Repeatenrich/permutation/results/Heatmap.pdf")
pheatmap(as.matrix(repeats.df3[,c("GBM","PFA","LSCp", "HF","hESC","HSC","Hemat.Progenitor", "LSCbulk","LSCn","Hemat.differentiated" ,"Brain")]),
        cluster_cols = F,
        cluster_rows = T,
        clustering_method="ward.D2",
        clustering_distance_rows = "binary",
        #fontsize_row=2,
        color=(scalered),
        border_color="black",
        trace="none",scale="none")
dev.off()
