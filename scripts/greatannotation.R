### R-3.4.1
### mordor
### Objective : generate annotation from GREAT API results

# ==============================================================================
# load dependencies
# ==============================================================================
library(ggplot2)

# ==============================================================================
# Run analysis
# ==============================================================================
anno <- read.delim("results/PCSC1/Cluster/ExtractClusters.Flexclust/GREAT/Compiled.Annotation.txt", header=T, sep="\t", stringsAsFactors=F)

library(ggplot2)

for ( f in names(table(anno$Group))) {
  
  df <- subset(anno, anno$Group==f & anno$Category=="GO Biological Process" )
  df <- df[order(-df$NegLog10FDR),]
  df <- df[1:5,]
  
  p <- ggplot(data=df, aes(x=reorder(Term,NegLog10FDR), y=NegLog10FDR)) +
    geom_bar(stat="identity") + coord_flip()
  
  pdf(paste0("results/PCSC1/Cluster/ExtractClusters.Flexclust/GREAT/Barplot.",f,".pdf"))
  print(p)
  dev.off()

}