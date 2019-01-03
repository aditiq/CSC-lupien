### R-3.4.1
### mordor
### Objective : Examine ATAC landscape (1/0 binary) at regions overlapping immune genes TSS

#---------------------------------
# load dependencies
#---------------------------------
library(data.table)
library(pheatmap)
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))

bed_to_granges <- function(file){
  df <- read.table(file,header=F,sep="\t",stringsAsFactors=F)
  
  if(length(df) > 3){ df <- df[,c(1:3)] }
  
  if(length(df)<3){ stop("File has less than 3 columns") }
  
  header <- c('chr','start','end')
  names(df) <- header[1:length(names(df))]
  gr <- with(df, GRanges(chr, IRanges(start, end)))
  return(gr)
}

#---------------------------------------------
# Get immune genes TSS
#---------------------------------------------
query <- bed_to_granges("data/immunelandscape/TSSofImmuneGenes.bed")
queryDF <- data.frame(query)

#---------------------------------------------
# Get immune tss overlap
#---------------------------------------------
totalOverlap <- data.frame(seqnames = queryDF$seqnames, start = queryDF$start, end = queryDF$end)    

files <- dir("data/immunelandscape/", pattern="intersectu", full.names=T)
cellName = basename(files)

for (i in 1:length(files)){
  
  subject = bed_to_granges(files[i])
  hits = findOverlaps(query, subject)
  hitsDF <- data.frame(hits)
  cellName[i] <- gsub(".ImmuneTSS.intersectu.bed", '', cellName[i])
  totalOverlap[hitsDF$queryHits, cellName[i]] <- 1
  totalOverlap[-hitsDF$queryHits, cellName[i]] <- 0
}



#---------------------------------------------
# Get gene annotation
#---------------------------------------------
immunegene <- read.table("data/immunelandscape/TSSofImmuneGenes.bed", header=F, sep="\t", stringsAsFactors=F)
immunegene$id <- paste(immunegene$V1, immunegene$V2, immunegene$V3, sep="_")
immunegene2 <- immunegene[,c(14,6,7)]
immunegene2 <- immunegene2[!duplicated(immunegene2),]

totalOverlap$id <- paste(totalOverlap$seqnames, totalOverlap$start, totalOverlap$end, sep="_")
totalOverlap2 <- merge(totalOverlap, immunegene2, by.x="id", by.y="id")
totalOverlap2 <- totalOverlap2[!duplicated(totalOverlap2), ]
write.table(totalOverlap2, "results/PCSC1/immune.landscape/ImmuneTSSOverlap.txt", row.names=FALSE, col.names=T,sep="\t", quote=F)    

rownames(totalOverlap2) <- paste0(totalOverlap2$V7, "_", totalOverlap2$V6)
totalOverlap3 <- subset(totalOverlap2, rowSums(totalOverlap2[,5:99]) >0 )
pdf("results/PCSC1/immune.landscape/Heatmap.Immunelandscape.Binary.pdf")
pheatmap(as.matrix(totalOverlap3[,5:99]), col=scalered,
         clustering_rows=T, clustering_cols=T,
         clustering_method="ward.D2", fontsize_row=0.5,
         clustering_distance_rows="binary",
         clustering_distance_cols="binary")
dev.off()



png("results/PCSC1/immune.landscape/Heatmap.Immunelandscape.Binary.png")
pheatmap(as.matrix(totalOverlap3[,5:99]), col=scalered,
         clustering_rows=T, clustering_cols=T,
         clustering_method="ward.D2", fontsize_row=0.5,
         clustering_distance_rows="binary",
         clustering_distance_cols="binary")
dev.off()
