## Upset plot for KitchenSink2


library(UpSetR)
library(data.table)


#----------------------------------------------------------------
## UPset Plot
#----------------------------------------------------------------

binarymat <- fread("data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t",data.table=F, stringsAsFactors=F, check.names=F)
binarymat[,1] <- paste(binarymat[,1], binarymat[,2], binarymat[,3],sep="_")
binarymat$start <- NULL
binarymat$end <- NULL
colnames(binarymat) <- gsub("_peaks","", colnames(binarymat))


mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(mapping) <- mapping$sample
mapping <- mapping[colnames(binarymat)[2:ncol(binarymat)],]

binarymat$LSCp <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="pos.LSC")])>0,1,0)
binarymat$LSCn <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="neg.LSC")])>0,1,0)
binarymat$LSCb <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="Bulk.LSC")])>0,1,0)
binarymat$HSC.Hemat <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="HSC.Hemat")])>0,1,0)
binarymat$diff.Hemat <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="diff.Hemat")])>0,1,0)
binarymat$prog.Hemat <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1=="prog.Hemat")])>0,1,0)
binarymat$Brain <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1 %in% c("brain"))])>0,1,0)
binarymat$HF <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1 %in% c("HF"))])>0,1,0)
binarymat$PFA <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1 %in% c("PFA.pos"))])>0,1,0)
binarymat$GBM <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1 %in% c("GBM.pos"))])>0,1,0)
binarymat$ESC <- ifelse(rowSums(binarymat[,subset(mapping$sample, mapping$group1 %in% c("esc"))])>0,1,0)


pdf("results/KitchenSink2/UpsetPlot.KitchenSink2.pdf", useDingbats = F) ;
upset(binarymat[,c(1,210:220)],nsets=11, order.by = "freq")
dev.off()

pdf("results/KitchenSink2/UpsetPlot.KitchenSink2.Blood.pdf", useDingbats = F) ;
upset(binarymat[,c(1,210:215, 216)],nsets=11, order.by = "freq")
dev.off()

pdf("results/KitchenSink2/UpsetPlot.KitchenSink2.Brain.pdf", useDingbats = F) ;
upset(binarymat[,c(1,216:220)],nsets=11, order.by = "freq")
dev.off()
