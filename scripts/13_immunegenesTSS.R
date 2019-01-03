### R-3.4.1
### mordor
### Objective : Get the TSS at immune genes to further check overlap with ATAC peaks

#---------------------------------
# load dependencies
#---------------------------------
library(data.table)

#---------------------------------
# Load immune genes
#---------------------------------
immune <- fread("data/nanostringimmunepanel.txt", header=T, sep="\t", stringsAsFactors = F, data.table=F)

#---------------------------------
# Load TSS
#---------------------------------
tss <- fread("../common.data/gencode.v24.annotation.TSS.wdgeneinfo.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)
colnames(tss) <- c("TSS.chr", "TSS.start", "TSS.stop","strand","geneid","txtid","genename","random")
tss$random <- NULL

#---------------------------------
# get TSS for immune genes
#---------------------------------
immune.tss <- merge(immune, tss, by.x="GeneName", by.y="genename", all.x=T)
immune.tss.na <- subset(immune.tss, is.na(immune.tss$txtid) == TRUE)
immune.tss.na2 <- merge(immune.tss.na[,1:7], tss, by.x="AltName", by.y="genename", all.x=T)

immune.tss2 <- rbind(subset(immune.tss, is.na(immune.tss$txtid) == FALSE), subset(immune.tss.na2, is.na(immune.tss.na2$txtid)== FALSE))
write.table(immune.tss2[,c(8:13,1:7)], file="data/immunelandscape/TSSofImmuneGenes.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(immune.tss2[,c(8:13,1:7)], file="data/immunelandscape/TSSofImmuneGenes.bed", sep="\t", row.names=F, col.names=F, quote=F)

