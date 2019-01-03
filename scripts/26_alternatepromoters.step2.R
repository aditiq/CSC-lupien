library(plyr)
library(data.table)
library(stringr)


binarymat <- fread("results/KitchenSink2/alternatepromoters//KitchenSink2.GencodeTSS.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
binarymat$id <- paste(binarymat[,1],binarymat[,2],binarymat[,3], sep="_")

mapping <- read.delim("data/ConsensusSet/KitchenSink1/Kitchensinkmapping.txt", header=T, sep="\t", stringsAsFactors=F)
'%!in%' <- function(x,y)!('%in%'(x,y))
mapping  <- subset(mapping, (mapping$sample %!in% c("ENCFF114WRE","ENCFF682RLN")))

brainnames <- subset(mapping$sample, mapping$group1 %in% c("brain"))
gbmnames <- subset(mapping$sample, mapping$group1 %in% c("GBM.pos"))
pfanames <- subset(mapping$sample, mapping$group1 %in% c("PFA.pos"))
lscpnames <- subset(mapping$sample, mapping$group1 %in% c("pos.LSC"))
hematdiffnames <- subset(mapping$sample, mapping$group1 %in% c("diff.Hemat"))

colnames(binarymat) <- gsub("_peaks","", colnames(binarymat))
gbm <- binarymat[,c("id",brainnames, gbmnames)]
pfa <- binarymat[,c("id",brainnames,pfanames )]
lscp <- binarymat[,c("id",hematdiffnames, lscpnames)]

gbm$csc <- rowSums(gbm[,gbmnames])
gbm$diff <- rowSums(gbm[,brainnames])
pfa$csc <- rowSums(pfa[,pfanames])
pfa$diff <- rowSums(pfa[,brainnames])
lscp$csc <- rowSums(lscp[,lscpnames])
lscp$diff <- rowSums(lscp[,hematdiffnames])

gbm <- gbm[,c("id","csc","diff")]
pfa <- pfa[,c("id","csc","diff")]
lscp <- lscp[,c("id","csc","diff")]

tss <- read.table("../common.data/c3d.anchors.gencodev24.500bparoundTSS.bed", header=F, sep="\t", stringsAsFactors=F)
tss$id <- paste(tss[,1],tss[,2],tss[,3], sep="_")

ap_function <- function(a2,opname,cscno,diffno,...){

    tss2 <- merge(tss, a2, by.x="id", by.y="id", all.x=T)
    tss2$csc2 <- ifelse(tss2$csc >0,1,0)
    tss2$diff2 <- ifelse(tss2$diff >0,1,0)

    tss2$gene <- data.frame(str_split_fixed(tss2$V5, "_",2), stringsAsFactors=F)[1]$X1
    tss2 <- tss2[order(tss2$gene, tss2$V5),]

    ## remove zero-zero and one-zero rows
    tss2$chk <- paste0(tss2$csc2, "_", tss2$diff2)
    tss3 <- rbind(subset(tss2, tss2$chk=="1_0" & tss2$csc >=cscno),  subset(tss2,tss2$chk=="1_1" & tss2$csc >=cscno & tss2$diff >=diffno))

    ## Remove single isoform genes
    ## count no. of isoforms per gene

    df2 <- tss3[,c("gene","V5")]
    df3 <- count(df2, c('gene'))
    df4 <- subset(df3, df3$freq >1)
    tss4 <- subset(tss3, tss3$gene %in% df4$gene)

    ## count no. of unique chk ids for each gene
    df3 <- ddply(tss4,~gene,summarise,number_of_distinct_ids=length(unique(chk)))
    df4 <- subset(df3, df3$number_of_distinct_ids >1)
    tss5 <- subset(tss4, tss4$gene %in% df4$gene)
    tss5 <- tss5[!duplicated(tss5),]
    tss5$V1 <- NULL
    tss5$V2 <- NULL
    tss5$V3 <- NULL
    write.table(names(table(tss5$gene)), file=paste0("results/KitchenSink2/alternatepromoters/", opname,".AP.txt"), sep="\t", row.names=F, col.names=F, quote=F)

}


ap_function(gbm,"GBM",11,9)
ap_function(pfa,"PFA",3,9)
ap_function(lscp,"LSCp", 9,16)
