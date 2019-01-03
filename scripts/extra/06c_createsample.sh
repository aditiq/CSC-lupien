#!/bin/bash
## NOT RUN

# a <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt", header=T, sep="\t", stringsAsFactors=F)
# a$lsc <- ifelse(rowSums(a[,4:21])>0,"lsc","zero")
# a$gbm <- ifelse(rowSums(a[,22:44])>0,"gbm","zero")
# a$pfa <- ifelse(rowSums(a[,45:51])>0,"pfa","zero")

# a$comb <- paste(a$lsc, a$gbm, a$pfa,sep="_")
# write.table(a[,c("seqnames","start","end", "comb")],
#             file="data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.State.txt",
#             sep="\t", row.names=F, col.names=T, quote=F)

module load bedops/2.4.14   

sh scripts/createsamplebed.v2.sh data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.State.txt 3330.18 data/ConsensusSet/PCSC1/sample "Enh"


