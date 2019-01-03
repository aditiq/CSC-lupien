### R-3.4.1
### mordor
### Objective : Split binary matrix

binarymat <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t", stringsAsFactors = F,check.names = F)
rownames(binarymat) <- paste(binarymat[,1],binarymat[,2],binarymat[,3],sep="_")

enh <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancers.bed", header=F, sep="\t", stringsAsFactors = F,check.names = F)
rownames(enh) <- paste(enh[,1],enh[,2],enh[,3],sep="_")

prom <- read.table("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoters.bed", header=F, sep="\t", stringsAsFactors = F,check.names = F)
rownames(prom) <- paste(prom[,1],prom[,2],prom[,3],sep="_")

enhmat <- (binarymat[rownames(enh),])
prommat <- (binarymat[rownames(prom),])

write.table(enhmat, file="data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt",row.names=F, col.names=T, sep="\t", quote=F)
write.table(prommat, file="data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Promoter.Binarymat.txt",row.names=F, col.names=T, sep="\t", quote=F)

