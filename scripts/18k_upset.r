#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(data.table)
library(UpSetR)

#----------------------------------------------------------------
# upset plot
#----------------------------------------------------------------
#
# gbm="results/PCSC1/cicero/run2.readcount/GBM.common.coaccessgt0.enhancer.bed"
# lsc="results/PCSC1/cicero/run2.readcount/LSCp.common.coaccessgt0.enhancer.bed"
# pfa="results/PCSC1/cicero/run2.readcount/PFA.common.coaccessgt0.enhancer.bed"
#
# cat $gbm $lsc $pfa | awk 'BEGIN {FS=OFS="\t"} { print $1,$2,$3}' - | sortBed -i stdin | mergeBed -i stdin >> tmp.bed
# mkdir tmp.cicero/
# cp $gbm   tmp.cicero/
# cp $lsc   tmp.cicero/
# cp $pfa   tmp.cicero/
#
# Rscript scripts/createbinarymat.R tmp.bed tmp.cicero/ ".bed" tmp.cicero/ common.coaccessgt0.enhancer.binarymat.txt

binarymat <- fread("tmp.cicero/common.coaccessgt0.enhancer.binarymat.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
colnames(binarymat) <- gsub(".common.coaccessgt0.enhancer","", colnames(binarymat))

pdf("results/PCSC1/cicero/run2.readcount/summary.plots/UpsetPlot.common.coaccessgt0.enhancer.pdf", useDingbats=F,onefile=F)
upset(binarymat[,-c(1:3)],order.by = "freq")
dev.off()
