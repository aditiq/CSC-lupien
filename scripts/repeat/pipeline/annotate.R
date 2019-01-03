args <- commandArgs(trailingOnly = TRUE)
ipfile <- as.character(args[1])
opfile <- as.character(args[2])

library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)

genecode.txdb<- makeTxDbFromGFF("/mnt/work1/users/lupiengroup/Projects/CommonData/gencode.v24.annotation.gtf", format="gtf", 
                                 organism = "Homo sapiens", chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene) )


runfunc=function(ipfile, opfile,...){
    bg <- annotatePeak(ipfile, 
             tssRegion=c(-3000, 3000), 
             TxDb=genecode.txdb, 
             level = "transcript", annoDb="org.Hs.eg.db",
             sameStrand = FALSE, ignoreOverlap = FALSE,
             ignoreDownstream = TRUE, overlap = "TSS")

    df <- data.frame(stringsAsFactors = F, chr=seqnames(bg@anno), start=start(bg@anno)-1, end=end(bg@anno), strand=strand(bg@anno), meta=elementMetadata(bg@anno)$annotation)
    write.table(df, file=opfile, quote=F, row.names=F, col.names=T, sep="\t" )
}

## Run only once for background generation
runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/Repeatanalysis/Background.merged.bed", "results/PCSC1/Repeatenrich/permutation/results/PFA.Consensus.Catalogue/backgroundfile.annotated")

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak", 
        "results/PCSC1/Repeatenrich/permutation/results/PFA.Consensus.Catalogue/queryfile.annotated")

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/LSCp.Consensus.Catalogue/queryfile.annotated")

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/GBM.Consensus.Catalogue/queryfile.annotated")       

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/PCSC1.Consensus.Catalogue/queryfile.annotated") 

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak", 
        "results/PCSC1/Repeatenrich/permutation/results/HF.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.Consensus/Brain.Consensus.sortBed", 
        "results/PCSC1/Repeatenrich/permutation/results/Brain.Consensus.Catalogue/queryfile.annotated")             

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/Consensus.ESC.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/hESC.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/LSC.neg.Consensus.Catalogue.narrowPeak", 
        "results/PCSC1/Repeatenrich/permutation/results/LSCn.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/LSC.bulk.Consensus.Catalogue.narrowPeak", 
        "results/PCSC1/Repeatenrich/permutation/results/LSCbulk.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/HSC.ConsensusSet.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/HSC.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.differentiated.ConsensusSet.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/Hemat.differentiated.Consensus.Catalogue/queryfile.annotated")    

runfunc("/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.Progenitor.ConsensusSet.bed", 
        "results/PCSC1/Repeatenrich/permutation/results/Hemat.Progenitor.Consensus.Catalogue/queryfile.annotated")    

        
