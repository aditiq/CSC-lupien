### R-3.4.1
### mordor
### not run
### Objective : Fragment length distribution of ATAC-seq


######################################
### load dependencies
######################################
library(ATACseqQC)

######################################
### Run analysis
######################################

## load bam files

bamfile <- dir("data/bams", pattern="", recursive = T, include.dirs = T, full.names = T)
bamfile.labels <- gsub(".bam", "", basename(bamfile))

## estimate fragment size
fragSize <- fragSizeDist(bamfile, bamfile.labels)
