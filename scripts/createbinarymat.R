### R-3.4.1
### mordor
### Objective : Generate binary matrix of peaks for PCSC1.Consensus.Catalogue for downstream analysis

## Implemented frm https://github.com/ManchesterBioinference/Scasat/blob/master/Deconvoluting_cell-types_Scasat_functions.ipynb


args <- commandArgs(trailingOnly = TRUE)

if (length(args)==5) {

  mergedPeakFile <- as.character(args[1])
  peakFolder <- as.character(args[2])
  peakFilePattern <- as.character(args[3])
  outputFolder <- as.character(args[4])
  outputPeakFileName <- as.character(args[5])

} else if (length(args)==0) {
  stop("Please provide the following arguments - Merged peak file, Directory containing all peak files, pattern to identify peak files,
       Output dir and Binary matrix file name", call.=FALSE)
}


######################################
### load dependencies
######################################
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicRanges))

######################################
### Create function
######################################

## from https://davetang.org/muse/2015/02/04/bed-granges/
bed_to_granges <- function(file){
  df <- read.table(file,header=F,sep="\t",stringsAsFactors=F)

  if(length(df) > 3){ df <- df[,-c(4:length(df))] }

  if(length(df)<3){ stop("File has less than 3 columns") }

  header <- c('chr','start','end')
  names(df) <- header[1:length(names(df))]
  gr <- with(df, GRanges(chr, IRanges(start, end)))
  return(gr)
}

files = dir(path=peakFolder,  pattern = peakFilePattern, full.names=TRUE)
query <- bed_to_granges(mergedPeakFile)
queryDF <- data.frame(query)
totalOverlap <- data.frame(seqnames = queryDF$seqnames, start = queryDF$start, end = queryDF$end)
cellName = basename(files)

for (i in 1:length(files)){
  print(i)
  subject = bed_to_granges(files[i])
  hits = findOverlaps(query, subject)
  print("overlapfinding done")
  hitsDF <- data.frame(hits)
  cellName[i] <- gsub(peakFilePattern, '', cellName[i])
  totalOverlap[hitsDF$queryHits, cellName[i]] <- 1
  totalOverlap[-hitsDF$queryHits, cellName[i]] <- 0
}

outputFile = paste0(outputFolder,"/",outputPeakFileName)
write.table(totalOverlap, outputFile, row.names=FALSE, col.names=T,sep="\t", quote=F)
