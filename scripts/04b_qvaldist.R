### R-3.4.1
### mordor
### Objective : plot qval distribution for the catalogue of ATAC seq peaks

#=======================================================
# Load dependencies
#=======================================================
library(data.table)
library(ggplot2)

#=======================================================
# Read in data
#=======================================================

files=dir("data/ConsensusSet/PCSC1/", pattern=".Catalogue.narrowPeak", full.names=F)
files <- files[-grep("Consensus", files)]

for ( f in 1:length(files)){
  
  dat <- fread(paste0("data/ConsensusSet/PCSC1/",files[f]), sep="\t", stringsAsFactors = F, header=F)
  dat$qval <- 10^(-1*dat$V9)
  p <- ggplot(dat, aes(x=V9)) + geom_histogram(bins = 100) + scale_x_continuous(breaks=seq(0, max(dat$V9), 5))
  pdf(paste0("data/ConsensusSet/PCSC1/Density",basename(files[f]), ".pdf"))
  print(p)
  dev.off()

}
