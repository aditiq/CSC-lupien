#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor cluster
# Objective : Run fisher test for repeat enrichment
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

files <- dir("results/PCSC1/Repeatenrich/fishertest/" , pattern=".jaccard.bed")

function(name,...){

  a <- fread(paste0("results/PCSC1/Repeatenrich/fishertest/a.", name,".jaccard.bed"), header=F, sep=" ", stringsAsFactors=F, data.table=F)
  c <- fread(paste0("results/PCSC1/Repeatenrich/fishertest/b.", name,".jaccard.bed"), header=F, sep=" ", stringsAsFactors=F, data.table=F)
  b <- fread(paste0("results/PCSC1/Repeatenrich/fishertest/c.", name,".jaccard.bed"), header=F, sep=" ", stringsAsFactors=F, data.table=F)
  d <- fread(paste0("results/PCSC1/Repeatenrich/fishertest/d.", name,".jaccard.bed"), header=F, sep=" ", stringsAsFactors=F, data.table=F)
  colnames(a)[1] <- "a"
  colnames(b)[1] <- "b"
  colnames(c)[1] <- "c"
  colnames(d)[1] <- "d"

  pval <- c()

  for (f in 1:nrow(a)) {
    pval1 <- fisher.test(matrix(c(round(a[f,1]*10000000000),
                                  round(b[f,1]*10000000000),
                                  round(c[f,1]*10000000000),
                                  round(d[f,1]*10000000000)), 2, 2), alternative='greater')$p.value
    pval <- c(pval, pval1)
  }


}


#
