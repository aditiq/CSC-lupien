#--------------------------------------------------------------
# Objective : Enrichment test2ing of repeat families ((No. of times Expected >= X)/N)
# mordor
#--------------------------------------------------------------


#--------------------------------------------------------------
# Load dependencies
#--------------------------------------------------------------
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
queryfilepath <- as.character(args[1])
jaccardfilesdir <- as.character(args[2])
opdir <- as.character(args[3])
#--------------------------------------------------------------
# Load query file jaccard results
#--------------------------------------------------------------
queryfile <- read.table(queryfilepath, header=F, sep=" ", stringsAsFactors=F)

#--------------------------------------------------------------
# Load background file jaccard results
#--------------------------------------------------------------
bgfiles <- dir(jaccardfilesdir, pattern="Bg.Jaccard.", full.names=T)

#--------------------------------------------------------------
# for each bg file, get pvalue
#--------------------------------------------------------------

enrich <- as.data.frame(str_split_fixed(gsub("Bg.Jaccard.","",gsub(".bed","",basename(bgfiles))),"_",2))
colnames(enrich) <- c("number","name")
enrich <- merge(enrich,queryfile, by.x="number",by.y="V2")
colnames(enrich)[3] <- "Query.jaccard"

pval <- c()
for (f in 1:1544){
  name <- dir(jaccardfilesdir, pattern=paste0("Bg.Jaccard.",f,"_"), full.names=T)
  bg <- read.table(name, sep=" ", stringsAsFactors = F, header=F)
  pval <- c(pval,sum(ifelse(bg$V1 >= subset(enrich$Query.jaccard, enrich$number==f),1,0))/10000)
}

pval <- as.data.frame(pval)
pval$number <- seq(1,1544,1)
colnames(pval)[1] <- "P.value"
pval$qvalue <- p.adjust(pval$P.value, method = "BH")

#--------------------------------------------------------------
# Get final results dataset
#--------------------------------------------------------------
enrich <- merge(enrich, pval, by.x="number", by.y="number")
write.table(enrich, file=paste0(opdir,"Enrichment.results.bed"), sep="\t", quote=F,row.names=F, col.names=T)
