#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.5.0
# mordor
# Objective : Differential test using VGAM binomial from package monocle. Similar to scABC but also provides a qvalue
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
binaryds <- as.character(args[1])
groupinfofile <- as.character(args[2])
opname <- as.character(args[3])


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
suppressMessages(library(data.table))
suppressMessages(library(monocle))
suppressMessages(library("Matrix"))
suppressMessages(library(proxy))
suppressMessages(library(gplots))
suppressMessages(library(Rtsne))
suppressMessages(library(densityClust))



a <- read.table(binaryds, sep="\t", stringsAsFactors=F, header=T, check.names=F)
colnames(a) <- gsub("_peaks","", colnames(a))
rownames(a) <- paste(a[,1], a[,2],a[,3], sep="_")
mat <- as(as.matrix(a[,4:ncol(a)]), "dgCMatrix")
rownames(mat) = NULL
colnames(mat) = NULL

groupinfo <- read.table(groupinfofile, sep="\t", header=F, stringsAsFactors=F)
rownames(groupinfo) <- groupinfo[,1]
groupinfo <- groupinfo[colnames(a[,4:ncol(a)]),]

## phenotype data frame
pda = data.frame(name=colnames(a[,4:ncol(a)]), size=colSums(a[,4:ncol(a)]),
                 group=factor(groupinfo[,2]),tissue=factor(groupinfo[,3]),stringsAsFactors=F)
names(pda) = c("Name", "ReadDepth","Group","Tissue")
pda = new("AnnotatedDataFrame", data = pda)

# We next set up a “feature” data frame - a framework for storing information about data features (in this case, sites).

fda = as.data.frame(as.character(rownames(a)))
names(fda) = "Peak"
fda = new("AnnotatedDataFrame", data = fda)

## Convert binary matrix to cds
submat_cds =  newCellDataSet(mat,
                              featureData = fda,
                              phenoData = pda,
                              expressionFamily=binomialff(),
                              lowerDetectionLimit=1)

pData(submat_cds)$Size_Factor = 1
differtest = differentialGeneTest(submat_cds, fullModelFormulaStr = "~Group + ReadDepth + Tissue",  reducedModelFormulaStr = "~ReadDepth+Tissue", cores=1)
write.table(differtest, file=opname, sep="\t", row.names=F, col.names=T)
