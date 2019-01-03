#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Cluster specific peaks obtained from scABC maybe shared. So assigning a priority system for it and extracting separate files
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------
# load dependencies
#----------------------------------------------------------------
library(gplots)
library(RColorBrewer)
library(pheatmap)
require(ggplot2)
library(reshape)
library(stringr)
library(data.table)


hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
#----------------------------------------------------------------
# Get scABC data
#----------------------------------------------------------------

scabc_assign=function(ds, opname,...){

  scabc <- read.table(ds, header=T, sep="\t", stringsAsFactors = F)
  scabc$GBM.grp <- ifelse(scabc$GBM<1e-6,1,0)
  scabc$LSC.grp <- ifelse(scabc$LSC<1e-6,1,0)
  scabc$PFA.grp <- ifelse(scabc$PFA<1e-6,1,0)
  scabc$GBM.PFA.grp <- ifelse(scabc$GBM.PFA<1e-6,1,0)
  scabc$LSC.GBM.grp <- ifelse(scabc$LSC.GBM<1e-6,1,0)
  scabc$LSC.PFA.grp <- ifelse(scabc$LSC.PFA<1e-6,1,0)

  scabc$combined <-  ifelse(scabc$GBM.grp>0,"GBM",
                                  ifelse(scabc$LSC.grp >0,"LSC", 
                                         ifelse(scabc$PFA.grp>0,"PFA", 
                                                ifelse( (scabc$GBM.PFA.grp >0 | scabc$LSC.GBM.grp >0 | scabc$LSC.PFA.grp >0 ) , "Shared","Else")))) 

  write.table(scabc, file=paste0("results/PCSC1/Cluster/scABC/scABC.Combined.", opname,".WdAnnotation.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  for (f in names(table(scabc$combined))){
    
    write.table(as.data.frame((str_split_fixed(subset(scabc$id,scabc$combined==f), "_",3))), 
                file=paste0("results/PCSC1/Cluster/scABC/", f, ".", opname, ".p0.05.bed"), 
                row.names=F, col.names=F, sep="\t", quote=F)
  }
}

# Run function on both enhancer and promoter sets
scabc_assign("results/PCSC1/Cluster/scABC/scABC.Combined.Promoter.txt", opname="Promoter")
scabc_assign("results/PCSC1/Cluster/scABC/scABC.Combined.Enhancer.txt", opname="Enhancer")



