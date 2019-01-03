### R-3.5.0
### mordor

# =====================================================================================================================================
# Objective : Compare C3D  Loops to those from GSCs HiC (SU2C) and OCIAML2 HiC (gbm and lsc resp)
# Only loops found overlapping HiC loops will be examined. (maxgap=1kb)
# Distance based analysis will be done (10kb ranges capped at 500kb)
# =====================================================================================================================================


# ========================================================================G======
library(data.table)
library(stringr)
library(reshape2)
library(cicero)

# ==============================================================================
# Run Analysis
# ==============================================================================
#
#
# ## Subset HiC for either end of the loop lying in a TSS
# cat ../common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr* > tmp2.txt
#
# loopfile="data/HiC/GBM/merged_loops_for_G523_aggr_ChIP_with_motifs.bedpe
# data/HiC/GBM/merged_loops_for_G583_aggr_ChIP_with_motifs.bedpe
# data/HiC/GBM/merged_loops_for_G567_aggr_ChIP_with_motifs.bedpe"
#
# for f in $loopfile ;
# do
#   name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/merged_loops_for_//g' | sed -e 's/_aggr_ChIP_with_motifs.bedpe//g')
#
#   awk 'BEGIN {FS=OFS="\t"} {  print "chr"$1,$2,$3}' $f | grep -v x1 | \
#   intersectBed -a stdin -b tmp2.txt -u -F 1 > results/PCSC1/C3D/$name.Hic.regions.in.TSS.bed
#
#   awk 'BEGIN {FS=OFS="\t"} {  print "chr"$4,$5,$6}' $f | grep -v y1 | \
#   intersectBed -a stdin -b tmp2.txt -u  -F 1  >> results/PCSC1/C3D/$name.Hic.regions.in.TSS.bed
#
# done


## Read in C3D Loops

c3d <- fread("results/PCSC1/C3D/finalresults/GBM.combined.results.txt", header=T, sep="\t", stringsAsFactors=F, data.table=F)
c3d <- subset(c3d, c3d[,1] !="COORD_1")
c3d <- c3d[,c(1,4,7,10)]
colnames(c3d) <- c("Peak1", "Peak2", "correlation","qvalue")
c3d$Peak1 <- gsub(":","_", gsub("-","_", c3d$Peak1))
c3d$Peak2 <- gsub(":","_", gsub("-","_", c3d$Peak2))

for  (name in c("G523","G583","G567")){

    ## Read in HiC data
    print(name)
    hicloops <- read.table(paste0("data/HiC/GBM/merged_loops_for_", name, "_aggr_ChIP_with_motifs.bedpe"), header=T, sep="\t", stringsAsFactors=F)
    hicloops$Peak1 <- paste(paste0("chr",hicloops[,1]), hicloops[,2], hicloops[,3], sep="_")
    hicloops$Peak2 <- paste(paste0("chr",hicloops[,4]), hicloops[,5], hicloops[,6], sep="_")
    hicloopsconns <- hicloops[,c("Peak1", "Peak2")]

    ##  Restrict HiC loops to those overlapping TSS
    hic.tss <- read.table(paste0("results/PCSC1/C3D/", name, ".Hic.regions.in.TSS.bed"), header=F, sep="\t", stringsAsFactors=F)
    hic.tss$id <- paste(hic.tss[,1], hic.tss[,2],hic.tss[,3], sep="_")
    hicloopsconns <- subset(hicloopsconns, hicloopsconns$Peak1 %in% hic.tss$id | hicloopsconns$Peak2 %in% hic.tss$id )

    ## Restrict HiC loops to those overlapping GBM ATAC peaks
    HicInAtac <- read.table(paste0("results/PCSC1/cicero/", name,".HiC.regions.in.", name, ".atac.bed"), header=F, sep="\t", stringsAsFactors=F)
    HicInAtac$id <- paste(HicInAtac[,1], HicInAtac[,2], HicInAtac[,3], sep="_")
    hicloopsconns <- subset(hicloopsconns, hicloopsconns$Peak1 %in% HicInAtac$id | hicloopsconns$Peak2 %in% HicInAtac$id )

    c3d2 <- c3d
    ## Overlap
    c3d2$in_hic <- compare_connections(c3d2, hicloopsconns, maxgap=1000)
    c3d2$corrflag <- as.character(cut(as.numeric(c3d2$correlation), seq(0,1, 0.1) ) )
    c3d2$corrflag <- ifelse(is.na(c3d2$corrflag )==TRUE, 0,c3d2$corrflag )
    c3d2$qvalflag <- ifelse(c3d2$qvalue < 0.05, "qvalsig","qvalnonsig")
    c3d2$sigflag <- ifelse(c3d2$qvalue < 0.05 & c3d2$correlation >=0.7, "sig","nonsig")

    freqdf <- as.data.frame(table(c3d2$in_hic,c3d2$corrflag,c3d2$qvalflag))
    freqdf.reshape <- dcast(data=freqdf, formula=Var2+Var1~Var3)
    write.table( freqdf.reshape, file=paste0("results/PCSC1/C3D/finalresults/GBM.", name, ".Hic.overlap.1kb.freqdf.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    ## Calculate percentage

    pctlist <- c()
    grplist <- c()
    grp2list <- c()
    for ( f in unique(  freqdf.reshape$Var2)){
      tmp <- subset(freqdf.reshape, freqdf.reshape$Var2==f)
      tmp <- t(tmp[,-c(1,2)])
      colnames(tmp) <- c("f","t")
      tmp <- as.data.frame(tmp)
      tmp$pct <- (tmp$t/(tmp$t + tmp$f))*100
      tmp$grp <- f
      colnames(tmp) <- paste0(colnames(tmp), ".",f)
      pctlist <- c(pctlist, tmp$pct)
      grplist <- c(grplist, tmp$grp)
      grp2list <- c(grp2list, rownames(tmp))
    }

    summ <- data.frame(stringsAsFactors=F, grplist, grp2list,pctlist)
    summ[is.na(summ)] <- 0
    summ2 <- dcast(summ, grplist~grp2list)
    write.table( summ2, file=paste0("results/PCSC1/C3D/finalresults/GBM.", name, ".Hic.overlap.1kb.dcast.Percentage.txt"), sep="\t", quote=F, row.names=F, col.names=T)


    ## Run for random regions
    cicerordsdedup2 <- c3d2
    cicerordsdedup2$chromflag <- substr(cicerordsdedup2$Peak1,1,5)

    combined.shuf=NULL
    for ( f in unique(cicerordsdedup2$chromflag)){
    chk <- subset(cicerordsdedup2, cicerordsdedup2$chromflag==f)
    for (k in 1){
          print(k)
          chk$end1shuf  <- chk[sample(nrow(chk)), "Peak1"]
          chk$end2shuf  <- chk[sample(nrow(chk)), "Peak2"]
          chk <- chk[,c("end1shuf","end2shuf", "in_hic", "corrflag", "qvalflag","sigflag")]
          colnames(chk)[1:2] <- c("Peak1","Peak2")
          chk$in_hic2 <- compare_connections(chk, hicloopsconns, maxgap=1000)
          combined.shuf <- rbind(combined.shuf, chk)
        }
    }


    freqdf <- as.data.frame(table(  combined.shuf$in_hic2,combined.shuf$corrflag,combined.shuf$qvalflag))
    freqdf.reshape <- dcast(data=freqdf, formula=Var2+Var1~Var3)

    ## Calculate percentage

    pctlist <- c()
    grplist <- c()
    grp2list <- c()
    for ( f in unique(  freqdf.reshape$Var2)){
      tmp <- subset(freqdf.reshape, freqdf.reshape$Var2==f)
      tmp <- t(tmp[,-c(1,2)])
      colnames(tmp) <- c("f","t")
      tmp <- as.data.frame(tmp)
      tmp$pct <- (tmp$t/(tmp$t + tmp$f))*100
      tmp$grp <- f
      colnames(tmp) <- paste0(colnames(tmp), ".",f)
      pctlist <- c(pctlist, tmp$pct)
      grplist <- c(grplist, tmp$grp)
      grp2list <- c(grp2list, rownames(tmp))
    }

    summ <- data.frame(stringsAsFactors=F, grplist, grp2list,pctlist)
    summ[is.na(summ)] <- 0
    summ2 <- dcast(summ, grplist~grp2list)
    write.table( summ2, file=paste0("results/PCSC1/C3D/finalresults/GBM.Shuffled.", name, ".Hic.overlap.1kb.dcast.Percentage.txt"), sep="\t", quote=F, row.names=F, col.names=T)
}


## Plot results of overlap with Hic
files <- dir("results/PCSC1/C3D/finalresults/", "dcast.Percentage.txt", full.names=T)
for (f in 1:length(files)){

  a <- read.table(files[f], sep="\t", stringsAsFactors=F, header=T)
  name <- gsub(".Hic.overlap.1kb.dcast.Percentage.txt", "", basename(files[f]))

  pdf(paste0("results/PCSC1/C3D/finalresults/" , name, ".pdf"))
  p <- ggplot(data=melt(a), aes(x=rep(seq(0,1,0.1),2), y=value, group=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  theme(text = element_text(size=12)) +
  labs(x = "Correlation" , y="% of C3D interactions overlapped by HiC")+
  scale_x_continuous(breaks=seq(0,1,0.1))+
  annotate("text",  x=Inf, y = Inf, label = gsub("Shuffled.","",gsub("GBM.","",name)), vjust=2, hjust=2)
  print(p)
  dev.off()

}
