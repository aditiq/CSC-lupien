### R-3.5.0
### mordor

# =====================================================================================================================================
# Objective : Compare Cicero  Loops to those from GSCs HiC (SU2C) and OCIAML2 HiC (gbm and lsc resp)
# Only loops found overlapping HiC loops will be examined. (maxgap=1kb)
# Distance based analysis will be done (10kb ranges capped at 500kb)
# We also want to know no. of cicero loops with coaccess gt0 and le 0 overlapping and not overlapping HiC
# =====================================================================================================================================


# ==============================================================================
# load dependencies
# ==============================================================================
library(cicero)
library(data.table)
library(stringr)
library(reshape2)

# ==============================================================================
# Run Analysis
# ==============================================================================


summarisefunction=function(cicerodsinput, opname, ...){

    cicerodsinput$coaccess2 <- ifelse(  cicerodsinput$coaccess > 0.5, 0.5,cicerodsinput$coaccess)
    cicerodsinput$coaccessranged <- as.character(cut(cicerodsinput$coaccess2, seq(0,0.5, 0.1) ) )
    cicerodsinput$coaccessranged <-  ifelse( is.na(cicerodsinput$coaccessranged)==TRUE, "zeronegative", cicerodsinput$coaccessranged)

    p1 <- (as.numeric(str_split_fixed(cicerodsinput$Peak1,"_", 3)[,3]) + as.numeric(str_split_fixed(cicerodsinput$Peak1,"_", 3)[,2])) /2
    p2 <- (as.numeric(str_split_fixed(cicerodsinput$Peak2,"_", 3)[,3]) + as.numeric(str_split_fixed(cicerodsinput$Peak2,"_", 3)[,2])) /2

    cicerodsinput$distance <- (abs(p2-p1))/1000
    cicerodsinput$distancerange <-   as.character(cut(cicerodsinput$distance, seq(0, 300, 50) ) )
    cicerodsinput$distancerange <- ifelse(is.na(cicerodsinput$distancerange)==TRUE, "(300,above]", cicerodsinput$distancerange)
    freqdf <- as.data.frame(table(cicerodsinput$in_hic, cicerodsinput$coaccessranged, cicerodsinput$distancerange))
    freqdf.reshape <- dcast(data=freqdf, formula=Var2+Var1~Var3)
    write.table(freqdf, file=paste0("results/PCSC1/cicero/run2.readcount/",opname ,".Hic.overlap.1kb.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    write.table(  freqdf.reshape, file=paste0("results/PCSC1/cicero/run2.readcount/",opname, ".Hic.overlap.1kb.dcast.txt"), sep="\t", quote=F, row.names=F, col.names=T)

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
    write.table( summ2, file=paste0("results/PCSC1/cicero/run2.readcount/", opname, ".Hic.overlap.1kb.dcast.Percentage.txt"), sep="\t", quote=F, row.names=F, col.names=T)

}


comparehicfunction=function(name,...){

## Read in cicero Loops

 cicerords <- readRDS(paste0("results/PCSC1/cicero/run2.readcount/", name, ".connsfull.rds"))
 #cicerords$pairID <- paste(cicerords$Peak1, cicerords$Peak2, sep=":")

## Read in deduped cicero loops
  # dedupcon <- fread(paste0("results/PCSC1/cicero/run2.readcount/", name, ".deduped.connsfull.pairs.bed"), header=F, sep="\t", stringsAsFactors=F, data.table=F)
  # colnames(dedupcon) <- c("Peak1", "Peak2")
  # dedupcon$pairID <- paste(dedupcon$Peak1, dedupcon$Peak2, sep=":")

# ## Restrict cicero loops to deduped ones
#  cicerordsdedup <- subset(cicerords,cicerords$pairID %in% dedupcon$pairID)

  cicerordsdedup <- cicerords

  ## Keep HiC loops to those overlapping ATAC

  if ( name=="GBM" ){

    for (hicfile in c("G567","G523","G583"))
    {
          hicloops <- read.table(paste0("data/HiC/GBM/merged_loops_for_", hicfile, "_aggr_ChIP_with_motifs.bedpe"), header=T, sep="\t", stringsAsFactors=F)
          hicloops$Peak1 <- paste(paste0("chr",hicloops[,1]), hicloops[,2], hicloops[,3], sep="_")
          hicloops$Peak2 <- paste(paste0("chr",hicloops[,4]), hicloops[,5], hicloops[,6], sep="_")
          hicloopsconns <- hicloops[,c("Peak1", "Peak2")]

          HicInAtac <- read.table(paste0("results/PCSC1/cicero/", hicfile,".HiC.regions.in.", hicfile, ".atac.bed"), header=F, sep="\t", stringsAsFactors=F)
          HicInAtac$id <- paste(HicInAtac[,1], HicInAtac[,2], HicInAtac[,3], sep="_")
          hicloopsconns <- subset(hicloopsconns, hicloopsconns$Peak1 %in% HicInAtac$id | hicloopsconns$Peak2 %in% HicInAtac$id )

          ## Cicero loops found in HiC
          #cicerordsdedup2.zeroneg <- subset(cicerordsdedup, cicerordsdedup$coaccess <=0)
          #cicerordsdedup2.pos <- subset(cicerordsdedup, cicerordsdedup$coaccess >0)

          #hicloopsconns$in_zeronegcicero <- compare_connections(hicloopsconns, cicerordsdedup2.zeroneg, maxgap=1000)
          #hicloopsconns$in_poscicero <- compare_connections(hicloopsconns, cicerordsdedup2.pos, maxgap=1000)
          #hicoverlappingcicerodf <- as.data.frame(table(hicloopsconns$in_zeronegcicero, hicloopsconns$in_poscicero  ))
          #colnames(hicoverlappingcicerodf)[1:2] <- c("Cicero.zeroneg", "Cicero.pos")
          #write.table(hicoverlappingcicerodf,
          #            file=paste0("results/PCSC1/cicero/run2.readcount/", name, ".Freq.of.HiC.overlapped.by.cicero.v2.txt"),
          #            sep="\t", row.names=F, col.names=T, quote=F)

          ## HiC loops in Cicero
          cicerordsdedup$in_hic <- compare_connections(cicerordsdedup, hicloopsconns, maxgap=1000)

          # Restrict cicero loops to deduped ones
          #cicerordsdedup3 <- subset(cicerordsdedup2,cicerordsdedup2$pairID %in% dedupcon$pairID)

            summarisefunction(cicerordsdedup,paste0(name,".",hicfile,".v2" ))
    }

  } else if ( name=="LSCp") {

    hicloops <- read.table(paste0("data/HiC/LSC/All_Chr_Merged_Loops.bedpe"), header=T, sep="\t", stringsAsFactors=F)
    hicloops$Peak1 <- paste(paste0("chr",hicloops[,1]), hicloops[,2], hicloops[,3], sep="_")
    hicloops$Peak2 <- paste(paste0("chr",hicloops[,4]), hicloops[,5], hicloops[,6], sep="_")
    hicloopsconns <- hicloops[,c("Peak1", "Peak2")]

    ## Cicero loops found in HiC
    cicerordsdedup2.zeroneg <- subset(cicerordsdedup2, cicerordsdedup2$coaccess <=0)
    cicerordsdedup2.pos <- subset(cicerordsdedup2, cicerordsdedup2$coaccess >0)

    hicloopsconns$in_zeronegcicero <- compare_connections(hicloopsconns, cicerordsdedup2.zeroneg, maxgap=1000)
    hicloopsconns$in_poscicero <- compare_connections(hicloopsconns, cicerordsdedup2.pos, maxgap=1000)
    hicoverlappingcicerodf <- as.data.frame(table(hicloopsconns$in_zeronegcicero, hicloopsconns$in_poscicero  ))
    colnames(hicoverlappingcicerodf)[1:2] <- c("Cicero.zeroneg", "Cicero.pos")
    write.table(hicoverlappingcicerodf,
                file=paste0("results/PCSC1/cicero/run2.readcount/", name,
                ".Freq.of.HiC.overlapped.by.cicero.v2.txt"),
                sep="\t", row.names=F, col.names=T, quote=F)

    ## HiC loops in Cicero
    cicerordsdedup2$in_hic <- compare_connections(cicerordsdedup2, hicloopsconns, maxgap=1000)
    # Restrict cicero loops to deduped ones
    #cicerordsdedup3 <- subset(cicerordsdedup2,cicerordsdedup2$pairID %in% dedupcon$pairID)
    summarisefunction(cicerordsdedup2,name, ".v2")

    ## Shuffle regions
    cicerordsdedup2$coaccessflag2 <- as.character(cut(cicerordsdedup2$coaccess, seq(0,1, 0.25) ) )
    cicerordsdedup2$coaccessflag2 <-  ifelse( is.na(cicerordsdedup2$coaccessflag2)==TRUE, "zeronegative", cicerordsdedup2$coaccessflag2)
    p1 <- (as.numeric(str_split_fixed(cicerordsdedup2$Peak1,"_", 3)[,3]) + as.numeric(str_split_fixed(cicerordsdedup2$Peak1,"_", 3)[,2])) /2
    p2 <- (as.numeric(str_split_fixed(cicerordsdedup2$Peak2,"_", 3)[,3]) + as.numeric(str_split_fixed(cicerordsdedup2$Peak2,"_", 3)[,2])) /2

    cicerordsdedup2$distance <- (abs(p2-p1))/1000
    cicerordsdedup2$distancerange <-   as.character(cut(cicerordsdedup2$distance, seq(0, 300, 50) ) )
    cicerordsdedup2$distancerange <- ifelse(is.na(cicerordsdedup2$distancerange)==TRUE, "(300,above]", cicerordsdedup2$distancerange)
    cicerordsdedup2$megaflag <- paste( cicerordsdedup2$distancerange, cicerordsdedup2$coaccessflag2, sep=";")
    cicerordsdedup2$chromflag <- substr(cicerordsdedup2$Peak1,1,5)

    combined.shuf=NULL
    for ( f in unique(cicerordsdedup2$chromflag)){
      chk1 <- subset(cicerordsdedup2, cicerordsdedup2$chromflag==f)
      for (mega in unique(cicerordsdedup2$megaflag)){
        chk <- subset(chk1, chk1$megaflag==mega)
        print(paste0(mega, "_",f))
         for (k in 1){
            print(k)
            chk$end1shuf  <- chk[sample(nrow(chk)), "Peak1"]
            chk$end2shuf  <- chk[sample(nrow(chk)), "Peak2"]
            chk <- chk[,c("end1shuf","end2shuf", "in_hic")]
            colnames(chk)[1:2] <- c("Peak1","Peak2")
            chk$in_hic2 <- compare_connections(chk, hicloopsconns, maxgap=1000)
            combined.shuf <- rbind(combined.shuf, chk)
          }
      }
    }

    write.table(combined.shuf, file="results/PCSC1/cicero/run2.readcount/LSCp.shuffled.overlap.with.HiC.txt", quote=F, row.names=F, col.names=T, sep="\t")

  }
}

comparehicfunction("GBM")
comparehicfunction("LSCp")



## Plot results of overlap with Hic
files <- dir("results/PCSC1/cicero/run2.readcount/", "dcast.Percentage.txt", full.names=T)
for (f in 1:length(files)){

  a <- read.table(files[f], sep="\t", stringsAsFactors=F, header=T, check.names=F)
  name <- gsub(".Hic.overlap.1kb.dcast.Percentage.txt", "", basename(files[f]))

  df <- melt(a)
  df$variable2 <- factor(df$variable, c("(0,50]", "(50,100]", "(100,150]", "(150,200]", "(200,250]" ,"(250,300]" ,"(300,above]"))

  pdf(paste0("results/PCSC1/cicero/run2.readcount/" , name, ".pdf"))
  p <- ggplot(data=df, aes(x=variable2, y=value, group=grplist)) +
  geom_line(aes(color=grplist))+
  geom_point(aes(color=grplist)) +
  theme(text = element_text(size=12)) +
  labs(x = "Distance" , y="% of Cicero interactions overlapped by HiC") +
  annotate("text",  x=Inf, y = Inf, label = gsub(".v2","",gsub("GBM.","",name)), vjust=2, hjust=2)


  print(p)
  dev.off()

}
