### R-3.5.0
### mordor

# =====================================================================================================================================
# Objective :  Proportion of cicero loops lying in TADs
              # Inter vs Intra TADs
              # Ratio of cicero loops in same TAD compared to unlinked pairs of same distance
# =====================================================================================================================================

#-------------------------------------------------------------------------------------------------
# load dependencies
#-------------------------------------------------------------------------------------------------
library(data.table)
library(stringr)
#-------------------------------------------------------------------------------------------------
  #cut -f1 results/PCSC1/cicero/run2.readcount/LSCp.deduped.connsfull.pairs.bed | cut -d"_" -f1,2,3 --output-delimiter=$'\t' | intersectBed -a stdin -b data/HiC/LSC/10kbTADs.bed -u -f 1 > end1.overlap.bed
#cut -f2 results/PCSC1/cicero/run2.readcount/LSCp.deduped.connsfull.pairs.bed | cut -d"_" -f1,2,3 --output-delimiter=$'\t' | intersectBed -a stdin -b data/HiC/LSC/10kbTADs.bed -u -f 1 > end2.overlap.bed

loop <- fread("results/PCSC1/cicero/run2.readcount/LSCp.deduped.connsfull.pairs.bed",
              header=F, sep="\t", stringsAsFactors=F, data.table=F)

end1 <- fread("end1.overlap.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)
end2 <- fread("end2.overlap.bed", header=F, sep="\t", stringsAsFactors=F, data.table=F)
end1$id <- paste(end1[,1], end1[,2], end1[,3], sep="_")
end2$id <- paste(end2[,1], end2[,2], end2[,3], sep="_")

loop$end1 <- ifelse(loop$V1 %in% end1$id,1,0)
loop$end2 <- ifelse(loop$V2 %in% end2$id,1,0)
loop$pair <- paste(loop$V1, loop$V2, sep=":")
loop$flag <- paste0(loop$end1, loop$end2)
loop$flag2 <- ifelse(loop$flag=="11","intad","outsidetad")

tmp <- readRDS("results/PCSC1/cicero/run2.readcount/LSCp.connsfull.rds")
tmp$pair <- paste(tmp$Peak1, tmp$Peak2, sep=":")
tmp$coaccess2 <-  ifelse( is.na(tmp$coaccess)==TRUE, -1, tmp$coaccess)

## Add coaccess info
loop2 <- merge(loop, tmp[,c("coaccess2","pair")], by.x="pair", by.y="pair",all.x=T)
loop2$coaccessflag <- as.character(cut(loop2$coaccess2, seq(0,1, 0.2) ) )
loop2$coaccessflag <-  ifelse( is.na(loop2$coaccessflag)==TRUE, "zeronegative", loop2$coaccessflag)

## Linked vs unlinked pairs at same distance
p1 <- (as.numeric(str_split_fixed(loop2$V1,"_", 3)[,3]) + as.numeric(str_split_fixed(loop2$V1,"_", 3)[,2])) /2
p2 <- (as.numeric(str_split_fixed(loop2$V2,"_", 3)[,3]) + as.numeric(str_split_fixed(loop2$V2,"_", 3)[,2])) /2

loop2$distance <- (abs(p2-p1))/1000
loop2$distancerange <-   as.character(cut(loop2$distance, seq(0, 300, 50) ) )
loop2$distancerange <- ifelse(is.na(loop2$distancerange)==TRUE, "(300,above]", loop2$distancerange)

## Remove those not in TAD at all
loop2 <- subset(loop2, loop2$flag !="00")


freqdf <- as.data.frame(table(loop2$flag2, loop2$distancerange, loop2$coaccessflag))
freqdf.wide <- dcast(freqdf, Var1+Var2~Var3)
opname="LSCp"
write.table(freqdf, file=paste0("results/PCSC1/cicero/run2.readcount/",opname ,".TAD.overlap.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(  freqdf.wide, file=paste0("results/PCSC1/cicero/run2.readcount/",opname, ".HTADic.overlap.dcast.txt"), sep="\t", quote=F, row.names=F, col.names=T)

pctlist <- c()
grplist <- c()
grp2list <- c()
for ( f in unique(  freqdf.wide$Var2)){
  tmp1 <- subset(freqdf.wide, freqdf.wide$Var2==f)
  tmp1 <- t(tmp1[,-c(1,2)])
  colnames(tmp1) <- c("intad","outsidetad")
  tmp1 <- as.data.frame(tmp1)
  tmp1$pct <- (tmp1$intad/(tmp1$intad + tmp1$outsidetad))*100
  tmp1$grp <- f
  colnames(tmp1) <- paste0(colnames(tmp1), ".",f)
  pctlist <- c(pctlist, tmp1$pct)
  grplist <- c(grplist, tmp1$grp)
  grp2list <- c(grp2list, rownames(tmp1))
}

summ <- data.frame(stringsAsFactors=F, grplist, grp2list,pctlist)
summ[is.na(summ)] <- 0
summ2 <- dcast(summ, grplist~grp2list)
write.table( summ2, file=paste0("results/PCSC1/cicero/run2.readcount/", "LSCp", ".Hic.TADoverlap.dcast.Percentage.txt"), sep="\t", quote=F, row.names=F, col.names=T)


pctlist <- c()
grplist <- c()
grp2list <- c()
for ( f in unique(  freqdf.wide$Var2)){
  tmp1 <- subset(freqdf.wide, freqdf.wide$Var2==f)
  tmp1 <- t(tmp1[,-c(1,2)])
  colnames(tmp1) <- c("intad","outsidetad")
  tmp1 <- as.data.frame(tmp1)
  tmp1$pct <- (tmp1$intad)/(tmp1$outsidetad)
  tmp1$grp <- f
  colnames(tmp1) <- paste0(colnames(tmp1), ".",f)
  pctlist <- c(pctlist, tmp1$pct)
  grplist <- c(grplist, tmp1$grp)
  grp2list <- c(grp2list, rownames(tmp1))
}

summ <- data.frame(stringsAsFactors=F, grplist, grp2list,pctlist)
summ[is.na(summ)] <- 0
summ2 <- dcast(summ, grplist~grp2list)
write.table( summ2, file=paste0("results/PCSC1/cicero/run2.readcount/", "LSCp", ".Hic.TADoverlap.dcast.OR.txt"), sep="\t", quote=F, row.names=F, col.names=T)

## Plot

summ2 <- summ2[c(1,7,2:6),c(7,1:6)]
summ3 <- melt(summ2)
summ3$variable <- as.character(summ3$variable)
summ3$grplist <- factor(summ3$grplist, levels=c("(0,50]","(50,100]","(100,150]","(150,200]","(200,250]","(250,300]","(300,above]") )
library(ggplot2)
pdf(paste0("results/PCSC1/cicero/run2.readcount/", "LSCp", ".Hic.TADoverlap.dcast.OR.pdf"))
ggplot(summ3, aes(x=grplist, y=value, group=variable)) +
geom_line(aes(color=variable))+
geom_point(aes(color=variable)) +
labs(y="Intra/inter TAD Cicero LSCp loops", x="Distance range")+
scale_y_continuous(breaks=seq(0,12.5,1.5))
#scale_color_manual(values=c("#CC6666", "#9999CC"))
dev.off()


## Define random regions
## For each distance bin randomly assigning two
## distance-matched peaks as linked and retaining the same total number of links for each coaccessibility

## split file by distance
## shuffle order
## calculate overlap
## Repeat 10K times
loop2$chromflag <- substr(loop2$V1,1,5)
loop2$coaccessflag2 <- as.character(cut(loop2$coaccess2, seq(0,1, 0.25) ) )
loop2$coaccessflag2 <-  ifelse( is.na(loop2$coaccessflag2)==TRUE, "zeronegative", loop2$coaccessflag2)
saveRDS(loop2, file="results/PCSC1/cicero/run2.readcount/LSCp.TAD.Overlap.dedup.connsfull.rds")

write.table(loop2, file="results/PCSC1/cicero/run2.readcount/LSCp.TAD.Overlap.dedup.connsfull.bed",
            sep="\t", row.names=F, col.names=T, quote=F)

loop2$megaflag <- paste( loop2$distancerange, loop2$coaccessflag2, sep=";")

loop3 <- subset(loop2, loop2$coaccessflag!="zeronegative")

loop3$megaflag <- paste( loop3$distancerange, loop3$coaccessflag2, sep=";")

shuflist <- c()
realist <- c()
namelist <- c()

for ( f in unique(loop3$chromflag)){

  chk1 <- subset(loop3, loop3$chromflag==f)

  for (mega in unique(loop3$megaflag)){

     chk <- subset(chk1, chk1$megaflag==mega)
     real <- table(chk$flag2)[1]/table(chk$flag2)[2]
     realist <- c(realist, real)
    namelist <- c(namelist, paste0(mega, "_",f))
    print(paste0(mega, "_",f))

     for (k in 1){
        print(k)
        chk$end1shuf  <- chk[sample(nrow(chk)), "V1"]
        chk$end2shuf  <- chk[sample(nrow(chk)), "V2"]
        chk$end1shufflag <- ifelse(chk$end1shuf %in% end1$id,1,0)
        chk$end2shufflag <- ifelse(chk$end2shuf %in% end2$id,1,0)
        chk$shufflag <- paste0(chk$end1shufflag, chk$end2shufflag)
        chk$shufflag2 <- ifelse(chk$shufflag=="11","intad","outsidetad")
        #chk <- subset(chk, chk$shufflag2 !="00")
        shuf <- table(chk$shufflag2)[1]/table(chk$shufflag2)[2]
        shuflist <- c(shuflist, shuf)
      }
  }
}
