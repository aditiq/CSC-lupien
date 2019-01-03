#===========================================================================================================================
# Objective : Plot distribution of co-accessibiity scores from cicero
# mordor
# R-3.5.0
#===========================================================================================================================

#===================================================================
# load dependencies
#===================================================================

library(data.table)
library(plyr)
library(ggplot2)
library(ggrepel)

#===================================================================
# Run summary plots
#===================================================================

plotfunc=function(ds, opname, ...) {

  df <- fread(ds, sep="\t", stringsAsFactors=F, header=T, data.table=F)
  df$anno1 <- ifelse( is.na(df$geneid.peak1 )==TRUE,"E","P")
  df$anno2 <- ifelse( is.na(df$geneid.peak2 )==TRUE,"E","P")
  df <- df[,c("Peak1", "Peak2","coaccess","anno1","anno2")]

  ep <- subset(df, paste0(df$anno1, df$anno2) %in% c("EP","PE"))
  switch <- subset(ep, ep$anno1=="E")
  switch2 <- switch[,c(2,1,3,5,4)]
  colnames(switch2) <- colnames(ep)
  ep2 <- rbind(subset(ep, ep$anno1=="P"), switch2)
  ep2 <- ep2[!duplicated(ep2),]

  promlist <- c()
  enhlist <- c()
  coaccesslist <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  for ( f in coaccesslist)
  {
    promlist <- c(promlist, print(median(count(subset(ep2,ep2$coaccess > f), "Peak1")$freq)))
    enhlist <- c(enhlist,print(median(count(subset(ep2,ep2$coaccess > f), "Peak2")$freq)))
  }

  ##---------------------------------------------------------------------------------------------
  ## plot no. of promoters/distal site and no. of distal sites/promoters by co-accessibiity score
  ##---------------------------------------------------------------------------------------------

  plotdf1 <- data.frame(stringsAsFactors=F, coaccess=coaccesslist, freq=c(promlist, enhlist),
                        name=c(rep("prom", length(coaccesslist)), rep("enh", length(coaccesslist))))

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".No.of.coaccess.sites.pdf"))
  p <- ggplot(plotdf1, aes(x=coaccess, y=freq, group=name)) +
  geom_line(aes(color=name))+
  geom_point(aes(color=name)) +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  print(p)
  dev.off()

  ##---------------------------------------------------------------------------------------------
  ## plot distribution of co-accessibiity scores
  ##---------------------------------------------------------------------------------------------

  df2 <- df[!duplicated(df),]

  df2$coaccess2 <- ifelse(df2$coaccess >0.5, 0.5, df2$coaccess)

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Distribution.coaccess.pdf"))
  histp <- ggplot(df2, aes(x=(coaccess))) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins= 1000)+
  geom_density(alpha=.2, fill="#FF6666")
  print(histp)
  dev.off()

  ##---------------------------------------------------------------------------------------------
  ## plot distribution of co-accessibiity scores as a function of pairs identified
  ##---------------------------------------------------------------------------------------------

  plotdf2 <- data.frame(Var1=c(), Var2=c(), Freq=c(), coaccess=c())

  for ( f in coaccesslist) {
    print(f)
    dftemp <- subset(df2, df2$coaccess > f)
    plotdf2  <- rbind(plotdf2 , cbind(as.data.frame(table(dftemp$anno1, dftemp$anno2)), coaccess=f))
  }

  plotdf2$grp <- paste(plotdf2$Var1, plotdf2$Var2, sep="-")
  plotdf2$grp <- gsub("P-E", "E-P" ,  plotdf2$grp  )

  plotdf2 <- aggregate(plotdf2$Freq, by=list(grp=plotdf2$grp, coaccess=plotdf2$coaccess), FUN=sum)
  #plotdf2$x <- ifelse(plotdf2$x > 1000000, 1000000, plotdf2$x )
  write.table(plotdf2,
              file=paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Freq.of.pairs.by.coaccess.bed"),
              sep="\t", row.names=F, col.names=T, quote=F)

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Freq.of.pairs.by.coaccess.pdf"))
  p <- ggplot(plotdf2, aes(x=coaccess, y=x, group=grp, label=x)) +
  geom_line(aes(color=grp))+
  geom_point(aes(color=grp)) +
    geom_text_repel() +
  scale_color_manual(values=c("#999999", "#E69F00", "#00e69f")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  print(p)
  dev.off()

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Freq.of.E-P.pairs.by.coaccess.pdf"))
  p3 <- ggplot(subset(plotdf2,plotdf2$grp == "E-P"), aes(x=coaccess, y=x, group=grp, label=x)) +
  geom_line(aes(color=grp))+
  geom_point(aes(color=grp)) +
    geom_text_repel() +
  scale_color_manual(values=c("#999999", "#E69F00", "#00e69f")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  print(p3)
  dev.off()

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Freq.of.E-E.pairs.by.coaccess.pdf"))
  p1 <- ggplot(subset(plotdf2,plotdf2$grp == "E-E"), aes(x=coaccess, y=x, group=grp, label=x)) +
  geom_line(aes(color=grp))+
  geom_point(aes(color=grp)) +
    geom_text_repel() +
  scale_color_manual(values=c("#999999", "#E69F00", "#00e69f")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  print(p1)
  dev.off()

  pdf(paste0("results/PCSC1/cicero/run2.readcount/summary.plots/", opname, ".Freq.of.P-P.pairs.by.coaccess.pdf"))
  p2 <- ggplot(subset(plotdf2,plotdf2$grp == "P-P"), aes(x=coaccess, y=x, group=grp, label=x)) +
  geom_line(aes(color=grp))+
  geom_point(aes(color=grp)) +
  geom_text_repel() +
  scale_color_manual(values=c("#999999", "#E69F00", "#00e69f")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  print(p2)
  dev.off()

}


plotfunc("results/PCSC1/cicero/run2.readcount/GBM.Deduplicated.coaccess.Annotated.txt", "GBM")
plotfunc("results/PCSC1/cicero/run2.readcount/PFA.Deduplicated.coaccess.Annotated.txt", "PFA")
plotfunc("results/PCSC1/cicero/run2.readcount/LSCp.Deduplicated.coaccess.Annotated.txt", "LSCp")


## Combine plots


files <- dir("results/PCSC1/cicero/run2.readcount/summary.plots/", pattern=".Freq.of.pairs.by.coaccess.bed", full.names=T)

epdf <- NULL

for ( f in 1:length(files)) {

    a <- read.table(files[f], header=T, sep="\t", stringsAsFactors=F)
    a <- subset(a,a$grp=="E-P")
    a$cancertype <- gsub(".Freq.of.pairs.by.coaccess.bed","",basename(files[f]))
    epdf <- rbind(epdf, a)
}

colnames(epdf)[3] <- "Freq"
epdf$grp <- NULL
epdf$Freq <- epdf$Freq/100000
epdf$Freqlabel <- ifelse(epdf$coaccess < 0.8, "", epdf$Freq)

pdf("results/PCSC1/cicero/run2.readcount/summary.plots/CombinedFreq.of.pairs.by.coaccess.pdf")
ggplot(data=epdf, aes(x=coaccess, y=Freq, group=cancertype)) +
  geom_line(aes(color=cancertype))+
  geom_point(aes(color=cancertype)) +  
  geom_text_repel(aes(label = epdf$Freqlabel)) +
  theme(text = element_text(size=12)) + 
  labs(x = "Coaccessibility Score" , y="No. of E-P interactions (x 10K)") + 
   scale_colour_manual(values=c(GBM="#1e71aa",PFA="#ef9c1f",LSCp="#1b9978",Cadenza="#CC0033",Cad_A="#FF6600",Cad_B="#FF9933"))+
  scale_y_continuous(breaks=seq(0,18,2)) +
  scale_x_continuous(breaks=seq(0,1,0.1))
dev.off()


