## Identify no. of E-P interactions mapping to common genes and different genes

library(data.table)
library(ggplot2)

gbmep <- fread("results/PCSC1/C3D/finalresults/GBM.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)
lscep <- fread("results/PCSC1/C3D/finalresults/LSCp.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)
pfaep <- fread("results/PCSC1/C3D/finalresults/PFA.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)

files <- dir("results/PCSC1/C3D/finalresults/", full.names=T, pattern=".interactions.genes.txt")
files <- files[-grep("GBM", files)]
files <- files[-grep("LSCp", files)]
files <- files[-grep("PFA", files)]

notshared <-  files[grep("notshared", files)]
common <- files[grep("common", files)]

namelist <- c()
gbmcommon <- c()
lsccommon <- c()
pfacommon <- c()

for (  f in 1:length(common)){

  namelist <- c(namelist,  gsub(".common.coaccess.regions.bed","",basename(common[f])))
  a <- fread(common[f], header=F, sep="\t", stringsAsFactors=F, data.table=F)
  gbmcommon <- c(gbmcommon, (nrow(subset(gbmep, gbmep$V4 %in% a$V1))))
  lsccommon <- c(lsccommon, (nrow(subset(lscep, lscep$V4 %in% a$V1))))
  pfacommon <- c(pfacommon, (nrow(subset(pfaep, pfaep$V4 %in% a$V1))))
}


namelist <- c()
gbmnotshared <- c()
lscnotshared <- c()
pfanotshared <- c()

for (  f in 1:length(notshared)){

  namelist <- c(namelist, gsub("notshared.coaccess","",basename(notshared[f])))
  a <- fread(notshared[f], header=F, sep="\t", stringsAsFactors=F, data.table=F)
  gbmnotshared <- c(gbmnotshared, (nrow(subset(gbmep, gbmep$V4 %in% subset(a$V1, a$V2=="GBM")))))
  lscnotshared <- c(lscnotshared, (nrow(subset(lscep, lscep$V4 %in% subset(a$V1, a$V2=="LSCp")))))
  pfanotshared <- c(pfanotshared, (nrow(subset(pfaep, pfaep$V4 %in% subset(a$V1, a$V2=="PFA")))))

}

commondf <- data.frame(coaccess=namelist,
                       stringsAsFactors=F,
                       epfreq=c(gbmcommon, lsccommon, pfacommon),
                      cancer=c(rep("GBM", length(gbmcommon)), rep("LSCp", length(lsccommon)), rep("PFA", length(pfacommon)) )
                      )
commondf$epfreqper100k <- commondf$epfreq/100000

## only edited till this poitn
pdf("results/PCSC1/C3D/finalresults/No.of.Common.EP.pairs.pdf")
ggplot(commondf, aes(x=coaccess, y=epfreqper100k, group=cancer)) +
geom_line(aes(color=cancer))+
geom_point(aes(color=cancer)) +
geom_text_repel(aes(label=commondf$label))+
labs(y="No.of E-P Pairs (x 100k)", x="Coaccessibility")+
scale_y_continuous(breaks=seq(0,15,1))
dev.off()


notshareddf <- data.frame(coaccess=namelist,
                       stringsAsFactors=F,
                       epfreq=c(gbmnotshared, lscnotshared, pfanotshared),
                      cancer=c(rep("GBM", length(gbmnotshared)), rep("LSCp", length(lscnotshared)), rep("PFA", length(pfanotshared)) )
                      )
notshareddf$epfreqper100k <- notshareddf$epfreq/100000
notshareddf$label <- ifelse(notshareddf$coaccess < 0.8, "", notshareddf$epfreq)

pdf("results/PCSC1/C3D/finalresults/summary.plots/No.of.notshared.EP.pairs.pdf")
ggplot(notshareddf, aes(x=coaccess, y=epfreqper100k, group=cancer)) +
geom_line(aes(color=cancer))+
geom_point(aes(color=cancer)) +
geom_text_repel(aes(label=notshareddf$label))+
labs(y="No.of E-P Pairs (x 100k)", x="Coaccessibility")+
scale_y_continuous(breaks=seq(0,15,1))
dev.off()



## Proportion


namelist <- c()
gbmcommon <- c()
lsccommon <- c()
pfacommon <- c()

for (  f in 1:length(common)){

  namelist <- c(namelist, gsub("common.coaccessgt","", gsub(".E-P.genes.txt","",basename(common[f]))))
  a <- fread(common[f], header=F, sep="\t", stringsAsFactors=F, data.table=F)
  gbmcommon <- c(gbmcommon, (nrow(subset(gbmep, gbmep$V4 %in% a$V1)))/nrow(gbmep) )
  lsccommon <- c(lsccommon, (nrow(subset(lscep, lscep$V4 %in% a$V1)))/nrow(lscep))
  pfacommon <- c(pfacommon, (nrow(subset(pfaep, pfaep$V4 %in% a$V1)))/nrow(pfaep))

}


namelist <- c()
gbmnotshared <- c()
lscnotshared <- c()
pfanotshared <- c()

for (  f in 1:length(notshared)){

  namelist <- c(namelist, gsub("notshared.coaccessgt","", gsub(".E-P.genes.txt","",basename(notshared[f]))))
  a <- fread(notshared[f], header=F, sep="\t", stringsAsFactors=F, data.table=F)
  gbmnotshared <- c(gbmnotshared, (nrow(subset(gbmep, gbmep$V4 %in% subset(a$V1, a$V2=="GBM"))))/nrow(gbmep))
  lscnotshared <- c(lscnotshared, (nrow(subset(lscep, lscep$V4 %in% subset(a$V1, a$V2=="LSCp"))))/nrow(lscep))
  pfanotshared <- c(pfanotshared, (nrow(subset(pfaep, pfaep$V4 %in% subset(a$V1, a$V2=="PFA"))))/nrow(pfaep))

}

df <- data.frame(coaccess=namelist,
                       stringsAsFactors=F,
                       freq=c(c((gbmcommon)*100, (lsccommon)*100, (pfacommon)*100), c((gbmnotshared)*100, (lscnotshared)*100, (pfanotshared)*100)),
                       grp=c(rep("common", 21 ), rep("notshared", 21)),
                      cancer=rep(c(rep("GBM", 7 ), rep("LSCp", 7), rep("PFA", 7)),2 )
                      )
 df$freqrounded <- round(df$freq, 2)

pdf("results/PCSC1/C3D/finalresults/summary.plots/Prop.of.GBM.common.notshared.EP.pairs.pdf")
ggplot() +
geom_bar(aes(y = freqrounded, x = coaccess, fill = grp), data = subset(df, df$cancer=="GBM"),stat="identity")+
scale_y_continuous(breaks=seq(0,100,20))+
labs(y="% of E-P pairs\n", x="\nCoaccessibility cutoff")+
#geom_text(data=subset(df, df$cancer=="GBM"), aes(x = coaccess, y = freqrounded,  label = paste0(freqrounded,"%")), size=4) +
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),text = element_text(size=16), legend.spacing.x = unit(1.0, 'cm')) +
annotate("text",  x=Inf, y = Inf, label = "GBM", vjust=2, hjust=2)+
scale_fill_manual(values=c("#5F9EA0", "#E1B378"))
dev.off()


pdf("results/PCSC1/C3D/finalresults/summary.plots/Prop.of.PFA.common.notshared.EP.pairs.pdf")
ggplot() +
geom_bar(aes(y = freqrounded, x = coaccess, fill = grp), data = subset(df, df$cancer=="PFA"),stat="identity")+
scale_y_continuous(breaks=seq(0,100,20))+
labs(y="% of E-P pairs\n", x="\nCoaccessibility cutoff")+
#geom_text(data=subset(df, df$cancer=="PFA"), aes(x = coaccess, y = freqrounded,  label = paste0(freqrounded,"%")), size=4) +
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),text = element_text(size=16), legend.spacing.x = unit(1.0, 'cm')) +
scale_fill_manual(values=c("#5F9EA0", "#E1B378"))+
annotate("text",  x=Inf, y = Inf, label = "PFA", vjust=2, hjust=2)
dev.off()


pdf("results/PCSC1/C3D/finalresults/summary.plots/Prop.of.LSCp.common.notshared.EP.pairs.pdf")
ggplot() +
geom_bar(aes(y = freqrounded, x = coaccess, fill = grp), data = subset(df, df$cancer=="LSCp"),stat="identity")+
scale_y_continuous(breaks=seq(0,100,20))+
labs(y="% of E-P pairs\n", x="\nCoaccessibility cutoff")+
#geom_text_repel(data=subset(df, df$cancer=="LSCp"), aes(x = coaccess, y = freqrounded,  label = paste0(freqrounded,"%")), size=4) +
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),text = element_text(size=16),legend.spacing.x = unit(1.0, 'cm')) +
scale_fill_manual(values=c("#5F9EA0", "#E1B378"))+
annotate("text",  x=Inf, y = Inf, label = "LSCp", vjust=2, hjust=2)
dev.off()



df <- subset(df, df$coaccess==0)
df$coaccess <- NULL
df <- df[,c("cancer","freq","grp","freqrounded")]
df$grp <- ifelse(df$grp=="common", "Common genes","Unique genes")
pdf("results/PCSC1/C3D/finalresults/summary.plots/Prop.of.common.notshared.EP.pairs.pdf")
ggplot() +
geom_bar(aes(y = freqrounded, x = cancer, fill = grp), data = df,stat="identity")+
scale_y_continuous(breaks=seq(0,100,20))+
labs(y="% of distal-Promoter pairs\n", x="\nCoaccessibility cutoff")+
#geom_text(data=df, aes(x = cancer, y = freqrounded,  label = paste0(freqrounded,"%")), size=4) +
theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(),
  text = element_text(size=20), legend.spacing.x = unit(1.0, 'cm')) +
#annotate("text",  x=Inf, y = Inf, label = "CancerType", vjust=2, hjust=2)+
scale_fill_manual(values=c("#5F9EA0", "#E1B378"))
dev.off()
