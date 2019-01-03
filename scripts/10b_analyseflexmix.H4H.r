### R-3.4.1
### h4h
### Objective : Examine clustering results from flexmix on h4h
### NOT RUN -- Needs to be optimised


######################################
### load dependencies
######################################
require(flexmix)          # For model bases clustering
require(RColorBrewer)     # For color palettes in plotting
require(ggplot2)          # For plotting
require(FactoMineR)       # For performing MCA
require(vegan)            # For calculating jaccard distance
require(gplots)           # For plotting heatmaps
require(Heatplus)         # For plotting heatmaps with annotations
require(reshape)          # For melting dataframe
library(RColorBrewer)

hmcols = colorpanel(100, "steelblue", "white", "tomato")
bluered299 <- colorRampPalette(c("blue","royalblue","aliceblue","brown1","red"))(n=299)
color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)


######################################
### load data
######################################
files=dir("results/PCSC1/Cluster/Flexmix/Binomial.h4h/", pattern=".Rdata")

## Determine best cluster no based on AIC and BIX

BIC=c()
AIC=c()
EIC=c()
nos_clus=c()

for (f in 1:length(files)){
  rm(enhmat_fmm)
  
  print(f)
  load(paste0("results/PCSC1/Cluster/Flexmix/Binomial.h4h/", files[f]))
  
  bic.temp <- BIC(enhmat_fmm)
  aic.temp <- AIC(enhmat_fmm)
  eic.temp <- EIC(enhmat_fmm)
  
  BIC <- c(BIC,bic.temp)
  AIC <- c(AIC,aic.temp)
  EIC <- c(EIC,eic.temp)
  nos_clus=c(nos_clus,max(enhmat_fmm@cluster))
  name=c(name, files[f])
  plotdf <- parameters(enhmat_fmm)
  pdf(paste0("results/PCSC1/Cluster/Flexmix/Binomial.h4h/Heatmap.K",max(enhmat_fmm@cluster),".",f,".",files[f],".pdf"))
  heatmap.2(t(plotdf),col=hmcols, scale="none",trace="none",Colv=NULL,cexRow=0.3,
            hclustfun=function(x) hclust(x,method="ward.D2"), 
            distfun=function(x) as.dist((1 - cor(  t(x), method="spearman"  ))))
  dev.off()
}

ic.df <- data.frame(stringsAsFactors = F, 
                    cluster=c(nos_clus),
                    name=name,
                    BIC=c(BIC ),
                    AIC=AIC)
ic.df$actualK <- as.numeric(unlist( str_split_fixed(gsub("Enhancer.Binarymat.Fmm.K", "", gsub(".Rdata","", ic.df$name)), "to",2))[,1])
ic.df <- ic.df[order(ic.df$actualK),]
ic.df <- unique(ic.df)
write.table(ic.df, file="results/PCSC1/Cluster/Flexmix/Binomial.h4h/IC.Criterion.txt", sep="\t", row.names=F, col.names=T,quote=F)

ic.df.k35 <- subset(ic.df, ic.df$actualK <= 35)
pdf("results/PCSC1/Cluster/Flexmix/Binomial.h4h/IC.Criterion.pdf", useDingbats = F)
plot(ic.df.k35$BIC, col="white")
lines(ic.df.k35$BIC, col="blue")
lines(ic.df.k35$AIC, col="green")
dev.off()

## Results can only be trusted till K=33 since the no. of clusters identified and the no. of K assigned is not the same


load(paste0("results/PCSC1/Cluster/Flexmix/Binomial.h4h/Enhancer.Binarymat.Fmm.K30to30.Rdata"))
plotdf <- parameters(enhmat_fmm)
## Ordered by me
# 
# ## Assign groups
# common <- c(30,27,6,11,12)
# shared <- c()
# lsc.sp <- c(22,20,25,9,24)
# gbm.sp <- c()
# pfa.sp <- c()
# 
# pdf(paste0("results/PCSC1/Cluster/Flexmix/Binomial.h4h/Heatmap.K30.v3.pdf"))
# heatmap.2(t(plotdf),col=hmcols, scale="none",trace="none",Colv=NULL,cexRow=0.3,
#           hclustfun=function(x) hclust(x,method="ward.D2"), 
#           distfun=function(x) as.dist((1 - cor(  t(x), method="spearman"  ))))
# dev.off()

