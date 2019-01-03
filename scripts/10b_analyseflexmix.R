### R-3.4.1
### mordor
### Objective : Examine clustering results from flexmix
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


######################################
### load data
######################################
files=dir("results/PCSC1/Cluster/Flexmix/Binomial/", pattern=".Rdata")

## Determine best cluster no based on AIC and BIX

BIC=c()
AIC=c()
nos_clus=c()
for (f in 1:length(files)){
  
  print(f)
  load(paste0("results/PCSC1/Cluster/Flexmix/Binomial/", files[f]))
  bic.temp <- BIC(enhmat_fmm)
  aic.temp <- AIC(enhmat_fmm)

  BIC <- c(BIC,bic.temp)
  AIC <- c(AIC,aic.temp)
  nos_clus=c(nos_clus,max(enhmat_fmm@cluster))
  
  plotdf <- parameters(enhmat_fmm)
  pdf(paste0("results/PCSC1/Cluster/Flexmix/Binomial/Heatmap.K",max(enhmat_fmm@cluster),".",f,".pdf"))
  heatmap.2(t(plotdf), Rowv=T,Colv=F,
           scale="none", trace="none",col=hmcols); 
  dev.off()
  
  rm(enhmat_fmm)
}

### Only till K23 convergence was reached

ic.df <- data.frame(stringsAsFactors = F, 
                    cluster=c(nos_clus), 
                    BIC=c(BIC ),
                    AIC=AIC)

ic.df <- ic.df[order(ic.df$cluster),]
write.table(ic.df, file="results/PCSC1/Cluster/Flexmix/Binomial/IC.Criterion.txt", sep="\t", row.names=F, col.names=T,quote=F)

pdf("results/PCSC1/Cluster/Flexmix/Binomial/IC.Criterion.pdf")
#plot(ic.df$BIC, col="white")
#lines(ic.df$BIC, col="blue")
#lines(ic.df$AIC, col="green")
plot(((ic.df$AIC)), col="white")
lines(((ic.df$AIC)), col="green")
lines(((ic.df$BIC)), col="blue", lwd=2, lty=3)
dev.off()

