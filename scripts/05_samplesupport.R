### R-3.4.1
### mordor
### Objective : Sample Support analysis


######################################
### load dependencies
######################################
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(UpSetR)

#####################################################################
#***********************VENN DIAGRAM************************
#####################################################################

samples<-list.files(path = "data/ConsensusSet/PCSC1/",recursive=TRUE,pattern="Consensus.Catalogue.narrowPeak")
samples <- samples[c(1,2,4)]

#create peaks list
peaks.l<-list()
number_peaks<-vector()
for (i in 1:length(samples)) {
  peaks.l[[i]]<-read.delim(paste("data/ConsensusSet/PCSC1/",samples[i],sep=""),sep="\t",header=F)
  peaks.l[[i]]<-peaks.l[[i]][,1:3]
  number_peaks<-c(number_peaks,nrow(peaks.l[[i]]))
}
names(peaks.l)<-samples

#transform the  bed files into GRanges objects
data2.gr<-list()
for (i in 1:length(samples)) {
  data2.gr[[i]]=GRanges(seqnames=peaks.l[[i]][,1],ranges=IRanges(start=peaks.l[[i]][,2],end=peaks.l[[i]][,3]))
}
names(data2.gr)<- samples

names(data2.gr[[1]]) <- paste0("GBM",seq(1,length(data2.gr[[1]]),1))
names(data2.gr[[2]]) <- paste0("LSCp",seq(1,length(data2.gr[[2]]),1))
names(data2.gr[[3]]) <- paste0("PFA",seq(1,length(data2.gr[[3]]),1))
ol <- findOverlapsOfPeaks(data2.gr[[1]], data2.gr[[2]], data2.gr[[3]],connectedPeaks="merge")

pdf("results/PCSC1/Venndiagram.GBM.PFA.LSCp.pdf") ;
makeVennDiagram(ol,NameOfPeaks=c("GBM","LSCp","PFA") , totalTest=1e+8, by="region",connectedPeaks = "merge")
dev.off()


#####################################################################
#***********************SAMPLE SUPPORT************************
#####################################################################
## list of samples
files <- dir("data/ConsensusSet/PCSC1/", pattern="Consensus.Catalogue.narrowPeak", full.names = T)
files <- files[c(1,2,4)]

runfunc=function(ds,...){
  peaks <- fread(ds, header=F, sep="\t", stringsAsFactors = F, data.table=F)
  dat <- data.frame(table(peaks$V4), stringsAsFactors = F)
  dat$Var1 <- as.numeric(as.character(dat$Var1))
  dat$cumsumfreq <-  rev(cumsum(rev(dat$Freq)))
  dat$sample.pct <- (dat$Var1 *100)/ max(dat$Var1)
  dat$peak.pct <- (dat$cumsumfreq *100)/  sum(dat$Freq)
  return(dat)
}


gbm <- runfunc(files[1])
lscp <- runfunc(files[2])
pfa <- runfunc(files[3])

## Plot all of them together
dat2 <- data.frame(stringsAsFactors = F, 
                   name=c(rep("LSCp", nrow(lscp)),rep("GBM", nrow(gbm)),rep("PFA", nrow(pfa)) ),
                   sample.pct = c(lscp$sample.pct, gbm$sample.pct, pfa$sample.pct),
                   peak.pct=c( lscp$peak.pct, gbm$peak.pct, pfa$peak.pct)
)

pdf("results/PCSC1/SampleSupport.pdf") ; 
ggplot(data=dat2, aes(x=peak.pct, y=sample.pct, color=name)) +  
  geom_line() +
  labs(x = "% of peaks covered", y= "% of samples") +
  geom_vline(xintercept = 50, linetype=3) +
  geom_hline(yintercept = 50, linetype=3); 
dev.off()

#####################################################################
#***********************UPSET PLOT************************
#####################################################################

binarymat <- fread("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Binarymat.txt", header=T, sep="\t",data.table=F, stringsAsFactors=F, check.names=F)
binarymat[,1] <- paste(binarymat[,1], binarymat[,2], binarymat[,3],sep="_")
binarymat$start <- NULL
binarymat$end <- NULL
colnames(binarymat) <- gsub("_peaks","", colnames(binarymat))

alignmeta <- fread("data/Alignment.Stats.txt", header=T, sep="\t", stringsAsFactors = F)
binarymat$LSC <- ifelse(rowSums(binarymat[,subset(alignmeta$Name, alignmeta$Cancer=="LSC" & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)])>0,1,0)
binarymat$GBM <- ifelse(rowSums(binarymat[,subset(alignmeta$Name, alignmeta$Cancer=="GBM" & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)])>0,1,0)
binarymat$PFA <- ifelse(rowSums(binarymat[,subset(alignmeta$Name, alignmeta$Cancer=="PFA" & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)])>0,1,0)

pdf("results/PCSC1/UpsetPlot.GBM.PFA.LSCp.pdf", useDingbats = F) ;
upset(binarymat[,c("seqnames","GBM","PFA","LSC")], order.by = "freq")
dev.off()



#####################################################################
#***********************PLOT SATURATION CURVE************************
#####################################################################
## list of samples
samples<-dir("data/peaks/", recursive=TRUE,pattern="_peaks.narrowPeak", full.names=T)
names <- gsub("_peaks.narrowPeak","",basename(samples))
names(samples) <- names


                    
#create peaks list
peaklistfunc=function(samplelist,...){
  peaks.l<-list()
  samplelist1=samplelist
  #number_peaks<-vector()
  for (i in 1:length(samplelist)) {
    peaks.l[[i]]<-fread(samplelist1[i],sep="\t",header=F, data.table=F, stringsAsFactors = F)
    peaks.l[[i]]<-peaks.l[[i]][,1:3]
    #number_peaks<-c(number_peaks,nrow(peaks.l[[i]]))
  }
  #names(peaks.l) <- namelist[f:j]
  #min_peaks<-min(number_peaks)
  return(peaks.l)
  #names(peaks.l) <- names(samplelist)
}

lsc.samples <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("LSC") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)]
gbm.samples <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("GBM") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)]
pfa.samples <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("PFA") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)]

lsc.peaks <- peaklistfunc(samplelist=lsc.samples)
gbm.peaks <- peaklistfunc(samplelist=gbm.samples)
pfa.peaks <- peaklistfunc(samplelist=pfa.samples)

names(lsc.peaks) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("LSC") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)
names(gbm.peaks) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("GBM") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)
names(pfa.peaks) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("PFA") & alignmeta$Group=="positive" & alignmeta$Inclusion.20kPeaks.15Mreads==1)


lsc.samples.all <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("LSC") & alignmeta$Group=="positive" )]
gbm.samples.all <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("GBM") & alignmeta$Group=="positive" )]
pfa.samples.all <- samples[subset(alignmeta$Name, alignmeta$Cancer %in% c("PFA") & alignmeta$Group=="positive")]

lsc.peaks.all <- peaklistfunc(samplelist=lsc.samples.all)
gbm.peaks.all <- peaklistfunc(samplelist=gbm.samples.all)
pfa.peaks.all <- peaklistfunc(samplelist=pfa.samples.all)

names(lsc.peaks.all) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("LSC") & alignmeta$Group=="positive")
names(gbm.peaks.all) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("GBM") & alignmeta$Group=="positive" )
names(pfa.peaks.all) <- subset(alignmeta$Name, alignmeta$Cancer %in% c("PFA") & alignmeta$Group=="positive" )


satfunc=function(n=100, maxsamps=100,samplelist,peaks.l,opname){

    numsam<-length(samplelist) #number of samples
    data.l<-list() #create the results list
    
    print("starting")
    for (j in 1:n) {  
    #Take a random peak list from the set of numsam lists, regardless of number of peaks
    sampled_peaks.l<-list()
    random_order<-sample(1:numsam, numsam, replace=F)
    for (i in 1:numsam) {
      sampled_peaks.l[[i]]<-peaks.l[[random_order[i]]]
    }
    
    print("Step2")
    #transform the sampled bed files into GRanges objects
    data.gr<-list()
    for (i in 1:numsam) {
      data.gr[[i]]=GRanges(seqnames=sampled_peaks.l[[i]][,1],
                           ranges=IRanges(start=sampled_peaks.l[[i]][,2],
                                          end=sampled_peaks.l[[i]][,3]))
    }
    
    print("Create Sat Curve")
    #create saturation curve data matrix
    data.l[[j]]<-matrix(rep(0,numsam),ncol=1,nrow=numsam)
    
    union_peaks<-data.gr[[1]]
    data.l[[j]][1,1]<-length(union_peaks)
    
    
    for (i in 2:numsam) {
      union_peaks2<-c(union_peaks,data.gr[[i]])
      union_peaks<-reduce(union_peaks2)
      data.l[[j]][i,1]<-length(union_peaks)
      
    }
    }
    
    #Calculate summary statistics
    data.m<-data.frame(as.numeric(seq(from=1,to=numsam)))
    for (j in 1:n) {data.m<-cbind(data.m,data.l[[j]][,1])} 
    names(data.m)<-c("SampleNumber",paste("Random_",seq(1:n),sep=""))
    
    for (i in 1:numsam) {data.m$avg[i]<-mean(as.numeric(data.m[i,2:n+1])) }
    for (i in 1:numsam) {data.m$sem[i]<-sd(as.numeric(data.m[i,2:n+1]))/sqrt(n) } #create the error bar data
    
    #Fit model
    data.m2<-data.frame(cbind(seq(from=1,to=length(samplelist),by=1),data.m$avg,data.m$sem))
    names(data.m2)<-c("sampleNum","uniquePeaks","sem")
    
    #use self starting model for asymptotic data SSasymp
    
    #it follows the formula: Asym+(R0-Asym)*exp(-exp(lrc)*input)
    # input  a numeric vector of values at which to evaluate the model.
    # Asym  a numeric parameter representing the upper asymptotic value of the model (as input goes to Inf).
    # R0	a numeric parameter representing the response when input is zero. This argument is not in SSasympOrig, in which it is assumed to be zero.
    # lrc	a numeric parameter representing the natural logarithm of the rate constant.
    
    nlsfitSS <- nls(uniquePeaks ~ SSasymp(sampleNum, Asym, R0, lrc),
                  data=data.m2)
    
    new.sams = seq(1,200,1)
    pred <- predict(nlsfitSS,list(sampleNum=new.sams)) #predict the future values with increasing numbers of samples
    
    #check that the model is a good fit for the data
    (RSS.p <- sum(residuals(nlsfitSS)^2))  # Residual sum of squares
    ## [1] 267158552
    (TSS <- sum((data.m2$uniquePeaks - mean(data.m2$uniquePeaks))^2))  # Total sum of squares
    ## [1] 66071809548
    1 - (RSS.p/TSS)  # R-squared measure
    ## [1] 0.9959565 #this should be close to 1
    
    #extract key parameters from the model
    Asym_coef<-summary(nlsfitSS)$coefficients[1] #100% saturation
    R0_coef<-summary(nlsfitSS)$coefficients[2]
    lrc_coef<-summary(nlsfitSS)$coefficients[3]
    
    #predict sample number for 95% and 99% saturation
    sat95<-Asym_coef*0.95 #peaks at 95% saturation
    sat99<-Asym_coef*0.99 #peaks at 99% saturation
    SN_sat95<--(log((sat95-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) #number of samples to reach 95% saturation
    SN_sat99<--(log((sat99-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) #number of samples to reach 99% saturation
    curpc<-round((data.m2$uniquePeaks[numsam]/Asym_coef)*100) #current percentage saturation achieved
    
    pdf(opname)
    plot(x=data.m2$sampleNum,y=data.m2$uniquePeaks,xlab="#Samples", ylab="#Peaks (x1000)",main="Saturation",xaxt="n",yaxt="n",
         ylim=c(0,max(data.m)+50000),xlim=c(1,ceiling(numsam/maxsamps)*maxsamps),type='p',col='dodgerblue4',pch=20,bty='n')
    arrows(data.m2$sampleNum, data.m2$uniquePeaks-data.m2$sem, data.m2$sampleNum, data.m2$uniquePeaks+data.m2$sem, length=0.1, angle=90, code=3)
    axis(side=1, at=seq(0,ceiling(numsam/maxsamps)*maxsamps,1),pos=0, las=2)
    axis(side=2, at=seq(0,max(data.m)+50000,20000),las=2,labels=seq(0,max(data.m/1000)+50,20),pos=0)
    abline(h=Asym_coef,col="brown2",lwd=1) #adds the asymptote, saturation=100%
    lines(pred, col="dodgerblue4", lwd=1,lty="dashed") #plots the fitted model,converging on the asymptote
    segments(x0=c(SN_sat95,SN_sat99,numsam),y0=c(0,0),y1=c(sat95,sat99,data.m2$uniquePeaks[numsam]),col="darkgrey",lwd=1,lty="dashed") #plots the vertical lines at 95% and 99% saturation and one where the current number is
    text(x=c(SN_sat95,SN_sat99,numsam)-1,y=c(20000,20000),srt=45,labels=c("95%","99%",paste(as.character(curpc),"%",sep="")),
         col="darkgrey",font=1) #adds the labels to the vertical lines for the saturation values
    dev.off()
    
    return(data.m)
}

lsc.sat <- satfunc(n=100, maxsamps=100,lsc.samples,lsc.peaks,"results/PCSC1/LSCp.SaturationCurve.pdf")
gbm.sat <- satfunc(n=100, maxsamps=100,gbm.samples,gbm.peaks,"results/PCSC1/GBM.SaturationCurve.pdf")
pfa.sat <- satfunc(n=100, maxsamps=100,pfa.samples,pfa.peaks,"results/PCSC1/PFA.SaturationCurve.pdf")


lsc.sat.all <- satfunc(n=100, maxsamps=100,lsc.samples,lsc.peaks,"results/PCSC1/LSCp.All.SaturationCurve.pdf")
gbm.sat.all <- satfunc(n=100, maxsamps=100,gbm.samples,gbm.peaks,"results/PCSC1/GBM.All.SaturationCurve.pdf")
pfa.sat.all <- satfunc(n=100, maxsamps=100,pfa.samples,pfa.peaks,"results/PCSC1/PFA.All.SaturationCurve.pdf")


#####################################################################
#**********************Peaks and Reads************************
#####################################################################
alignmeta <- fread("data/Alignment.Stats.txt", header=T, sep="\t", stringsAsFactors = F, data.table=F)
alignmeta.pos <- subset(alignmeta, alignmeta$Group=="positive" )

df1.r <- subset(alignmeta.pos[,c(1,4)], alignmeta.pos$Cancer=="LSC")
df2.r <- subset(alignmeta.pos[,c(1,4)], alignmeta.pos$Cancer=="GBM")
df3.r <- subset(alignmeta.pos[,c(1,4)], alignmeta.pos$Cancer=="PFA")

df1.p <- subset(alignmeta.pos[,c(1,8)], alignmeta.pos$Cancer=="LSC")
df2.p <- subset(alignmeta.pos[,c(1,8)], alignmeta.pos$Cancer=="GBM")
df3.p <- subset(alignmeta.pos[,c(1,8)], alignmeta.pos$Cancer=="PFA")

colnames(df1.r) <- c("Sample","Value")
colnames(df2.r) <- c("Sample","Value")
colnames(df3.r) <- c("Sample","Value")
colnames(df1.p) <- c("Sample","Value")
colnames(df2.p) <- c("Sample","Value")
colnames(df3.p) <- c("Sample","Value")
df1.p$Value <- (-1*df1.p$Value)/1000
df2.p$Value <- (-1*df2.p$Value)/1000
df3.p$Value <- (-1*df3.p$Value)/1000
df1.r$Value <- (1*df1.r$Value)/1000000
df2.r$Value <- (1*df2.r$Value)/1000000
df3.r$Value <- (1*df3.r$Value)/1000000

df1.r$group <- rep("Reads", nrow(df1.r))
df2.r$group <- rep("Reads", nrow(df2.r))
df3.r$group <- rep("Reads", nrow(df3.r))
df1.p$group <- rep("Peaks", nrow(df1.p))
df2.p$group <- rep("Peaks", nrow(df2.p))
df3.p$group <- rep("Peaks", nrow(df3.p))

df1 <- rbind(df1.r, df1.p)
df2 <- rbind(df2.r, df2.p)
df3 <- rbind(df3.r, df3.p)

df1$rank <- c(rank(df1.r$Value),rank(df1.r$Value))
df2$rank <- c(rank(df2.r$Value),rank(df2.r$Value))
df3$rank <- c(rank(df3.r$Value),rank(df3.r$Value))

pdf(paste0("results/PCSC1/LSC.ReadsPeaks.pdf"))
ggplot(df1, aes(x = reorder(Sample, rank), y = Value, fill = group)) + 
  geom_bar(stat="identity", position="identity") +
  ylim(round(min(df1$Value)-1), round(max(df1$Value)+1) ) +
  scale_y_continuous(breaks = seq(round(min(df1$Value)-1), round(max(df1$Value)+1),5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+ coord_flip()
dev.off()


pdf(paste0("results/PCSC1/GBM.ReadsPeaks.pdf"))
ggplot(df2, aes(x = reorder(Sample, rank), y = Value, fill = group)) + 
  geom_bar(stat="identity", position="identity") +
  ylim(round(min(df2$Value)-1), round(max(df2$Value)+1) ) +
  scale_y_continuous(breaks = seq(round(min(df2$Value)-1), round(max(df2$Value)+1),5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15)) + coord_flip()
dev.off()


pdf(paste0("results/PCSC1/PFA.ReadsPeaks.pdf"))
ggplot(df3, aes(x = reorder(Sample, rank), y = Value, fill = group)) + 
  geom_bar(stat="identity", position="identity") +
  ylim(round(min(df3$Value)-1), round(max(df3$Value)+1) ) +
  scale_y_continuous(breaks = seq(-125,40,5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+ coord_flip()
dev.off()