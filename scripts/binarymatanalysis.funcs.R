lawhoops = c('#371377','#7700FF','#9E0142',"#CB2314",'#FF0080', '#DC494C',
             '#F88D51', "#FAD510","#FFFF5F",'#88CFA4','#238B45',"#02401B",
             "#0AD7D3", "dodgerblue", "#046C9A", "#273046", "#A2A475", "#354823", 'grey35', "#1E1E1E")

#https://github.com/ManchesterBioinference/Scasat/blob/master/ScAsAT_functions_Buenrostro.ipynb
getJaccardDist <- function(cdBinary){
  
  if(colnames(cdBinary[,2:3])[1] == 'start' && colnames(cdBinary[,2:3])[2] == 'end'){
    SingleCell.Binary <- cdBinary[,4:(dim(cdBinary)[2])]
  }
  else
    SingleCell.Binary <- cdBinary
  
  SingleCell.Binary.Jaccard <- jaccard(as.matrix(SingleCell.Binary))
  
  return(SingleCell.Binary.Jaccard)
}


plotMDS <- function(cdBinary, k, groups=NULL, cellName, ret.val=FALSE, text.label=FALSE, 
                    title=""){
  
  if(missing(ret.val)){
    ret.val = FALSE
  }
  if(text.label==TRUE & missing(cellName)){
    stop("ERROR: Please give cellName")
  }
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  fit <- cmdscale(as.dist(SingleCell.Binary.Jaccard),eig=TRUE, k=k)
  
  if(is.null(groups)){
    df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard))
    p <- ggplot(df, aes_string(x="x",y ="y"))
  }
  else{
    df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard), Batch=groups)
    p <- ggplot(df, aes_string(x="x",y ="y", color="Batch"))+ 
      ##scale_colour_hue(l=40) + 
      scale_fill_manual(values=lawhoops) 
  }
  
  
  p <- p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p <- p + geom_point(size = 4)
  p <- p + xlab("Coordinate 1") 
  p <- p + ylab("Codinate 2")+
    theme_light(base_size=24) +
    theme(strip.background = element_blank(),
          panel.border     = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  if(text.label==TRUE)
    p <- p + geom_text(data=df,aes(label=cellName),alpha=0.5,size=4, vjust=1,hjust=0.5,angle=45)
  
  #p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
  #legend.text = element_text(size = 14), legend.title = element_text(size = 14))
  
  print(p)
  
  if(ret.val == TRUE)
    return(fit)
}


plotMDSClust <- function(cdBinary, k, Clusters=NULL, cellName, ret.val=FALSE, text.label=FALSE, 
                         title=""){
  
  if(missing(k)){
    stop("ERROR: Number of Clusters \"k\" is missing")
  }
  if(missing(ret.val)){
    ret.val = FALSE
  }
  if(text.label==TRUE & missing(cellName)){
    stop("ERROR: Please give cellName")
  }
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  fit <- cmdscale(as.dist(SingleCell.Binary.Jaccard),eig=TRUE, k=k)
  
  if(is.null(Clusters)){
    df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard))
    p <- ggplot(df, aes_string(x="x",y ="y"))
  }
  else{
    df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard), Clusters=Clusters)
    p <- ggplot(df, aes_string(x="x",y ="y", color="Clusters"))
  }
  
  
  p <- p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p <- p + geom_point(size = 4)
  p <- p + xlab("Coordinate 1") 
  p <- p + ylab("Coordinate 2")+
    theme_light(base_size=24) +
    theme(strip.background = element_blank(),
          panel.border     = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  if(length(levels(Clusters)) < 12)
    p <- p + scale_fill_manual(values=lawhoops)        
  if(text.label==TRUE)
    p <- p + geom_text(data=df,aes(label=cellName),alpha=0.5,size=4, vjust=1,hjust=0.5,angle=45)
  
  #p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
  #legend.text = element_text(size = 14), legend.title = element_text(size = 14))
  
  print(p)
  
  if(ret.val == TRUE)
    return(fit)
}


plotPCAJaccard <- function(cdBinary, groups=NULL, ret.val=FALSE , cellName, text.label=FALSE, title){
  
  
  if(missing(ret.val)){
    ret.val = FALSE
  }
  if(text.label==TRUE & missing(cellName)){
    stop("ERROR: Please give cellName")
  }
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  FinalPCAData <- t(SingleCell.Binary.Jaccard)
  PCx=1
  PCy=2
  pcaPRComp <- prcomp(FinalPCAData)
  percentVar <- pcaPRComp$sdev^2/sum(pcaPRComp$sdev^2)
  
  if(is.null(groups)){
    df<-data.frame(PCX=pcaPRComp$x[,PCx],PCY=pcaPRComp$x[,PCy])
    p1<-ggplot(df, aes_string(x="PCX",y ="PCY"))        
  }
  else{
    df<-data.frame(PCX=pcaPRComp$x[,PCx],PCY=pcaPRComp$x[,PCy], groups=groups)
    p1<-ggplot(df, aes_string(x="PCX",y ="PCY", color="groups"))             
  }          
  if(missing(title))
    p1<-p1+ggtitle("PCA with Jacard Matrix")
  else
    p1<-p1+ggtitle(title)
  
  p1<-p1+geom_point(size = 3)
  p1<-p1+xlab(paste("PC",PCx,": ", round(percentVar[PCx] * 100), "% variance"))
  p1<-p1+ylab(paste("PC",PCy,": ", round(percentVar[PCy] * 100), "% variance"))
  if(text.label == TRUE){
    p1<-p1+geom_text(data=df,aes(label=cellName),
                     alpha=0.5,size=4, vjust=1,hjust=0.5,angle=45, 
                     color="black")
  }
  
  p1<-p1 +
    theme_light(base_size=15)+     
    theme(plot.title = element_text(hjust = 0.5),
          panel.border     = element_blank())
  
  if(length(levels(groups)) < 12)
    p1 <- p1 + scale_fill_manual(values=lawhoops)
  
  print(p1)
  
  if(ret.val == TRUE)
    return(pcaPRComp)
}

plotMultiplePCAJaccard <- function(cdBinary, nPCAToDisplay, groups=NULL, ret.val=FALSE , text.label=FALSE, title=""){
  
  if(missing(nPCAToDisplay)){
    stop("ERROR: Number of PCA's \"nPCAToDisplay\" to display is missing")
  }
  if(missing(ret.val)){
    ret.val = FALSE
  }
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  FinalPCAData <- t(SingleCell.Binary.Jaccard)
  PCx=1
  PCy=2
  pcaPRComp <- prcomp(FinalPCAData)
  
  nmax = nPCAToDisplay
  plotTitle = paste0('Plotting first ',nmax,' PCAs')
  
  txt1 <- paste("Percent_PC_Var_onfirst",nmax,"PCs",sep="")
  pca_var = pcaPRComp$sdev ^ 2
  pca_var_percent <- 100 * pca_var / sum(pca_var)
  pca_var_percent_first10 <- NA * pca_var
  pca_var_percent_first10[1:nmax] <- 100 * pca_var[1:nmax] / sum(pca_var[1:nmax])
  
  pca_var_out <- data.frame(round(pca_var,3),round(pca_var_percent,1),
                            round(pca_var_percent_first10,1))
  rownames(pca_var_out) <- colnames(pcaPRComp$x)
  colnames(pca_var_out) <- c("PC_Var","PC_Var_percent",txt1)
  
  
  if(is.null(groups)){
    df <- as.data.frame(pcaPRComp$x)
    p <- ggpairs(df, columns=1:nPCAToDisplay, upper=list(continuous="points"), 
                 title=plotTitle,              
                 columnLabels = as.character(paste0(colnames(df[,1:nPCAToDisplay]), ' : ', 
                                                    pca_var_out$PC_Var_percent[1:nPCAToDisplay], '% variance')))+
      theme_light(base_size=15)+     
      theme(plot.title = element_text(hjust = 0.5))        
  }
  else{
    df <- as.data.frame(pcaPRComp$x)
    df$Cell=groups        
    p <- ggpairs(df, columns=1:nPCAToDisplay, upper=list(continuous="points"), 
                 title=plotTitle, 
                 mapping = aes_string(color="Cell"),
                 legend = c(1,nPCAToDisplay),
                 columnLabels = as.character(paste0(colnames(df[,1:nPCAToDisplay]), ' : ', 
                                                    pca_var_out$PC_Var_percent[1:nPCAToDisplay], '% variance')))+
      theme_light(base_size=15)+     
      theme(plot.title = element_text(hjust = 0.5))
  }          
  
  print(p)
}

plotVarExplained <- function(cdBinary, nPCAToDisplay, groups=NULL, ret.val=FALSE , text.label=FALSE, title=""){
  
  if(missing(nPCAToDisplay)){
    stop("ERROR: Number of PCA's \"nPCAToDisplay\" to display is missing")
  }
  if(missing(ret.val)){
    ret.val = FALSE
  }
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  FinalPCAData <- t(SingleCell.Binary.Jaccard)
  PCx=1
  PCy=2
  pcaPRComp <- prcomp(FinalPCAData)
  
  nmax = nPCAToDisplay
  plotTitle = paste0('Plotting first ',nmax,' PCAs')
  
  txt1 <- paste("Percent_PC_Var_onfirst",nmax,"PCs",sep="")
  pca_var = pcaPRComp$sdev ^ 2
  pca_var_percent <- 100 * pca_var / sum(pca_var)
  pca_var_percent_first10 <- NA * pca_var
  pca_var_percent_first10[1:nmax] <- 100 * pca_var[1:nmax] / sum(pca_var[1:nmax])
  
  pca_var_out <- data.frame(round(pca_var,3),round(pca_var_percent,1),round(pca_var_percent_first10,1))
  rownames(pca_var_out) <- colnames(pcaPRComp$x)
  colnames(pca_var_out) <- c("PC_Var","PC_Var_percent","PC_Var_n_percent")
  
  
  p <- ggplot(pca_var_out[1:nmax,]) + #scale_colour_hue(l=40) 
  p <- p + geom_line(aes(x=c(1:nmax), y = PC_Var_percent, color="PC_Var_percent"), size=1.5)+
    geom_point(aes(x=c(1:nmax), y = PC_Var_percent, color="PC_Var_percent"), size = 2)+
    geom_line(aes(x=c(1:nmax), y = PC_Var_n_percent, color="PC_Var_n_percent"), size=1.5)+
    geom_point(aes(x=c(1:nmax), y = PC_Var_n_percent, color="PC_Var_n_percent"), size = 2)+
    ggtitle("Total variance by the PCAs considered")+
    xlab(paste("PC1 to PC10"))+
    ylab(paste("% of Variance"))+     
    theme_light(base_size=15)+     
    theme(plot.title = element_text(hjust = 0.5),
          panel.border     = element_blank())+
    scale_x_continuous(breaks=c(1:nmax))+
    scale_fill_manual(values=lawhoops)
  
  print(p)
}

plot_tSNE <- function(cdBinary, nPCAToUSE, groups=NULL, cellName, perplexity_division, ret.val=FALSE , text.label=FALSE, title=""){
  
  if(missing(nPCAToUSE)){
    stop("ERROR: Number of PCA's \"nPCAToUSE\" to use is missing")
  }
  if(missing(ret.val)){
    ret.val = FALSE
  }
  
  if(text.label==TRUE & missing(cellName)){
    stop("ERROR: Please give cellName")
  }
  
  if(missing(perplexity_division)){
    stop("ERROR: Please enter the number with which you want to divide the dimension of your data for perplexity setting")
  }
  
  
  SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
  FinalPCAData <- t(SingleCell.Binary.Jaccard)
  
  pcaPRComp <- prcomp(FinalPCAData)
  
  
  #rtsne_pca_out <- Rtsne(as.matrix(pcaPRComp$x[,1:nPCAToUSE]), perplexity = dim(pcaPRComp$x)[1]/perplexity_division, check_duplicates = FALSE)
  rtsne_pca_out <- Rtsne(as.matrix(pcaPRComp$x[,1:nPCAToUSE]), perplexity = perplexity_division, 
                         check_duplicates = FALSE, pca=FALSE, theta=0.01, max_iter=3000)
  
  if(is.null(groups)){
    if(ret.val==TRUE)
      df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], cellName=cellName)
    else
      df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2])
    p1 <- ggplot(df, aes_string(x="X",y ="Y"))
  }
  else{
    if(ret.val==TRUE)
      df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], Cell=colnames(SingleCell.Binary.Jaccard), Batch=groups, cellName=cellName)
    else
      df<-data.frame(X=rtsne_pca_out$Y[,1],Y=rtsne_pca_out$Y[,2], Cell=colnames(SingleCell.Binary.Jaccard), Batch=groups)
    p1 <- ggplot(df, aes_string(x="X",y ="Y", color="Batch"))
  }
  
  
  p1<-p1 + ggtitle("t-SNE plot")
  p1<-p1 + geom_point(size = 2) 
  p1<-p1 + xlab(paste("Dim-1"))
  p1<-p1 + ylab("Dim-2")+
    theme_light(base_size=15) +
    theme(strip.background = element_blank(),
          panel.border     = element_blank(),
          plot.title = element_text(hjust = 0.5))
  if(length(levels(groups)) < 12)
    p1 <- p1 + scale_fill_manual(values=lawhoops)
  if(text.label==TRUE)
    p1<-p1 + geom_text(data=df,aes(label=cellName),alpha=0.5,size=4, vjust=1,hjust=0.5,angle=45)
  print(p1)
  return(df)
}


getDiffAccessInformationGain <- function(cdBinary, groups=NULL){
  
  if (length(levels(groups)) != 2) {
    stop(paste("ERROR: wrong number of levels in the grouping factor (", 
               paste(levels(groups), collapse = " "), "), but must be two.", 
               sep = ""))
  }
  if (is.null(groups)){
    stop("ERROR: groups factor is not provided")
  }
  
  
  SingleCell.Group1.CellNames <- names(groups[groups==levels(groups)[1]])
  SingleCell.Group2.CellNames <- names(groups[groups==levels(groups)[2]])
  
  SingleCell.Group1.Binary <- cdBinary[,SingleCell.Group1.CellNames]
  SingleCell.Group2.Binary <- cdBinary[,SingleCell.Group2.CellNames]
  
  CellType <- data.frame(CellType=c(rep(levels(groups)[1], ncol(SingleCell.Group1.Binary)), rep(levels(groups)[2], ncol(SingleCell.Group2.Binary)) ) )
  
  SingleCell.Group1VsGroup2 <- cbind(SingleCell.Group1.Binary,SingleCell.Group2.Binary)
  SingleCell.Group1VsGroup2 <- t(SingleCell.Group1VsGroup2)
  #SingleCell.Group1VsGroup2 <- cbind(SingleCell.Group1VsGroup2,CellType)   
  
  
  dataDim = (dim(SingleCell.Group1VsGroup2)[2])
  #dataDim = 5
  
  information.gain = vector(mode="numeric", length=(ncol(SingleCell.Group1VsGroup2)))
  SingleCell.Group1VsGroup2.res <- data.frame(Chr = cdBinary[1:dataDim,1],
                                              Start = cdBinary[1:dataDim,2],
                                              end = cdBinary[1:dataDim,3])        
  
  
  for(i in 1:dataDim)
  {
    jointData = cbind(SingleCell.Group1VsGroup2[,i], CellType)
    
    gain <- calcEntropy(jointData[,2]) - 
      (sum(table(jointData)[c(1,3)])/length(jointData[,2]))*calcEntropy(jointData[jointData[,1]==0,]) - 
      (sum(table(jointData)[c(2,4)])/length(jointData[,2]))*calcEntropy(jointData[jointData[,1]==1,])
    SingleCell.Group1VsGroup2.res$information.gain[i] <- gain        
  }
  
  
  SingleCell.Group1.rawAvg <- rowSums(SingleCell.Group1.Binary==1)/dim(SingleCell.Group1.Binary)[2]
  SingleCell.Group2.rawAvg <- rowSums(SingleCell.Group2.Binary==1)/dim(SingleCell.Group2.Binary)[2]
  
  SingleCell.Group1VsGroup2.res$rawMeanGroup1 <- SingleCell.Group1.rawAvg
  SingleCell.Group1VsGroup2.res$rawMeanGroup2 <- SingleCell.Group2.rawAvg
  SingleCell.Group1VsGroup2.res$log2FoldChange <- log2(SingleCell.Group1VsGroup2.res$rawMeanGroup1/
                                                         SingleCell.Group1VsGroup2.res$rawMeanGroup2)
  
  
  #write.csv(SingleCell.Group1VsGroup2[order(SingleCell.Group1VsGroup2$information.gain, decreasing=TRUE),], 
  #          paste0(levels(groups)[1],'_vs_',levels(groups)[2],'InformationGain.csv'), row.names=FALSE)
  
  SingleCell.Group1VsGroup2.res <- SingleCell.Group1VsGroup2.res[!is.nan(SingleCell.Group1VsGroup2.res$information.gain) &
                                                                   !is.na(SingleCell.Group1VsGroup2.res$information.gain),]
  
  return(SingleCell.Group1VsGroup2.res)    
}



calcEntropy <- function(data){
  freqs <- table(data)/sum(table(data))
  return(-sum(freqs * log(freqs)))
}