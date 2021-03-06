#!/usr/bin/Rscript

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi
# Princess Margaret Cancer Centre - University Health Network, July 27, 2016
# Tested on R 3.3.1

# Takes correlations between open regions of chromatin based on DNaseI hypersensitivity signals
# Regions with high correlations are candidates for 3D interactions
# Performs association tests on each candidate & adjusts p-values
# Identifies transcription factor motifs which overlap open regions
# Produces interaction landscapes and motif tracks in PDF format

args <- commandArgs(trailingOnly = TRUE)

refMapDir <- args[1] # directory of mapped files
print(refMapDir)

outDir <- sub('/$', '', args[2]) # output directory
print(outDir)

anchor <- args[3] # anchor file (bed)
print(anchor)

bg <- args[4] # file with list of bg files
print(bg)

date <- as.character(args[5]) # timestamp for output pdf
print(date)

signalMatrixFile <- (args[6]) # signal data (if available)
print(signalMatrixFile)

window <- 500000 # flanking bps from anchor to search for open regions
rcut <- 0.5 # correlation threshold
pcut <- 0.05 # p-value threshold
qcut <- 0.05 # q-value threshold
corMethod <- "pearson" # correlation coefficient (pearson, spearman, kendall)
motifFile <- " " # a file with a list of motifs
fimo_pcut <- "0.0001" # p-value cut-off for FIMO
fimo_qcut <- "1" # q-value cut-off for FIMO
figures <- "n" # 'y' if you want interaction landscapes outputted
figureWidth <- "500000"# flanking bps from anchor to display on figure
zoom <- "0" # how many flanking bps for zoomed landscapes
colour <- "#bdd7e7,#6baed6,#3182bd,#08519c" # colours for the interaction plots
tracks <- "y" # 'y' if you want tracks outputted
sampleName <- " "
trackNumber <- 1
numSamples <- 1
assembly <- "hg38"



setwd("/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts/C3D/")
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Sushi"))
suppressMessages(library("data.table"))
suppressMessages(library("preprocessCore"))
suppressMessages(library("dynamicTreeCut"))

# From http://davetang.org/muse/2015/02/04/bed-granges/
# converts bed file to a granges object
bed_to_granges = function(file) {
  df <- read.table(file, header=F, stringsAsFactors=F)[,1:3]
  names(df) <- c('chr','start','end')
  gr <- with(df, GRanges(chr, IRanges(start, end)))
  return(gr)
}

# PRE: int[>=1] int[>=1] matrix(num)
# POST: vector(num)
# returns correlation and p-value of xth & yth rows of A
getCor = function(x, y, A) {
  if (sd(A[x,])==0 || sd(A[y,])==0) { # avoid div by 0
    return(c(0,NA))
  }
  corr <- cor.test(A[x,], A[y,], method=corMethod, alternative="two.sided")
  return(c(unname(corr$estimate), corr$p.value))
}

# PRE: df(chr1,start1,end1,chr2,start2,end2,score,q_Value) str int[>=1] int[>=1]
# POST: 0<=num<=1
# returns max score of int on chromosome chr between start & end
get_max_cor = function(int, chr, start, end) {
  filteredInt <- int[with(int, chrom1==chr & end1>=start & start1<=end &
                            chrom2==chr & end2>=start & start2<=end),]
  return(max(filteredInt$score))
}

# PRE: str int[>=1]
# POST: vect
# Breaks a string into a vector of strings separated at commas
commaSepStr_to_vector = function(commaSepStr, desiredLength) {
  vect <- unlist(strsplit(commaSepStr, split=","))
  # if the string didn't have desiredLength items,
  # then just make a vector with the first element repeated desiredLength times
  if (length(vect)!=desiredLength) {
    vect <- rep(vect[1], desiredLength)
  }
  return(vect)
}

# PRE: str
# POST: int or str
# converts strands into numbers for the Sushi package
strand_to_num = function(strand) {
  if (strand=="+") {
    return(1)
  } else if (strand=="-") {
    return(-1)
  } else {
    return(".")
  }
}

# PRE: int[>=0] vect[hex]
# POST: str
# Converts the xth element of pal (a palette) to RGB values
getRGB <- function(x, pal) {
  rgbValues <- col2rgb(pal[x])
  return(paste(rgbValues[1],rgbValues[2],rgbValues[3], sep=","))
}

# PRE: matrix or vector
# POST: vector
# Calculates the mean of each row
rowAvg = function(A) {
  if (is.matrix(A)) {
    return(rowMeans(A))
  } else {
    return(A)
  }
}

# if signalMatrixFile not provided, generate matrix from mapped files
if (signalMatrixFile=="n") {
  # list of files in refMapDir ending in .map.bed
  print("Signal matrix file is not present. it will be generated")
  refMapDir <-  sub('/$', '', refMapDir)
  refMapFiles <- list.files(refMapDir, pattern="*.map.bed", full.names=TRUE)
  # granges object with regions from reference file
  ref.bed <- bed_to_granges(refMapFiles[1])
  # format names of regions from reference file
  regionNames <- paste(as.character(seqnames(ref.bed)),":",start(ranges(ref.bed)),
                       "-",end(ranges(ref.bed)), sep="")

  cat("Merging Mapped Files...\n")
  # merge scores from bg files
  signals <- do.call(cbind, lapply(refMapFiles,
                                   function(f) read.table(f,header=FALSE, sep="\t")[,4]))
  rownames(signals) <- regionNames
  write.table(signals, paste(outDir,"signalMatrix.txt", sep="/"),
              row.names=T, col.names=F, quote=F) # write matrix to file
} else { # if signalMatrixFile provided, use that as the matrix
  # load matrix
  print("Using given signal matrix")
  signals <- as.matrix(read.table(signalMatrixFile, header=FALSE, row.names=1))
  regionNames <- rownames(signals)
  # create granges object from matrix file
  partsOfBed <- strsplit(regionNames, ":|-")
  chr <- sapply(partsOfBed, function(x) x[1])
  start <- sapply(partsOfBed, function(x) as.numeric(x[2]))
  end <- sapply(partsOfBed, function(x) as.numeric(x[3]))
  df <- data.frame(chr, start, end)
  ref.bed <- with(df, GRanges(chr, IRanges(start, end)))
}
#
# if (! is.numeric(signals)) { # if matrix has NAs
#   stop("There are null values in the signal matrix. Check mapped files.\n")
# }
#
# signals.norm <- normalize.quantiles(signals) # quantile normalize the signals
# rownames(signals.norm) <- regionNames
# # calculate distances between samples (1-correlation)
# crossGenomeCors <- as.dist(1-cor(signals.norm, method="pearson"))
# dendro <- hclust(crossGenomeCors) # hierarchically cluster the samples
# # run the dynamic tree cut algorithm to detect clusters
# clusters <- as.character(cutreeDynamic(dendro, minClusterSize=1,
#                                        distM=as.matrix(crossGenomeCors), deepSplit=4))
# # average normalized signal at each open region for each cluster
# avg.norm.signals <- do.call(cbind, lapply(unique(clusters),
#                                           function(p) rowAvg(signals.norm[,which(clusters==p)])))
#
# # vector of gene names in same order as anchor file
# genes <- read.table(anchor, header=F, stringsAsFactors=F)[,5]
# anchorStrands <- read.table(anchor, header=F, stringsAsFactors=F)[,4]
# # convert anchor to GRanges object
# anchor.bed <- bed_to_granges(anchor)
# # store motif occurances of TFs from motifFile from FIMO
# motifs <- data.frame(chr=c(),start=c(),end=c(), strand=c(), score=c(),
#                      p_Value=c(), q_Value=c(), name=c(), row=c())
# TFs <- c() # list of TFs
# if (motifFile!="") {
#   # load JASPAR ID-Name key
#   jasparTbl <- read.table("JASPAR_IDs.txt", header=T, stringsAsFactors=F)
#   # read each line of motif file
#   TFs <- scan(motifFile, what=character(), sep="\n")
#   counter <- 1 # initialize counter for rows
#   for (motifID in TFs) {
#     # read FIMO results from motif database
#     fileLoc <- paste("motifs/",motifID,".meme.output/fimo.txt.bed", sep="")
#     df <- read.table(fileLoc, header=F, stringsAsFactors=F)
#     names(df) <- c('chr','start','end', 'strand','score','p_Value','q_Value')
#     df$name <- rep(motifID, nrow(df)) # populate name column with TF name
#     df$row <- rep(counter, nrow(df)) # each TF gets a different row number
#     motifs <- rbind(motifs, df)
#     counter <- counter + 1
#   }
#   # filter motifs based on p & q-values & sort
#   motifs <- motifs[with(motifs, p_Value<=fimo_pcut & q_Value<=fimo_qcut),]
#   motifs <- motifs[with(motifs, order(chr, start)),]
#   # convert motifs to GRanges object
#   motifs.bed <- with(motifs, GRanges(chr, IRanges(start, end), name=name,
#                                      strand=strand, score=score, p_Value=p_Value,
#                                      q_Value=q_Value, row=row))
#   # filter out motifs not in DHSs
#   motifs.bed <- motifs.bed[unique(queryHits(findOverlaps(motifs.bed, ref.bed)))]
# } else {
#   motifs.bed <- GRanges()
# }
# # find promoters which overlap reference regions
# promoterOverlaps <- findOverlaps(anchor.bed, ref.bed)
# if (window != "genome") {
#   window <- as.numeric(window)
#   # find reference regions within window of each promoter
#   leftOverlaps <- findOverlaps(flank(anchor.bed, window), ref.bed)
#   rightOverlaps <- findOverlaps(flank(anchor.bed, window, start=FALSE), ref.bed)
# }
#
# # find motifs which overlap reference regions
# motifOverlaps <- findOverlaps(motifs.bed, ref.bed)
# # initialize list of interaction candidate
# anchorStats <- vector(mode="list", length=length(anchor.bed))
# file <- paste(outDir,"/results_",date,".txt", sep="")
# resultsHeader <- paste("COORD_1\tCOORD_1_Motifs\tCOORD_1_Strands\tCOORD_2\tCOORD_2_Motifs\tCOORD_2_Strands\tR_",
#                          corMethod,"\tGENE\tp_Value\tq_Value" , sep="")
# junk <- cat(resultsHeader, file=file, append=FALSE, sep="\n")
#
# for (r in 1:length(anchor.bed)) { # iterate through promoters
#   COORD_1 <- c() # regions in promoter
#   COORD_1_Motifs <- c() # TFs whose motifs overlap COORD_1
#   COORD_1_Strands <- c() # strands of the TFs
#   COORD_2 <- c() # distal regions from promoter
#   COORD_2_Motifs <- c() # TFs whose motifs overlap COORD_2
#   COORD_2_Strands <- c() # strands of the TFs
#   Correlation <- c() # correlation
#   p_Value <- c() # p-value
#   q_Value <- c() # q-value
#   GENE <- c() # gene name
#   # regions in rth promoter
#   regionsInPromoters <- subjectHits(promoterOverlaps[which(queryHits(promoterOverlaps)==r)])
#   if (window=="genome") { # search the entire genome
#     candidateRegions <- setdiff(1:length(ref.bed), regionsInPromoters)
#   } else {
#     regionsInLeft <- subjectHits(leftOverlaps[which(queryHits(leftOverlaps)==r)])
#     regionsInRight <- subjectHits(rightOverlaps[which(queryHits(rightOverlaps)==r)])
#     # regions in window
#     candidateRegions <- unique(c(regionsInLeft, regionsInRight))
#     candidateRegions <- candidateRegions[! candidateRegions %in% regionsInPromoters]
#   }
#   if (length(regionsInPromoters)==0 || length(candidateRegions)==0) next
#   # calculate correlations between open regions in promoter & distal regions
#   regionIndices <- cbind(rep(regionsInPromoters,each=length(candidateRegions)),
#                          rep(candidateRegions, length(regionsInPromoters)))
#   corStats <- apply(regionIndices, 1, function(x) getCor(x[1], x[2], avg.norm.signals))
#   # record regionNames of anchor DHSs
#   COORD_1 <- regionNames[regionIndices[,1]]
#   # motifs which overlap COORD_1 regions
#   COORD_1_Motifs <- sapply(regionIndices[,1], function (x) paste(motifs.bed[queryHits(motifOverlaps[which(subjectHits(motifOverlaps)==x)])]$name, collapse=","))
#   COORD_1_Motifs[COORD_1_Motifs==""] <- "NA"
#   # strand of each motif
#   COORD_1_Strands <- sapply(regionIndices[,1], function (x) paste(as.character(strand(motifs.bed[queryHits(motifOverlaps[which(subjectHits(motifOverlaps)==x)])])), collapse=","))
#   COORD_1_Strands[COORD_1_Strands==""] <- "NA"
#   # record regionNames of distal DHSs
#   COORD_2 <- regionNames[regionIndices[,2]]
#   # motifs which overlap COORD_2 regions
#   COORD_2_Motifs <- sapply(regionIndices[,2], function (x) paste(motifs.bed[queryHits(motifOverlaps[which(subjectHits(motifOverlaps)==x)])]$name, collapse=","))
#   COORD_2_Motifs[COORD_2_Motifs==""] <- "NA"
#   # strands of each motif
#   COORD_2_Strands <- sapply(regionIndices[,2], function (x) paste(as.character(strand(motifs.bed[queryHits(motifOverlaps[which(subjectHits(motifOverlaps)==x)])])), collapse=","))
#   COORD_2_Strands[COORD_2_Strands==""] <- "NA"
#   # record correlations, p-values & q-values
#   Correlation <- corStats[1,]
#   p_Value <- corStats[2,]
#   GENE <- rep(genes[r], nrow(regionIndices)) # gene name
#   q_Value <- p.adjust(p_Value, method="BH")
#   interactionCandidates <- data.frame(COORD_1, COORD_1_Motifs, COORD_1_Strands,
#                                       COORD_2, COORD_2_Motifs, COORD_2_Strands,
#                                       Correlation, GENE, p_Value, q_Value)
#   if (nrow(interactionCandidates)>0) { # stop if genomeStats is empty
#     junk <- write.table(interactionCandidates, file=file, sep="\t", col.names=F,row.names=F,quote=F, append=TRUE)
#   }
#   anchorStats[[r]] <- interactionCandidates
#   cat('\rCalculating Correlations: Processed anchor', r, 'of', length(anchor.bed))
# }
# # concatenate anchorStats data.frames
# genomeStats <- rbindlist(anchorStats)
# # filter genomeStats
# genomeStats <- genomeStats[with(genomeStats, Correlation>=rcut & p_Value<=pcut & q_Value<=qcut),]
# # get ref.bed indices corresponding to COORD_1 & COORD_2
# regionIndices.Coord1 <- sapply(genomeStats$COORD_1,
#                                function(x) which(regionNames==x))
# regionIndices.Coord2 <- sapply(genomeStats$COORD_2,
#                                function(x) which(regionNames==x))
# if (length(regionIndices.Coord1)<=0 && figures=="y") {
#   cat("No interaction candidates passed the filters; No figures generated\n")
# }
# # get widths for figures & tracks
# figureWidth <- as.numeric(commaSepStr_to_vector(figureWidth, length(anchor.bed)))
#
# # graphics --------------------------------------------------------------------
# if (figures=="y" && length(regionIndices.Coord1)>0) {
#   # store motifs.bed into motifs
#   chr <- as.character(seqnames(motifs.bed))
#   start <- start(ranges(motifs.bed))
#   end <- end(ranges(motifs.bed))
#   # convert motifs strands to numbers for Sushi package
#   strand <- sapply(as.character(strand(motifs.bed)), function(x) strand_to_num(x))
#   name <- mcols(motifs.bed)$name
#   score <- mcols(motifs.bed)$score
#   p_Value <- mcols(motifs.bed)$p_Value
#   q_Value <- mcols(motifs.bed)$q_Value
#   row <- mcols(motifs.bed)$row
#
#   motifs <- data.frame(chr, start, end, strand, score, p_Value, q_Value, name, row)
#   # create data frame for anchors
#   chr <- as.character(seqnames(anchor.bed))
#   start <- start(ranges(anchor.bed))
#   end <- end(ranges(anchor.bed))
#   name <- sapply(strsplit(genes, "_"), function(x) head(x, n=1))
#   strand <- sapply(anchorStrands, function(x) strand_to_num(x))
#   score <- rep(".", length(name))
#   p_Value <- q_Value <- rep(0, length(name))
#   # assign row number to anchors
#   if (nrow(motifs)==0) {
#     row <- rep(1, length(name))
#   } else {
#     row <- rep(max(motifs$row)+1, length(name))
#   }
#
#   anchor.df <- data.frame(chr,start,end, strand, score, p_Value, q_Value, name, row)
#   # bind anchors to motifs, sort and assign each TF a colour
#   motifs <- rbind(motifs, anchor.df)
#
#   # color of each motif in the figure
#   color <- rep("#5e3c99", nrow(motifs)) # colour motifs purple if no strand info given
#   color[which(motifs$strand==1)] <- "#0571b0" # + strand
#   color[which(motifs$strand==-1)] <- "#ca0020" # - strand
#   motifs$color <- as.character(color) # convert colors to characters
#   motifs <- motifs[with(motifs, order(chr, start)),] # sort motifs
#
#   # extract colours for interactions
#   colour <- commaSepStr_to_vector(colour, 4)
#   # make interaction data frame
#   chrom1 <- as.character(seqnames(ref.bed[regionIndices.Coord1]))
#   start1 <- start(ranges(ref.bed[regionIndices.Coord1]))
#   end1 <- end(ranges(ref.bed[regionIndices.Coord1]))
#   chrom2 <- as.character(seqnames(ref.bed[regionIndices.Coord2]))
#   start2 <- start(ranges(ref.bed[regionIndices.Coord2]))
#   end2 <- end(ranges(ref.bed[regionIndices.Coord2]))
#   color <- rep("black", length(genomeStats$q_Value))
#   # make colours based on q-values
#   color[which(genomeStats$q_Value>0.05)] <- colour[1]
#   color[which(0.01<genomeStats$q_Value & genomeStats$q_Value<=0.05)] <- colour[2]
#   color[which(0.001<genomeStats$q_Value & genomeStats$q_Value<=0.01)] <- colour[3]
#   color[which(genomeStats$q_Value<=0.001)] <- colour[4]
#   interactions.bedpe <- data.frame(chrom1, start1, end1, chrom2, start2, end2,
#                                    score=genomeStats$Correlation, q_Value=genomeStats$q_Value, color)
#   # extract zoom lengths
#   zoom <- as.numeric(commaSepStr_to_vector(zoom, length(anchor.bed)))
#   zoomOverWidth <- which(zoom>figureWidth)
#   zoom[zoomOverWidth] <- figureWidth[zoomOverWidth] # truncate values over window to window
#   noZoom <- which(zoom<=0) # figures with no zoom
#
#   corName <- paste(toupper(substr(corMethod, 1, 1)),
#                    substr(corMethod, 2, nchar(corMethod)), sep="")
#   # make pdf file with figures
#   pdf(file=paste(outDir,"/figures_",date,".pdf", sep=""), height=11, width=8.5)
#   for (r in 1:length(anchor.bed)) { # make figure for each anchor
#     # set up the region
#     chrom = as.character(seqnames(anchor.bed[r]))
#     chromstart <- start(ranges(anchor.bed[r])) - figureWidth[r]
#     chromend <- end(ranges(anchor.bed[r])) + figureWidth[r]
#     maxCor <- get_max_cor(interactions.bedpe, chrom,chromstart,chromend)
#     if (r %in% noZoom) { # figure with no zoom
#       layout(matrix(c(1,1,1,1,2,2), 3,2, byrow=TRUE))
#       par(mar=c(4,6,1,8), oma=rep(0.75,4), xpd=TRUE)
#       # landscape plot
#       intLandscape = plotBedpe(interactions.bedpe, chrom,chromstart,chromend,
#                                heights=interactions.bedpe$score,
#                                lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
#                                plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
#       labelgenome(chrom,chromstart,chromend, n=5,scale="bp",
#                   chromline=0.25,scaleline=0.25)
#       legend("right",inset=-0.16,title=expression(bold("q-value")),
#              legend=c(expression(""<="1"),expression(""<="0.05"),
#                       expression(""<="0.01"),expression(""<="0.001")),
#              col=colour,lty=1,lwd=2.5,text.font=2,cex=1.25)
#       axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
#       mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
#       # tracks for anchors & motifs
#       par(mar=c(1,6,1,8))
#       motifsTrack = plotBed(beddata=motifs,chrom=chrom,chromstart=chromstart,chromend=chromend,
#                             rownumber=motifs$row, type="region", color=motifs$color,row="given",
#                             rowlabels=c(TFs, "Anchors"), rowlabelcol="black", rowlabelcex=1.25)
#       if (motifFile!="") {
#         legend("bottomright",inset=c(-0.16,0.1),title=expression(bold("Strand")),
#                legend=c("forward","reverse","unknown"),fill=c("#0571b0","#ca0020","#5e3c99"),
#                border=c("#0571b0","#ca0020","#5e3c99"),cex=1.25)
#       }
#     } else { # figure with zoom
#       layout(matrix(c(1,1,1,1,2,2,2,2,3,3), 5,2, byrow=TRUE))
#       par(mar=c(4,6,1,8), oma=rep(0.75,4), xpd=TRUE)
#       # landscape plot
#       intLandscape = plotBedpe(interactions.bedpe, chrom,chromstart,chromend,
#                                heights=interactions.bedpe$score,
#                                lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
#                                plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
#       labelgenome(chrom,chromstart,chromend, n=5,scale="bp",
#                   chromline=0.25,scaleline=0.25)
#       legend("right",inset=-0.16,title=expression(bold("q-value")),
#              legend=c(expression(""<="1"),expression(""<="0.05"),
#                       expression(""<="0.01"),expression(""<="0.001")),
#              col=colour,lty=1,lwd=2.5,text.font=2,cex=1.25)
#       axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
#       mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
#       # zoomed region
#       par(mar=c(1,6,1,8))
#       regionZoom <- c(start(ranges(anchor.bed[r]))-zoom[r], end(ranges(anchor.bed[r]))+zoom[r])
#       zoomsregion(regionZoom,extend=c(-1,0.13),wideextend=0,offsets=c(0,0))
#       maxCor <- get_max_cor(interactions.bedpe, chrom,regionZoom[1],regionZoom[2])
#       # zoomed landscape plot
#       intLandscapeZoom = plotBedpe(interactions.bedpe, chrom,chromstart=regionZoom[1],chromend=regionZoom[2],
#                                    heights=interactions.bedpe$score,
#                                    lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
#                                    plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
#       zoombox()
#       labelgenome(chrom,chromstart=regionZoom[1],chromend=regionZoom[2], n=5,scale="bp",
#                   chromline=0.25,scaleline=0.25)
#       axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
#       mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
#       # tracks for anchors & motifs
#       par(mar=c(1,6,1,8))
#       motifsTrack = plotBed(beddata=motifs,chrom=chrom,chromstart=regionZoom[1],chromend=regionZoom[2],
#                             rownumber=motifs$row, type="region", color=motifs$color,row="given",
#                             rowlabels=c(TFs, "Anchors"), rowlabelcol="black", rowlabelcex=1.25)
#       if (motifFile!="") {
#         legend("bottomright",inset=c(-0.16,0.1),title=expression(bold("Strand")),
#                legend=c("forward","reverse","unknown"),fill=c("#0571b0","#ca0020","#5e3c99"),
#                border=c("#0571b0","#ca0020","#5e3c99"),cex=1.25)
#       }
#     }
#     mtext(paste("Figure ",r,": Interaction Landscape of ",genes[r], sep=""), side=1)
#     cat('\rGenerating Figures: Created figure', r, 'of', length(anchor.bed))
#   }
#   cat("\n")
#   turnDevOff <- dev.off()
# }
#
# # tracks ----------------------------------------------------------------------
# if (tracks=="y") {
#   if (length(TFs) > 0) {
#     # format the motif occurances so the UCSC genome browser can read it
#     motifBeds <- vector(mode="list", length=length(TFs))
#     motifBeds <- lapply(TFs, function(x) motifs.bed[which(mcols(motifs.bed)$name==x)])
#     names(motifBeds) <- TFs
#     motifBeds <- lapply(motifBeds, function(b)
#                                    data.frame(as.character(seqnames(b)), start(b), end(b), ".", "0", as.character(strand(b))))
#     # the palette of colours for the motif tracks, each TF gets a colour between green and magenta
#     motifPal <- colorRampPalette(c("#4dac26", "#d01c8b"))(length(TFs))
#
#     for (t in 1:length(TFs)) { # write motifs to track if motifFile was passed
#       file <- paste(outDir,"/",TFs[t],".bed", sep="")
#       trackHeader <- paste("track name=\"",TFs[t],
#                            "\" description=\" \" visibility=1 color=",getRGB(t, motifPal)," db=", assembly, sep="")
#       junk <- cat(trackHeader, file=file, append=FALSE, sep="\n")
#       junk <- write.table(motifBeds[[t]], file=file, sep="\t", col.names=F,row.names=F,quote=F, append=TRUE)
#     }
#   }
#
#
#   trackPal <- colorRampPalette(c("blue", "purple", "red", "orange"))(numSamples)
#   for (r in 1:length(anchor.bed)) {
#     if (is.null(anchorStats[[r]])) next
#     file <- paste(outDir,"/",genes[r],".anchor.bedGraph", sep="")
#     # create & print the browser header to a file
#     browserHeader <- paste("browser position ",as.character(seqnames(anchor.bed[r])),":",
#                            start(anchor.bed[r])-figureWidth[r],"-",end(anchor.bed[r])+figureWidth[r],
#                            "\nbrowser hide all\nbrowser pack refGene\nbrowser full altGraph", sep="")
#     junk <- cat(browserHeader, file=file, append=FALSE, sep="\n")
#     # create & print the track header for the anchor
#     trackHeader <- paste("track type=bedGraph name=\"Anchor\" description=\"",gsub("_.*", "", genes[r]),
#                          "\" visibility=full color=0,0,0 altColor=0,0,0 viewLimits=0:1 autoScale=off gridDefault=on db=", assembly, sep="")
#     junk <- cat(trackHeader, file=file, append=TRUE, sep="\n")
#     # write the anchor track
#     anchorLine <- paste(as.character(seqnames(anchor.bed[r])),start(anchor.bed[r]),end(anchor.bed[r]),1, sep="\t")
#     junk <- cat(anchorLine, file=file, append=TRUE, sep="\n")
#
#     file <- paste(outDir,"/",genes[r],".bedGraph", sep="")
#     # create & print the track header for the distal DHSs
#     trackHeader <- paste("track type=bedGraph name=\"",sampleName,"\" description=\" \" visibility=full color=",getRGB(trackNumber, trackPal),
#                          " altColor=0,0,0 viewLimits=0:1 autoScale=off gridDefault=on yLineMark=",rcut," yLineOnOff=on db=", assembly, sep="")
#     junk <- cat(trackHeader, file=file, append=FALSE, sep="\n")
#     # create a temporary data frame for storing distal DHSs & their correlations with the anchor
#     tmp.df <- data.frame(distal=anchorStats[[r]]$COORD_2, corr=anchorStats[[r]]$Correlation)
#     # if a distal DHS appears more than once, pick the one with the highest correlation
#     # this will only happen if the promoter has multiple open regions
#     tmp.df <- tmp.df[with(tmp.df, order(distal, -corr)), ]
#     tmp.df <- tmp.df[ !duplicated(tmp.df$distal), ]
#     tmp.df <- tmp.df[tmp.df$corr>=0, ] # filter out negative correlations
#     # split the distal DHSs by chr, start & stop
#     distalDHS <- strsplit(as.character(tmp.df$distal), ":|-")
#     chr <- sapply(distalDHS, function(x) x[1])
#     start <- sapply(distalDHS, function(x) x[2])
#     end <- sapply(distalDHS, function(x) x[3])
#     correlation <- round(tmp.df$corr,2) # round correlations to 2 decimals
#     # write track to file
#     track <- data.frame(chr, start, end, correlation)
#     junk <- write.table(track, file=file, sep="\t", col.names=F,row.names=F,quote=F, append=TRUE)
#     cat('\rGenerating Tracks: Created track', r, 'of', length(anchor.bed))
#   }
#   cat("\n")
# }
# cat("Done\n")
