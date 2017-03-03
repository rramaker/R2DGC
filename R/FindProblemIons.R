#' This function scans over a range of potential ions and looks for ions that are not present in any peak or are common enough to decrease alignment quality.
#' Can use as input to the consensus align function to avoid including these ions during alignment to speed up processing time and improve alignments


#' @param inputFile The file path of a representative chromatof file to use in searching for ions to filter
#' @param possibleIons A numeric vector of possible ions to search. Make sure each ion listed is present in the input file. Defaults to 70 through 600.
#' @param numCores The number of cores to use for parallel processing. Defaults to 1
#' @param absentIonThreshold Numeric indicating the fraction of total ion intensity an ion has to greater than in at least one peak to not be filtered as an absent ion. Defaults to 0.01.
#' @param commonIonThreshold Numeric indicating the number of standard deviations below the mean an ion has to decrease the global metabolite similarity score to be filtered as a common ion. Defaults to 2.
#' @param plotData Boolean. If true, relative ion impact scores will be plotted.

#' @return Two column data frame identifying filtered ions and the reason they were filtered (absent or common). If plotData is TRUE, plots common ion scores. Y-axis is the z-scored number of pairwise metabolite comparisons with a similarity score greater than 50. X-axis is the ion with the filtered ions labeled in red.
#' @examples
#' FindProblemIons(inputFile=system.file("extdata", "SampleA.txt", package="R2DGC"),
#'     possibleIons = 70:78)
#'
#' @import parallel
#' @import graphics
#' @export

FindProblemIons<-function(inputFile, possibleIons=c(70:600), numCores=1, absentIonThreshold=0.01, commonIonThreshold=2, plotData=T){

  #Make empty vector to store similarity score impact of dropping each ion
  Sum50<-c()

  #Read in and format input file
  currentRawFile<-read.table(inputFile, header=T, sep="\t", fill=T, quote="",strip.white = T, stringsAsFactors = F)
  currentRawFile[,4]<-as.character(currentRawFile[,4])
  currentRawFile<-currentRawFile[which(!is.na(currentRawFile[,3])&nchar(currentRawFile[,4])!=0),]
  currentRawFile[,2]<-as.character(currentRawFile[,2])

  #Parse retention times
  RTSplit<-data.frame(strsplit(currentRawFile[,2], " , "), stringsAsFactors = F)
  RTSplit[1,]<-gsub("\"", "", RTSplit[1,])
  RTSplit[2,]<-gsub("\"", "", RTSplit[2,])
  currentRawFile[,"RT1"]<-as.numeric(t(RTSplit[1,]))
  currentRawFile[,"RT2"]<-as.numeric(t(RTSplit[2,]))

  #Remove duplicate metabolites
  uniqueIndex<-data.frame(paste(currentRawFile[,1], currentRawFile[,2], currentRawFile[,3]))
  currentRawFile<-currentRawFile[which(!duplicated(uniqueIndex)),]
  row.names(currentRawFile)<-c(1:nrow(currentRawFile))

  #Parse ion spectra column into list of numeric vectors
  currentRawFileSplit<-split(currentRawFile,1:nrow(currentRawFile))
  spectraSplit<-lapply(currentRawFileSplit, function(a) strsplit(a[[4]]," "))
  spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
  spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
  spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),])
  spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))

  #Identify absent ions
  spectraFrame<-do.call(cbind,spectraSplit)
  spectraFrame<-spectraFrame[,-c(seq(1,ncol(spectraFrame),2))]
  spectraSums<-colSums(spectraFrame)
  ionProportions<-apply(spectraFrame,1,function(x) sum((x/spectraSums)>absentIonThreshold))
  AbsentIons<-spectraSplit[[1]][,1][which(ionProportions==0)]

  #Calculate pairwise retention time comparisons
  RT1Index<-matrix(unlist(lapply(currentRawFile[,"RT1"],function(x) abs(x-currentRawFile[,"RT1"]))),nrow=length(spectraSplit))
  RT2Index<-matrix(unlist(lapply(currentRawFile[,"RT2"],function(x) abs(x-currentRawFile[,"RT2"])*100)),nrow=length(spectraSplit))

  #Calculate change global metabolite spectral similarity after dropping each ion not previously identified as an absent ion
  for(Ion in possibleIons[!possibleIons%in%AbsentIons]){
    CommonIons<-c(Ion,AbsentIons)
    spectraSplitMask<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%CommonIons),])
    spectraSplitMask<-lapply(spectraSplitMask, function(d) d[order(d[,1]),])
    spectraSplitMask<-lapply(spectraSplitMask, function(d) apply(d,2,as.numeric))
    SimilarityScores<- parallel::mclapply(spectraSplitMask, function(e) lapply(spectraSplitMask, function(d) ((e[,2]%*%d[,2])/(sqrt(sum(e[,2]*e[,2]))*sqrt(sum(d[,2]*d[,2]))))*100), mc.cores=numCores)
    SimilarityMatrix<-matrix(unlist(SimilarityScores), nrow=length(spectraSplit))
    SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
    SimilarityMatrix[lower.tri(SimilarityMatrix, diag = T)]<-0
    Sum50[Ion]<-sum(rowSums(SimilarityMatrix>50))
  }
  Sum50<-Sum50[which(!is.na(Sum50))]
  names(Sum50)<-possibleIons[!possibleIons%in%AbsentIons]
  ProblemIons<-possibleIons[!possibleIons%in%AbsentIons][which(scale(Sum50)<(-1*commonIonThreshold))]

  #Plot commmon ion results
  if(plotData==T){
    graphics::plot(possibleIons[!possibleIons%in%AbsentIons],scale(Sum50), pch=16,cex=0.5, ylab="Sum 50 (Std. Dev.)", xlab="Ion")
    if(length(ProblemIons)>0){
      text(ProblemIons,scale(Sum50)[which(scale(Sum50)<(-2))], labels = ProblemIons, pos=3, cex=0.75, col="red")
    }
  }

  #Combine absent ions and common ions into a dataframe to output
  ionsToRemove<-possibleIons[possibleIons%in%AbsentIons|possibleIons%in%ProblemIons]
  resultFrame<-data.frame("Ions"=possibleIons,Status=rep(NA,length(possibleIons)))
  resultFrame[which(resultFrame[,1]%in%AbsentIons),2]<-"Absent"
  resultFrame[which(resultFrame[,1]%in%ProblemIons),2]<-"Common"
  return(resultFrame[which(!is.na(resultFrame[,2])),])
}

