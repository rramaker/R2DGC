#' This function is an optional pre-prossessing step before running consensus align to identify peaks that likely need to be combined prior to running consensus align and will perform a rough combine of these peaks depending on the quant method as an output.

#' @param inputFileList A character vector with full file paths to chromatof files for processing.
#' @param RT1Penalty A numeric indicating penalty used for first retention time differences. Defaults to 1
#' @param RT2Penalty A numeric indicating penalty used for second retention time differences. Defaults to 100
#' @param similarityCutoff A numeric indicating the similarity threshold (max=100) to use for declaring peaks to combine. Defaults to 95
#' @param numCores Number of cores used to parallelize alignment. See parallel package for details. Defaults to 1
#' @param commonIons A numeric vector of ions to exclude from alignment scores. Can provide first column of output from FindProblemIons function.
#' @param quantMethod Character indicating the quant method used in computing peak areas on chromatof. Accepts "U", "T", or "A" for unique mass, total ion chromatograph or apexing mass. Defaults to "T". If "T" or "A", peaks meeting similarity thresholds will simply be summed. If "U", peaks with the same unique mass with be summed and a proportional conversion will be used before combining peaks with different unique masses.
#' @param outputFiles A boolean indicating if putative peak combinations should be outputted. Will be present at the same path as the input file with _Processed.txt appended to the end.

#' @import parallel
#' @import utils
#' @export
#' @return Returns a data frame with peaks recommended to be combined. If outputFiles is TRUE, peaks returned will be combined and new sample files will be written to the original directory with "_Processed.txt" added to the file name.
#' @examples
#' PrecompressFiles(inputFileList=system.file("extdata", "SampleA.txt", package="R2DGC"))

PrecompressFiles<-function(inputFileList, RT1Penalty=1, RT2Penalty=10,similarityCutoff=95, numCores=1, commonIons=c(), quantMethod="T", outputFiles=F){

  #Create empty list to store processed dataframes for output
  processedFileList<-list()

  #Create empty list to store record of all peaks that should be combined
  combinedList<-list()


  for(File in inputFileList){

    #Read in file
    currentRawFile<-read.table(File, sep="\t", fill=T, quote="",strip.white = T, stringsAsFactors = F,header=T)
    currentRawFile[,4]<-as.character(currentRawFile[,4])
    currentRawFile<-currentRawFile[which(!is.na(currentRawFile[,3])&nchar(currentRawFile[,4])!=0),]
    currentRawFile[,2]<-as.character(currentRawFile[,2])

    #Parse retention times
    RTSplit<-data.frame(strsplit(currentRawFile[,2], " , "), stringsAsFactors = F)
    RTSplit[1,]<-gsub("\"", "", RTSplit[1,])
    RTSplit[2,]<-gsub("\"", "", RTSplit[2,])
    currentRawFile[,"RT1"]<-as.numeric(t(RTSplit[1,]))
    currentRawFile[,"RT2"]<-as.numeric(t(RTSplit[2,]))

    #Remove identical metabolite rows
    uniqueIndex<-data.frame(paste(currentRawFile[,1], currentRawFile[,2], currentRawFile[,3]))
    currentRawFile<-currentRawFile[which(!duplicated(uniqueIndex)),]
    row.names(currentRawFile)<-c(1:nrow(currentRawFile))

    #Parse metabolite spectra into a list
    currentRawFileSplit<-split(currentRawFile,1:nrow(currentRawFile))
    spectraSplit<-lapply(currentRawFileSplit, function(a) strsplit(a[[4]]," "))
    spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),])
    spectraSplit<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
    spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))

    #Calculate pair wise similarity scores between all metabolite spectras
    SimilarityScores<- parallel::mclapply(spectraSplit, function(e) lapply(spectraSplit, function(d) ((e[,2]%*%d[,2])/(sqrt(sum(e[,2]*e[,2]))*sqrt(sum(d[,2]*d[,2]))))*100), mc.cores=numCores)
    SimilarityMatrix<-matrix(unlist(SimilarityScores), nrow=length(SimilarityScores))

    #Subtract retention time difference penalties from similarity scores
    RT1Index<-matrix(unlist(lapply(currentRawFile[,"RT1"],function(x) abs(x-currentRawFile[,"RT1"])*RT1Penalty)),nrow=length(SimilarityScores))
    RT2Index<-matrix(unlist(lapply(currentRawFile[,"RT2"],function(x) abs(x-currentRawFile[,"RT2"])*RT2Penalty)),nrow=length(SimilarityScores))
    SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
    diag(SimilarityMatrix)<-0

    #Find metabolites to with similarity scores greater than similarityCutoff to combine
    MatchList<-apply(SimilarityMatrix,1,function(x) which(x>=similarityCutoff))

    #Initialize number of times to loop through match list in case more than two peaks need to be combined
    NumReps<-0
    if(length(MatchList)>0){
      NumReps<-max(unlist(lapply(MatchList, length)))-1

      if(quantMethod=="U"){

        #Find mates to combine
        Mates<-lapply(MatchList,function(x) x[1])
        BindingQMs<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),5]
        BindingAreas<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),3]
        BindingSpectra<-spectraSplit[unlist(Mates[which(!is.na(Mates))])]

        #Find mate partner to combine
        toBind<-currentRawFile[which(!is.na(Mates)),]
        #Add peak info to combinedList to for output
        combinedList[[File]]<-cbind(toBind,currentRawFile[unlist(Mates[which(!is.na(Mates))]),],File)
        toBind[,"Bound"]<-rep(NA, nrow(toBind))
        toBindQMs<-toBind[,5]
        toBindSpectra<-spectraSplit[which(!is.na(Mates))]

        #Perform proportional conversion to adjust peak areas with differing unique masses
        ConvNumerator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(toBindSpectra[[x]][,1]==toBindQMs[x]),2]))
        ConvDenominator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(toBindSpectra[[x]][,1]==BindingQMs[x]),2]))
        ConvDenominator[which(ConvDenominator==0)]<-NA
        toBind[,3]<-(toBind[,3]*(ConvNumerator/ConvDenominator))+BindingAreas

        #Make sure only one combination (mate to partner) is included in output dataframe
        toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
        toBind<-toBind[which(!duplicated(toBind$Bound)),]

        #Modify sample metabolite frame to include only combined peak
        currentRawFile<-currentRawFile[which(is.na(Mates)),]
        currentRawFile<-rbind(currentRawFile,toBind[,-ncol(toBind)])
      }

      if(quantMethod=="A"|quantMethod=="T"){

        #Find mates to combine
        Mates<-lapply(MatchList,function(x) x[1])
        BindingAreas<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),3]

        #Find mates partners to combine
        toBind<-currentRawFile[which(!is.na(Mates)),]
        #Add peak info to combined list for output
        combinedList[[File]]<-cbind(toBind,currentRawFile[unlist(Mates[which(!is.na(Mates))]),],File)
        toBind[,"Bound"]<-rep(NA, nrow(toBind))

        #Sum peak areas
        toBind[,3]<-toBind[,3]+BindingAreas

        #Ensure only one peak combination gets included in output
        toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
        toBind<-toBind[which(!duplicated(toBind$Bound)),]

        #Update sample metabolite file to include on combined peak
        currentRawFile<-currentRawFile[which(is.na(Mates)),]
        currentRawFile<-rbind(currentRawFile,toBind[,-ncol(toBind)])
      }
    }

    #If any metabolites had greater than two peaks to combine, loop through and make those combinations iteratively
    if(NumReps>0){
      for(Rep in 1:NumReps){

        #Repeat similarity scores with combined peaks
        row.names(currentRawFile)<-c(1:nrow(currentRawFile))
        currentRawFileSplit<-split(currentRawFile,1:nrow(currentRawFile))
        spectraSplit<-lapply(currentRawFileSplit,function(a) strsplit(a[[4]]," "))
        spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
        spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
        spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),])
        spectraSplit<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
        spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))
        SimilarityScores<- parallel::mclapply(spectraSplit, function(e) lapply(spectraSplit, function(d) ((e[,2]%*%d[,2])/(sqrt(sum(e[,2]*e[,2]))*sqrt(sum(d[,2]*d[,2]))))*100), mc.cores=numCores)
        SimilarityMatrix<-matrix(unlist(SimilarityScores), nrow=length(SimilarityScores))
        RT1Index<-matrix(unlist(lapply(currentRawFile[,"RT1"],function(x) abs(x-currentRawFile[,"RT1"])*RT1Penalty)),nrow=length(SimilarityScores))
        RT2Index<-matrix(unlist(lapply(currentRawFile[,"RT2"],function(x) abs(x-currentRawFile[,"RT2"])*RT2Penalty)),nrow=length(SimilarityScores))
        SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
        diag(SimilarityMatrix)<-0

        #Repeat peak combination if more combinations are necessary
        MatchList<-apply(SimilarityMatrix,1,function(x) which(x>=similarityCutoff))
        if(length(MatchList)>0){
          if(quantMethod=="U"){
            Mates<-lapply(MatchList,function(x) x[1])
            BindingQMs<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),5]
            BindingAreas<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),3]
            BindingSpectra<-spectraSplit[unlist(Mates[which(!is.na(Mates))])]
            toBind<-currentRawFile[which(!is.na(Mates)),]
            combinedList[[File]]<-rbind(combinedList[[File]],cbind(toBind,currentRawFile[unlist(Mates[which(!is.na(Mates))]),],File))
            toBind[,"Bound"]<-rep(NA, nrow(toBind))
            toBindQMs<-toBind[,5]
            toBindSpectra<-spectraSplit[which(!is.na(Mates))]
            ConvNumerator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(toBindSpectra[[x]][,1]==toBindQMs[x]),2]))
            ConvDenominator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(toBindSpectra[[x]][,1]==BindingQMs[x]),2]))
            ConvDenominator[which(ConvDenominator==0)]<-NA
            toBind[,3]<-toBind[,3]*(ConvNumerator/ConvDenominator)
            toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
            toBind<-toBind[which(!duplicated(toBind$Bound)),]
            currentRawFile<-currentRawFile[which(is.na(Mates)),]
            currentRawFile<-rbind(currentRawFile,toBind[,-ncol(toBind)])
          }
          if(quantMethod=="A"|quantMethod=="T"){
            Mates<-lapply(MatchList,function(x) x[1])
            BindingAreas<-currentRawFile[unlist(Mates[which(!is.na(Mates))]),3]
            toBind<-currentRawFile[which(!is.na(Mates)),]
            combinedList[[File]]<-rbind(combinedList[[File]],cbind(toBind,currentRawFile[unlist(Mates[which(!is.na(Mates))]),],File))
            toBind[,"Bound"]<-rep(NA, nrow(toBind))
            toBind[,3]<-toBind[,3]+BindingAreas
            toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
            toBind<-toBind[which(!duplicated(toBind$Bound)),]
            currentRawFile<-currentRawFile[which(is.na(Mates)),]
            currentRawFile<-rbind(currentRawFile,toBind[,-ncol(toBind)])
          }
        }
      }
    }

    #Add processed file with putative peak area combinations to processedFileList for output
    processedFileList[[File]]<-currentRawFile[,1:5]
  }

  #Make data frame with all combined peak pair info
  combinedFrame<-do.call(rbind,combinedList)
  if(length(combinedFrame)>0){
    row.names(combinedFrame)<-1:nrow(combinedFrame)
  }
  #If outputFiles==TRUE, write processed files out to the input file directory
  if(outputFiles==TRUE){
    for(File in inputFileList){
      utils::write.table(processedFileList[[File]], paste0(substr(File,1,nchar(File)-4),"_Processed.txt"),sep="\t",quote=F,row.names=F)
    }
  }
  return(combinedFrame)
}
