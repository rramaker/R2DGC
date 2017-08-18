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

  #Create empty list to store record of all peaks that should be combined
  combinedList<-list()

  ImportFile<-function(File){

    #Read in file
    currentRawFile<-read.table(File, sep="\t", fill=T, quote="",strip.white = T, stringsAsFactors = F,header=T)
    currentRawFile[,5]<-as.character(currentRawFile[,5])
    currentRawFile<-currentRawFile[which(!is.na(currentRawFile[,3])&nchar(currentRawFile[,5])!=0),]
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
    spectraSplit<-lapply(currentRawFileSplit, function(a) strsplit(a[[5]]," "))
    spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    spectraSplit<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
    spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))
    ionNames<-spectraSplit[[1]][order(spectraSplit[[1]][,1]),1]
    spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),2,drop=F])
    return(list(currentRawFile,spectraSplit, ionNames))
  }
  importedFiles<-mclapply(inputFileList, ImportFile, mc.cores=numCores)

  #Calculate pair wise similarity scores between all metabolite spectras
  FindMatches<-function(Sample){
    spectraFrame<-do.call(cbind,Sample[[2]])
    spectraFrame<-t(spectraFrame)
    spectraFrame<-as.matrix(spectraFrame)/sqrt(apply((as.matrix(spectraFrame))^2,1,sum))
    SimilarityMatrix<-(spectraFrame %*% t(spectraFrame))*100

    #Subtract retention time difference penalties from similarity scores
    RT1Index<-matrix(unlist(lapply(Sample[[1]][,"RT1"],function(x) abs(x-Sample[[1]][,"RT1"])*RT1Penalty)),nrow=nrow(SimilarityMatrix))
    RT2Index<-matrix(unlist(lapply(Sample[[1]][,"RT2"],function(x) abs(x-Sample[[1]][,"RT2"])*RT2Penalty)),nrow=nrow(SimilarityMatrix))
    SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
    diag(SimilarityMatrix)<-0

    #Find metabolites to with similarity scores greater than similarityCutoff to combine
    return(apply(SimilarityMatrix,1,function(x) which(x>=similarityCutoff)))
  }
  MatchList<-mclapply(importedFiles, FindMatches, mc.cores=numCores)

  #Initialize number of times to loop through match list in case more than two peaks need to be combined
  for(SampNum in 1:length(MatchList)){
    NumReps<-0
    if(length(MatchList[[SampNum]])>0){
      NumReps<-max(unlist(lapply(MatchList[[SampNum]], length)))-1

      if(quantMethod=="U"){

        #Find mates to combine
        Mates<-lapply(MatchList[[SampNum]],function(x) x[1])
        BindingQMs<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),4]
        BindingAreas<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),3]
        BindingSpectra<-importedFiles[[SampNum]][[2]][unlist(Mates[which(!is.na(Mates))])]
        ionNames<-importedFiles[[SampNum]][[3]]

        #Find mate partner to combine
        toBind<-importedFiles[[SampNum]][[1]][which(!is.na(Mates)),]
        #Add peak info to combinedList to for output
        combinedList[[inputFileList[SampNum]]]<-cbind(toBind,importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),],inputFileList[SampNum])
        toBind[,"Bound"]<-rep(NA, nrow(toBind))
        toBindQMs<-toBind[,4]
        toBindSpectra<-importedFiles[[SampNum]][[2]][which(!is.na(Mates))]

        #Perform proportional conversion to adjust peak areas with differing unique masses
        ConvNumerator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(ionNames==toBindQMs[x])]))
        ConvDenominator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(ionNames==BindingQMs[x])]))
        ConvDenominator[which(ConvDenominator==0)]<-NA
        toBind[,3]<-(toBind[,3]*(ConvNumerator/ConvDenominator))+BindingAreas

        #Make sure only one combination (mate to partner) is included in output dataframe
        toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
        toBind<-toBind[which(!duplicated(toBind$Bound)),]

        #Modify sample metabolite frame to include only combined peak
        importedFiles[[SampNum]][[1]]<-importedFiles[[SampNum]][[1]][which(is.na(Mates)),]
        importedFiles[[SampNum]][[1]]<-rbind(importedFiles[[SampNum]][[1]],toBind[,-ncol(toBind)])
      }

      if(quantMethod=="A"|quantMethod=="T"){

        #Find mates to combine
        Mates<-lapply(MatchList[[SampNum]],function(x) x[1])
        BindingAreas<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),3]

        #Find mates partners to combine
        toBind<-importedFiles[[SampNum]][[1]][which(!is.na(Mates)),]
        #Add peak info to combined list for output
        combinedList[[inputFileList[SampNum]]]<-cbind(toBind,importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),],inputFileList[SampNum])
        toBind[,"Bound"]<-rep(NA, nrow(toBind))

        #Sum peak areas
        toBind[,3]<-toBind[,3]+BindingAreas

        #Ensure only one peak combination gets included in output
        toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
        toBind<-toBind[which(!duplicated(toBind$Bound)),]

        #Update sample metabolite file to include on combined peak
        importedFiles[[SampNum]][[1]]<-importedFiles[[SampNum]][[1]][which(is.na(Mates)),]
        importedFiles[[SampNum]][[1]]<-rbind(importedFiles[[SampNum]][[1]],toBind[,-ncol(toBind)])
      }
    }


    #If any metabolites had greater than two peaks to combine, loop through and make those combinations iteratively
    if(NumReps>0){
      for(Rep in 1:NumReps){

        #Repeat similarity scores with combined peaks
        row.names(importedFiles[[SampNum]][[1]])<-c(1:nrow(importedFiles[[SampNum]][[1]]))
        currentRawFileSplit<-split(importedFiles[[SampNum]][[1]],1:nrow(importedFiles[[SampNum]][[1]]))
        spectraSplit<-lapply(currentRawFileSplit,function(a) strsplit(a[[5]]," "))
        spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
        spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
        spectraSplit<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
        spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))
        spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),2,drop=F])
        spectraFrame<-do.call(cbind,spectraSplit)
        spectraFrame<-t(spectraFrame)
        spectraFrame<-as.matrix(spectraFrame)/sqrt(apply((as.matrix(spectraFrame))^2,1,sum))
        SimilarityMatrix<-(spectraFrame %*% t(spectraFrame))*100

        #Subtract retention time difference penalties from similarity scores
        RT1Index<-matrix(unlist(lapply(importedFiles[[SampNum]][[1]][,"RT1"],function(x) abs(x-importedFiles[[SampNum]][[1]][,"RT1"])*RT1Penalty)),nrow=nrow(SimilarityMatrix))
        RT2Index<-matrix(unlist(lapply(importedFiles[[SampNum]][[1]][,"RT2"],function(x) abs(x-importedFiles[[SampNum]][[1]][,"RT2"])*RT2Penalty)),nrow=nrow(SimilarityMatrix))
        SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
        diag(SimilarityMatrix)<-0

        #Repeat peak combination if more combinations are necessary
        NewMatchList<-apply(SimilarityMatrix,1,function(x) which(x>=similarityCutoff))
        if(length(NewMatchList)>0){
          if(quantMethod=="U"){
            Mates<-lapply(NewMatchList,function(x) x[1])
            BindingQMs<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),4]
            BindingAreas<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),3]
            BindingSpectra<-spectraSplit[unlist(Mates[which(!is.na(Mates))])]
            toBind<-importedFiles[[SampNum]][[1]][which(!is.na(Mates)),]
            combinedList[[inputFileList[SampNum]]]<-rbind(combinedList[[inputFileList[SampNum]]],cbind(toBind,importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),],inputFileList[SampNum]))
            toBind[,"Bound"]<-rep(NA, nrow(toBind))
            toBindQMs<-toBind[,4]
            toBindSpectra<-spectraSplit[which(!is.na(Mates))]
            ConvNumerator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(ionNames==toBindQMs[x])]))
            ConvDenominator<-unlist(lapply(1:length(toBindQMs), function(x) toBindSpectra[[x]][which(ionNames==BindingQMs[x])]))
            ConvDenominator[which(ConvDenominator==0)]<-NA
            toBind[,3]<-toBind[,3]*(ConvNumerator/ConvDenominator)
            toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
            toBind<-toBind[which(!duplicated(toBind$Bound)),]
            importedFiles[[SampNum]][[1]]<-importedFiles[[SampNum]][[1]][which(is.na(Mates)),]
            importedFiles[[SampNum]][[1]]<-rbind(importedFiles[[SampNum]][[1]],toBind[,-ncol(toBind)])
          }
          if(quantMethod=="A"|quantMethod=="T"){
            Mates<-lapply(NewMatchList,function(x) x[1])
            BindingAreas<-importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),3]
            toBind<-importedFiles[[SampNum]][[1]][which(!is.na(Mates)),]
            combinedList[[inputFileList[SampNum]]]<-rbind(combinedList[[inputFileList[SampNum]]],cbind(toBind,importedFiles[[SampNum]][[1]][unlist(Mates[which(!is.na(Mates))]),],inputFileList[SampNum]))
            toBind[,"Bound"]<-rep(NA, nrow(toBind))
            toBind[,3]<-toBind[,3]+BindingAreas
            toBind$Bound<-paste(toBind$Bound,apply(cbind(unlist(Mates[which(!is.na(Mates))]),which(!is.na(Mates))),1,min),sep="_")
            toBind<-toBind[which(!duplicated(toBind$Bound)),]
            importedFiles[[SampNum]][[1]]<-importedFiles[[SampNum]][[1]][which(is.na(Mates)),]
            importedFiles[[SampNum]][[1]]<-rbind(importedFiles[[SampNum]][[1]],toBind[,-ncol(toBind)])
          }
        }
      }
    }
  }

  #Make data frame with all combined peak pair info
  combinedFrame<-do.call(rbind,combinedList)
  if(length(combinedFrame)>0){
    row.names(combinedFrame)<-1:nrow(combinedFrame)
  }
  #If outputFiles==TRUE, write processed files out to the input file directory
  if(outputFiles==TRUE){
    for(SampNum in 1:length(importedFiles)){
      utils::write.table(importedFiles[[SampNum]][[1]][,1:5], paste0(substr(inputFileList[[SampNum]],1,nchar(inputFileList[[SampNum]])-4),"_Processed.txt"),sep="\t",quote=F,row.names=F)
    }
  }
  return(combinedFrame)
}
