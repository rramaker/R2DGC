#' This function takes input chromatof files from metabolite standards and parses them into a dataframe of retention time indexed standards that can be used as an input for ConsensusAlign function
#'
#' @param inputFileList A character vector of full file paths to the metabolite standard chromatof files to include in library.
#' @param RT1_Standards A character vector with the name of all first retention time standards to use to index metabolites. Defaults to NULL.
#' @param RT2_Standards A character vector with the name of all second retention time standards to use to index metabolites. Defaults to NULL.
#' @return Returns a list with a retention standard indexed metabolite library that can be used in the standard library argument of the ConsensusAlign function
#' @examples
#' MakeReference(c(system.file("extdata", "Alanine_150226_1.txt", package="R2DGC"),
#'     system.file("extdata", "Serine_022715_1.txt", package="R2DGC")))
#' @export

MakeReference<-function(inputFileList, RT1_Standards=NULL, RT2_Standards=NULL){

  #Create empty list to add RT-indexed standards
  StandardLibList<-list()

  #Create empty list for missing retention indices
  MissingRTIndices<-list()

  for(File in inputFileList){

    #Read in file
    currentFile<-read.table(File, sep="\t", fill=T, quote="",strip.white = T, stringsAsFactors = F, header=T)

    #Parse retention time
    RTSplit<-data.frame(strsplit(currentFile[,2], " , "), stringsAsFactors = F)
    RTSplit[1,]<-gsub("\"", "", RTSplit[1,])
    RTSplit[2,]<-gsub("\"", "", RTSplit[2,])
    currentFile[,"RT1"]<-as.numeric(t(RTSplit[1,]))
    currentFile[,"RT2"]<-as.numeric(t(RTSplit[2,]))

    #Find RT1 differences from all RT1 standards
    if(!is.null(RT1_Standards)){
      if(sum(RT1_Standards%in%currentFile[,1])!=length(RT1_Standards)){
        MissingRTIndices[[File]]<-RT1_Standards[which(!RT1_Standards%in%currentFile[,1])]
        next
      }
      RT1_Length<-max(currentFile[which(currentFile[,1]%in%RT1_Standards),4])-min(currentFile[which(currentFile[,1]%in%RT1_Standards),4])
      for(Standard in RT1_Standards){
        currentFile[,paste(Standard,"RT1",sep="_")]<-(currentFile[,4]-currentFile[grep(Standard,currentFile[,1],perl = T),4])/RT1_Length
      }
    }

    #Find RT2 differences from all RT2 standards
    if(!is.null(RT2_Standards)){
      if(sum(RT2_Standards%in%currentFile[,1])!=length(RT2_Standards)){
        MissingRTIndices[[File]]<-c(MissingRTIndices[[File]],RT2_Standards[which(!RT1_Standards%in%currentFile[,1])])
        next
      }
      RT2_Length<-max(currentFile[which(currentFile[,1]%in%RT2_Standards),4])-min(currentFile[which(currentFile[,1]%in%RT2_Standards),4])
      for(Standard in RT2_Standards){
        currentFile[,paste(Standard,"RT2",sep="_")]<-(currentFile[,4]-currentFile[grep(Standard,currentFile[,1],perl = T),4])/RT2_Length
      }
    }

    #Add indexed metabolite only to StandardLibList
    StandardLibList[[File]]<-currentFile[which(!currentFile[,1]%in%RT1_Standards&!currentFile[,1]%in%RT2_Standards),]
  }

  #Convert library to dataframe for output
  if(length(MissingRTIndices)==0){
    StandardLibrary<-do.call(rbind,StandardLibList)
    row.names(StandardLibrary)<-1:nrow(StandardLibrary)
    return(StandardLibrary)
  }
  if(length(MissingRTIndices)>0){
    message("Error: Missing RT indices detected. See output list")
   return(MissingRTIndices)
  }
}

