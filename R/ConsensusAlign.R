#' Takes a vector of paths to input files and aligns common metabolites into a final table.  Will also identify metabolites if a reference library is provided

#' @param inputFileList Vector of file paths to align
#' @param RT1_Standards Vector of standard names used to adjust first retention time. All names must be found in input files. Defaults to NULL.
#' @param RT2_Standards Vector of standard names used to adjust second retention time. All names must be found in input files. Defaults to NULL.
#' @param seedFile File number in inputFileList to initialize alignment. Can also input a vector of different seed files (3 is usually sufficient) to prevent bias from seed file.  Defaults to 1.
#' @param RT1Penalty Penalty used for first retention time errors.  Defaults to 1.
#' @param RT2Penalty Penalty used for first retention time errors.  Defaults to 10.
#' @param autoTuneMatchStringency Will automatically find optimal match threshold. If TRUE, will ignore similarityCutoff. Defaults to TRUE.
#' @param similarityCutoff Adjusts peak similarity threshold required for alignment. Adjust in concordance with RT1 and RT2 penalties. Will be ignored if autoTuneMatchStrigency is TRUE. Defaults to 90.
#' @param disimilarityCutoff Defaults to similarityCutoff-90. Sets the threshold for including a new peak in the alignment table to ensure new metabolites aren't just below alignment thresholds
#' @param numCores Number of cores used to parallelize alignment. See parallel package. Defaults to 4.
#' @param commonIons Provide a vector of ions to ignore from the FindProblemIons function. Defaults to empty vector.
#' @param missingValueLimit Maximum fraction (Numeric between 0 and 1) of missing values acceptable for retaining a metabolite in the final alignment table. Defaults to 0.25.
#' @param missingPeakFinderSimilarityLax Fraction of Similarity Cutoff to use to find missing alignments just below threshold. Set to 1 to prevent searching for missing peaks. Defaults to 0.85.
#' @param quantMethod Set to U, A, or T to indicate if unique mass (U), appexing masses (A), or total ion chormatograph (T) was used to quantify peak areas. Defaults to T.  If "T" or "A", peaks meeting similarity thresholds will simply be summed. If "U", peaks with the same unique mass with be summed and a proportional conversion will be used before combining peaks with different unique masses.
#' @param standardLibrary Defaults to NULL. Provide standard library generated from MakeReference function to ID metabolites with retention index.

#' @return A list with three items: AlignmentMatix - A dataframe with peak areas for all metabolites matched in sufficient number of samples. MetaboliteInfo - An info file with RT, spectra, and metabolite ID info for each metabolite in the AlignmentMatrix. UnmatchedQuantMasses- Info on metabolites combined that had different unique masses (if quantMethod="U") or greater than 50% different apexing masses (if quantMethod="A")
#' @import parallel
#' @import stats
#' @export
#' @examples
#' ConsensusAlign(c(system.file("extdata", "SampleA.txt", package="R2DGC"),
#'     system.file("extdata", "SampleB.txt", package="R2DGC")), RT1_Standards= c())

ConsensusAlign<-function(inputFileList,
                         RT1_Standards=NULL,
                         RT2_Standards=NULL,
                         seedFile=1, #Change to 5 for test case reproducibility
                         RT1Penalty=1,
                         RT2Penalty=10, #1 for relatively unstable spectras
                         autoTuneMatchStringency=TRUE,
                         similarityCutoff=90,
                         disimilarityCutoff=similarityCutoff-90,
                         numCores=1,
                         commonIons=c(),
                         missingValueLimit=0.75,
                         missingPeakFinderSimilarityLax=0.85,
                         quantMethod="T",
                         standardLibrary=NULL){

  #Function to import files
  ImportFile<-function(File){
    MissingStandards<-c()
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

    #Index metabolites by RT1 Standards
    if(!is.null(RT1_Standards)){
      #Check if all RT1 standards are present in each file
      if(sum(RT1_Standards%in%currentRawFile[,1])!=length(RT1_Standards)){
        MissingStandards<-c(MissingStandards,paste(File,"file missing RT1 standards:",RT1_Standards[which(!RT1_Standards%in%currentRawFile[,1])]," ",sep=" "))
        #break
      }
      #Index each metabolite by RT1 Standards
      RT1_Length<-max(currentRawFile[which(currentRawFile[,1]%in%RT1_Standards),6])-min(currentRawFile[which(currentRawFile[,1]%in%RT1_Standards),6])
      for(Standard in RT1_Standards){
        currentRawFile[,paste(Standard,"RT1",sep="_")]<-(currentRawFile[,6]-currentRawFile[grep(Standard,currentRawFile[,1],perl = T),6])/RT1_Length
      }
    }
    #Index metabolites by RT2 Standards
    if(!is.null(RT2_Standards)){
      #Check if all RT2_Standards are present
      if(sum(RT2_Standards%in%currentRawFile[,1])!=length(RT2_Standards)){
        MissingStandards<-c(MissingStandards,paste(File,"file missing RT2 standards:",RT2_Standards[which(!RT2_Standards%in%currentRawFile[,1])]," ",sep=" "))
        #break
      }
      #Index each metabolite by RT2 standards
      RT2_Length<-max(currentRawFile[which(currentRawFile[,1]%in%RT2_Standards),6])-min(currentRawFile[which(currentRawFile[,1]%in%RT2_Standards),6])
      for(Standard in RT2_Standards){
        currentRawFile[,paste(Standard,"RT2",sep="_")]<-(currentRawFile[,6]-currentRawFile[grep(Standard,currentRawFile[,1],perl = T),6])/RT2_Length
      }
    }
    #Parse metabolite spectra into a list
    currentRawFileSplit<-split(currentRawFile,1:nrow(currentRawFile))
    spectraSplit<-lapply(currentRawFileSplit, function(a) strsplit(a[[5]]," "))
    spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    spectraSplitFilt<-lapply(spectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
    spectraSplitFilt<-lapply(spectraSplitFilt, function(d) apply(d,2,as.numeric))
    spectraSplitFilt<-lapply(spectraSplitFilt, function(d) d[order(d[,1]),2,drop=F])
    spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))
    ionNames<-spectraSplit[[1]][order(spectraSplit[[1]][,1]),1]
    spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),2,drop=F])
    return(list(currentRawFile,spectraSplitFilt, MissingStandards, ionNames, spectraSplit))
  }
  ImportedFiles<-mclapply(inputFileList, ImportFile, mc.cores=numCores)

  MissingFileList<-c()
  for(File in ImportedFiles){
    MissingFileList<-c(MissingFileList,File[3])
  }
  if(length(unlist(MissingFileList))>0){
    stop(unlist(MissingFileList), call.=FALSE)
  }

  #Function to calculate pair wise similarity scores between all metabolite spectras
  GenerateSimFrames<-function(Sample, SeedSample){
    seedSpectraFrame<-do.call(cbind,SeedSample[[2]])
    seedSpectraFrame<-t(seedSpectraFrame)
    seedSpectraFrame<-as.matrix(seedSpectraFrame)/sqrt(apply((as.matrix(seedSpectraFrame))^2,1,sum))
    sampleSpectraFrame<-do.call(cbind,Sample[[2]])
    sampleSpectraFrame<-t(sampleSpectraFrame)
    sampleSpectraFrame<-as.matrix(sampleSpectraFrame)/sqrt(apply((as.matrix(sampleSpectraFrame))^2,1,sum))
    SimilarityMatrix<-(seedSpectraFrame %*% t(sampleSpectraFrame))*100

    #Calculate pairwise RT penalties for each current file metabolite and seed file metabolites
    RT1Index<-matrix(unlist(lapply(Sample[[1]][,"RT1"],function(x) abs(x-SeedSample[[1]][,"RT1"])*RT1Penalty)),nrow=nrow(SimilarityMatrix))
    RT2Index<-matrix(unlist(lapply(Sample[[1]][,"RT2"],function(x) abs(x-SeedSample[[1]][,"RT2"])*RT2Penalty)),nrow=nrow(SimilarityMatrix))

    #Use RT indices to calculate RT penalties if necessary
    if(!is.null(RT1_Standards)){
      #Compute list of metabolite to RT1 standard differences between current file and seed file for each metabolite
      RT1Index<-list()
      RT1_Length<-max(Sample[[1]][which(Sample[[1]][,1]%in%RT1_Standards),6])-min(Sample[[1]][which(Sample[[1]][,1]%in%RT1_Standards),6])
      for(Standard in RT1_Standards){
        RT1Index[[Standard]]<-matrix(unlist(lapply(Sample[[1]][,paste(Standard,"RT1",sep="_")],function(x) abs(x-SeedSample[[1]][,paste(Standard,"RT1",sep="_")])*(RT1Penalty/length(RT1_Standards)))),nrow=nrow(SimilarityMatrix))*RT1_Length
      }
      #Sum all relative standard differences into a final score
      RT1Index<-Reduce("+",RT1Index)
    }
    if(!is.null(RT2_Standards)){
      #Compute list of metabolite to RT2 standard differences between current file and seed file for each metabolite
      RT2Index<-list()
      RT2_Length<-max(Sample[[1]][which(Sample[[1]][,1]%in%RT2_Standards),6])-min(Sample[[1]][which(Sample[[1]][,1]%in%RT2_Standards),6])
      for(Standard in RT2_Standards){
        RT2Index[[Standard]]<-matrix(unlist(lapply(Sample[[1]][,paste(Standard,"RT2",sep="_")],function(x) abs(x-SeedSample[[1]][,paste(Standard,"RT2",sep="_")])*(RT2Penalty/length(RT2_Standards)))),nrow=nrow(SimilarityMatrix))*RT2_Length
      }
      #Sum all relative standard differences into a final score
      RT2Index<-Reduce("+",RT2Index)
    }
    return(SimilarityMatrix-RT1Index-RT2Index)
  }

  AlignmentTableList<-list()

  for(seed in seedFile){
    message(paste0("seed is ",seed))
    SeedSample<- ImportedFiles[[seed]]

    #Calculate pairwise similarity scores
    SimCutoffs<-mclapply(ImportedFiles[-seed], function(x) GenerateSimFrames(x, SeedSample), mc.cores=numCores)
    names(SimCutoffs)<-as.character(1:length(ImportedFiles))[-seed]

    #Calculate optimal similarity score cutoff if desired
    if(autoTuneMatchStringency==TRUE){
      SimScores<-mclapply(SimCutoffs, function(y) unlist(lapply(1:100,function(x) sum(rowSums(y>x)>0)/(sum(y>x)^0.5))),mc.cores = numCores)
      SimScores<-matrix(unlist(SimScores),ncol=length(SimScores))
      similarityCutoff<-which.max(rowSums(SimScores))
      disimilarityCutoff<-similarityCutoff-90
    }

    #Find Metabolites to add to seed file
    for(SampNum in (1:length(ImportedFiles))[-seed]){
      #Find best matches and mate pairs for each metabolite and remove inferior matches if metabolite is matched twice
      MatchScores<-apply(SimCutoffs[[as.character(SampNum)]],2,function(x) max(x,na.rm=T))
      Mates<-apply(SimCutoffs[[as.character(SampNum)]],2,function(x) which.max(x))
      names(MatchScores)<-1:length(MatchScores)
      names(Mates)<-1:length(Mates)
      Mates<-Mates[order(-MatchScores)]
      MatchScores<-MatchScores[order(-MatchScores)]
      MatchScores[which(duplicated(Mates))]<-NA
      Mates<-Mates[order(as.numeric(names(Mates)))]
      MatchScores<-MatchScores[order(as.numeric(names(MatchScores)))]

      #Find metabolites in current file sufficiently dissimilar to add to alignment matrix
      SeedSample[[1]]<-rbind(SeedSample[[1]],ImportedFiles[[SampNum]][[1]][which(MatchScores<disimilarityCutoff),])
      if(length(which(MatchScores<disimilarityCutoff))>0){
        SeedSample[[2]][as.character((length(SeedSample[[2]])+1):(length(SeedSample[[2]])+length(which(MatchScores<disimilarityCutoff))))]<- ImportedFiles[[SampNum]][[2]][which(MatchScores<disimilarityCutoff)]
      }
    }

    #Repeat pairwise similarity score calculation on full seed sample file
    SimCutoffs<-mclapply(ImportedFiles, function(x) GenerateSimFrames(x, SeedSample), mc.cores=numCores)

    #Calculate optimal similarity score cutoff if desired
    if(autoTuneMatchStringency==TRUE){
      message("Computing peak similarity threshold")
      SimScores<-mclapply(SimCutoffs[-seed], function(y) unlist(lapply(1:100,function(x) sum(rowSums(y>x)>0)/(sum(y>x)^0.5))),mc.cores = numCores)
      SimScores<-matrix(unlist(SimScores),ncol=length(SimScores))
      similarityCutoff<-which.max(rowSums(SimScores))
      disimilarityCutoff<-similarityCutoff-90
    }

    #Establish alignment matrix
    FinalMatrix<-matrix(nrow=nrow(SeedSample[[1]]),ncol=length(inputFileList))
    row.names(FinalMatrix)<-SeedSample[[1]][,1]
    colnames(FinalMatrix)<-inputFileList

    #Establish emptly list to store incongruent quant matches if using quantMethod "U" or "A"
    MissingQMList<-list()
    #Loop back through input files and find matches above similarityCutoff threshold
    for(SampNum in 1:length(ImportedFiles)){
      ionNames<-ImportedFiles[[SampNum]][[4]]
      #Find best seed sample match scores
      MatchScores<-apply(SimCutoffs[[SampNum]],2,function(x) max(x,na.rm=T))
      #Find current sample metabolites with best match
      Mates<-apply(SimCutoffs[[SampNum]],2,function(x) which.max(x))
      names(MatchScores)<-1:length(MatchScores)
      names(Mates)<-1:length(Mates)
      Mates<-Mates[order(-MatchScores)]
      MatchScores<-MatchScores[order(-MatchScores)]
      MatchScores[which(duplicated(Mates))]<-NA
      Mates<-Mates[order(as.numeric(names(Mates)))]
      MatchScores<-MatchScores[order(as.numeric(names(MatchScores)))]

      if(quantMethod=="U"){
        #Find quant masses for each match pair
        MatchedSeedQMs<- SeedSample[[1]][,4][Mates[which(MatchScores>=similarityCutoff)]]
        currentFileQMs<- ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),4]
        #Add incongruent quant mass info to MissingQMList for output
        MissingQMList[[inputFileList[SampNum]]]<-cbind(inputFileList[SampNum],which(MatchScores>=similarityCutoff),currentFileQMs,inputFileList[seed],Mates[which(MatchScores>=similarityCutoff)],MatchedSeedQMs)[which(currentFileQMs!=MatchedSeedQMs),]
        #Convert areas proportionally for incongruent quant masses
        currentFileAreas<- ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),3]
        currentFileSpectra<- ImportedFiles[[SampNum]][[5]][which(MatchScores>=similarityCutoff)]
        MatchedSeedSpectra<- SeedSample[[5]][Mates[which(MatchScores>=similarityCutoff)]]
        ConvNumerator<-unlist(lapply(1:length(currentFileQMs), function(x) currentFileSpectra[[x]][which(ionNames==currentFileQMs[x])]))
        ConvDenominator<-unlist(lapply(1:length(currentFileQMs), function(x) currentFileSpectra[[x]][which(ionNames==MatchedSeedQMs[x])]))
        ConvDenominator[which(ConvDenominator==0)]<-NA
        #Add matched peaks to final alignment matrix
        FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-currentFileAreas*(ConvNumerator/ConvDenominator)
      }

      if(quantMethod=="A"){
        #Make function to parse apexing masses and test whether 50% are in common with seed file
        TestQMOverlap<-function(x){
          SeedQMs<- strsplit(x[1],"\\+")
          FileQMs<- strsplit(x[2],"\\+")
          sum(unlist(SeedQMs)%in%unlist(FileQMs))/min(length(unlist(SeedQMs)),length(unlist(FileQMs)))<0.5
        }
        #Test apexing mass overlap for each metabolite match
        MatchedSeedQMs<- SeedSample[[1]][,4][Mates[which(MatchScores>=similarityCutoff)]]
        currentFileQMs<- ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),4]
        QM_Bind<-cbind(MatchedSeedQMs,currentFileQMs)
        QM_Match<-apply(QM_Bind, 1, function(x) TestQMOverlap(x))
        #Add incongruent apexing masses to MissingQMList for output
        MissingQMList[[inputFileList[SampNum]]]<-cbind(inputFileList[SampNum],which(MatchScores>=similarityCutoff),currentFileQMs,inputFileList[seed],Mates[which(MatchScores>=similarityCutoff)],MatchedSeedQMs)[which(QM_Match==TRUE),]
        #Add matched peaks to final alignment matrix
        FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),3]
      }

      if(quantMethod=="T"){
        #Add newly aligned peaks to alignment matrix
        FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),3]
      }
    }

    #Filter final alignment matrix to only peaks passing the missing value limit
    SeedSample[[1]]<-SeedSample[[1]][which(rowSums(is.na(FinalMatrix))<=round(length(inputFileList)*(1-missingValueLimit))),]
    SeedSample[[2]]<-SeedSample[[2]][which(rowSums(is.na(FinalMatrix))<=round(length(inputFileList)*(1-missingValueLimit)))]
    FinalMatrix<-FinalMatrix[which(rowSums(is.na(FinalMatrix))<=round(length(inputFileList)*(1-missingValueLimit))),]

    #Compute relaxed similarity cutoff with missingPeakFinderSimilarityLax
    similarityCutoff<-similarityCutoff*missingPeakFinderSimilarityLax
    message("Searching for missing peaks")

    #Loop through each file again and check for matches in high probability missing metabolites meeting relaxed similarity cutoff
    SimCutoffs<-mclapply(ImportedFiles, function(x) GenerateSimFrames(x, SeedSample), mc.cores=numCores)

    #Find peaks with missing values
    for(SampNum in 1:length(ImportedFiles)){
      ionNames<-ImportedFiles[[SampNum]][[4]]
      MissingPeaks<-which(is.na(FinalMatrix[,inputFileList[SampNum]]))
      if(length(MissingPeaks)>0){
        MatchScores<-apply(SimCutoffs[[SampNum]],2,function(x) max(x,na.rm=T))
        Mates<-apply(SimCutoffs[[SampNum]],2,function(x) which.max(x))
        names(MatchScores)<-1:length(MatchScores)
        names(Mates)<-1:length(Mates)
        Mates<-Mates[order(-MatchScores)]
        MatchScores<-MatchScores[order(-MatchScores)]
        MatchScores[which(duplicated(Mates))]<-NA
        Mates<-Mates[order(as.numeric(names(Mates)))]
        MatchScores<-MatchScores[order(as.numeric(names(MatchScores)))]
        MatchScores<-MatchScores[which(Mates%in%MissingPeaks)]
        Mates<-Mates[which(Mates%in%MissingPeaks)]

        #If matches are greater than relaxed simlarity cutoff add to final alignment table
        if(length(which(MatchScores>=similarityCutoff))>0){
          if(quantMethod=="U"){
            #Find quant masses for each match pair
            MatchedSeedQMs<- SeedSample[[1]][,4][Mates[which(MatchScores>=similarityCutoff)]]
            currentFileQMs<- ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),4]
            #Add incongruent quant mass info to MissingQMList for output
            MissingQMList[[paste0(inputFileList[SampNum],"_MPF")]]<-cbind(inputFileList[SampNum],which(MatchScores>=similarityCutoff),currentFileQMs,inputFileList[seed],Mates[which(MatchScores>=similarityCutoff)],MatchedSeedQMs)[which(currentFileQMs!=MatchedSeedQMs),]
            #Convert areas proportionally for incongruent quant masses
            currentFileAreas<- ImportedFiles[[SampNum]][[1]][names(which(MatchScores>=similarityCutoff)),3]
            currentFileSpectra<- ImportedFiles[[SampNum]][[5]][names(which(MatchScores>=similarityCutoff))]
            MatchedSeedSpectra<- SeedSample[[5]][Mates[which(MatchScores>=similarityCutoff)]]
            ConvNumerator<-unlist(lapply(1:length(currentFileQMs), function(x) currentFileSpectra[[x]][which(ionNames==currentFileQMs[x])]))
            ConvDenominator<-unlist(lapply(1:length(currentFileQMs), function(x) currentFileSpectra[[x]][which(ionNames==MatchedSeedQMs[x])]))
            ConvDenominator[which(ConvDenominator==0)]<-NA
            #Add matched peaks to final alignment matrix
            FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-currentFileAreas*(ConvNumerator/ConvDenominator)
          }
          if(quantMethod=="A"){
            #Make function to parse apexing masses and test whether 50% are in common with seed file
            TestQMOverlap<-function(x){
              SeedQMs<- strsplit(x[1],"\\+")
              FileQMs<- strsplit(x[2],"\\+")
              return(sum(unlist(SeedQMs)%in%unlist(FileQMs))/min(length(unlist(SeedQMs)),length(unlist(FileQMs)))<0.5)
            }
            #Test apexing mass overlap for each metabolite match
            MatchedSeedQMs<- SeedSample[[1]][,4][Mates[which(MatchScores>=similarityCutoff)]]
            currentFileQMs<- ImportedFiles[[SampNum]][[1]][which(MatchScores>=similarityCutoff),4]
            QM_Bind<-cbind(MatchedSeedQMs,currentFileQMs)
            QM_Match<-apply(QM_Bind, 1, function(x) TestQMOverlap(x))
            #Add incongruent apexing masses to MissingQMList for output
            MissingQMList[[paste0(inputFileList[SampNum],"_MPF")]]<-cbind(inputFileList[SampNum],which(MatchScores>=similarityCutoff),currentFileQMs,inputFileList[seed],Mates[which(MatchScores>=similarityCutoff)],MatchedSeedQMs)[which(currentFileQMs!=MatchedSeedQMs),]
            #Add matched peaks to final alignment matrix
            FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-ImportedFiles[[SampNum]][[1]][names(which(MatchScores>=similarityCutoff)),3]
          }
          if(quantMethod=="T"){
            #Add matched peaks to final alignment matrix
            FinalMatrix[Mates[which(MatchScores>=similarityCutoff)],inputFileList[SampNum]]<-ImportedFiles[[SampNum]][[1]][names(which(MatchScores>=similarityCutoff)),3]
          }
        }
      }
    }
    #Make MissingQMList into dataframe for output
    MissingQMFrame<-do.call(rbind,MissingQMList)
    AlignmentTableList[[as.character(seed)]]<-FinalMatrix
  }

  #If only one seed file is provided just output alignment matrix
  if(length(seedFile)==1){
    FinalMatrix<-AlignmentTableList[[1]]
  }

  #If multiple seed files provided find peaks with >50% overlap across all seed files
  if(length(seedFile)>1){

    #Find all peaks with at >50% alignment overlap
    ConsensusPeaks<-list()
    ConsensusMatches<-list()
    for(i in 1:(length(AlignmentTableList)-1)){
      Overlaps<-apply(AlignmentTableList[[i]], 1, function(y) apply(AlignmentTableList[[length(AlignmentTableList)]], 1, function(x) sum(x%in%y, na.rm = T)))
      Indexes<-arrayInd(which(Overlaps>(length(inputFileList)/2)), dim(Overlaps))
      row.names(Indexes)<-Indexes[,1]
      ConsensusPeaks[[i]]<-Indexes[,1]
      ConsensusMatches[[i]]<-Indexes
    }
    ConsensusPeaks<-Reduce(intersect, ConsensusPeaks)

    #Filter all alignments to consensus peaks (>50% overlap)
    AlignmentTableListFilt<-AlignmentTableList
    for(i in 1:(length(AlignmentTableList)-1)){
      AlignmentTableListFilt[[i]]<-AlignmentTableListFilt[[i]][ConsensusMatches[[i]][as.character(ConsensusPeaks),2],inputFileList]
    }
    AlignmentTableListFilt[[length(AlignmentTableListFilt)]]<-AlignmentTableListFilt[[length(AlignmentTableListFilt)]][ConsensusPeaks,inputFileList]

    #Find median value of peak areas across all alignments
    AlignVector<-data.frame(lapply(AlignmentTableListFilt, as.vector))
    FinalMatrix<-matrix(apply(AlignVector,1,function(x) stats::median(x, na.rm=T)), nrow = length(ConsensusPeaks))
    colnames(FinalMatrix)<-inputFileList
    row.names(FinalMatrix)<- row.names(AlignmentTableListFilt[[length(AlignmentTableListFilt)]])

    #Filter metabolite info file by consensus peaks
    SeedSample[[1]]<-SeedSample[[1]][ConsensusPeaks,]
  }

  #Add metabolite IDs if standardLibrary is used
  if(!is.null(standardLibrary)){
    message("Matching peaks to standard library")

    #Parse seed file spectras
    peakSplit<-split(SeedSample[[1]],1:nrow(SeedSample[[1]]))
    peakSpectraSplit<-lapply(peakSplit, function(a) strsplit(a[[5]]," "))
    peakSpectraSplit<-lapply(peakSpectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    peakSpectraSplit<-lapply(peakSpectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    peakSpectraSplit<-lapply(peakSpectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
    peakSpectraSplit<-lapply(peakSpectraSplit, function(d) apply(d,2,as.numeric))
    peakSpectraSplit<-lapply(peakSpectraSplit, function(d) d[order(d[,1]),2,drop=F])
    peakSpectraFrame<-do.call(cbind,peakSpectraSplit)
    peakSpectraFrame<-t(peakSpectraFrame)
    peakSpectraFrame<-as.matrix(peakSpectraFrame)/sqrt(apply((as.matrix(peakSpectraFrame))^2,1,sum))

    #Parse standard library spectras
    standardSplit<-split(standardLibrary,1:nrow(standardLibrary))
    standardSpectraSplit<-lapply(standardSplit, function(a) strsplit(a[[3]]," "))
    standardSpectraSplit<-lapply(standardSpectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    standardSpectraSplit<-lapply(standardSpectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    standardSpectraSplit<-lapply(standardSpectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
    standardSpectraSplit<-lapply(standardSpectraSplit, function(d) apply(d,2,as.numeric))
    standardSpectraSplit<-lapply(standardSpectraSplit, function(d) d[order(d[,1]),2,drop=F])
    standardSpectraFrame<-do.call(cbind,standardSpectraSplit)
    standardSpectraFrame<-t(standardSpectraFrame)
    standardSpectraFrame<-as.matrix(standardSpectraFrame)/sqrt(apply((as.matrix(standardSpectraFrame))^2,1,sum))
    SimilarityMatrix<-(standardSpectraFrame %*% t(peakSpectraFrame))*100

    #Compute RT differences
    RT1Index<-matrix(unlist(lapply(SeedSample[[1]][,6],function(x) abs(x-standardLibrary[,4])*RT1Penalty)),nrow=nrow(standardLibrary))
    RT2Index<-matrix(unlist(lapply(SeedSample[[1]][,7],function(x) abs(x-standardLibrary[,5])*RT2Penalty)),nrow=nrow(standardLibrary))

    #Use RT indexes to compute RT differences
    if(!is.null(RT1_Standards)){
      RT1Index<-list()
      RT1_Length<-max(SeedSample[[1]][which(SeedSample[[1]][,1]%in%RT1_Standards),6])-min(SeedSample[[1]][which(SeedSample[[1]][,1]%in%RT1_Standards),6])
      for(Standard in RT1_Standards){
        RT1Index[[Standard]]<-matrix(unlist(lapply(SeedSample[[1]][,paste(Standard,"RT1",sep="_")],function(x) abs(x-standardLibrary[,paste(Standard,"RT1",sep="_")])*(RT1Penalty/length(RT1_Standards)))),nrow=nrow(standardLibrary))*RT1_Length
      }
      RT1Index<-Reduce("+",RT1Index)
    }
    if(!is.null(RT2_Standards)){
      RT2Index<-list()
      RT2_Length<-max(SeedSample[[1]][which(SeedSample[[1]][,1]%in%RT2_Standards),6])-min(SeedSample[[1]][which(SeedSample[[1]][,1]%in%RT2_Standards),6])
      for(Standard in RT2_Standards){
        RT2Index[[Standard]]<-matrix(unlist(lapply(SeedSample[[1]][,paste(Standard,"RT2",sep="_")],function(x) abs(x-standardLibrary[,paste(Standard,"RT2",sep="_")])*(RT2Penalty/length(RT2_Standards)))),nrow=nrow(standardLibrary))*RT2_Length
      }
      RT2Index<-Reduce("+",RT2Index)
    }

    #Subtract RT penalties from Similarity Scores
    SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
    row.names(SimilarityMatrix)<-standardLibrary[,1]

    #Append top three ID matches to each metabolite and scores to seedRaw for output
    SeedSample[[1]]<-cbind(t(apply(SimilarityMatrix,2,function(x) paste(names(x[order(-x)])[1:3],round(x[order(-x)][1:3],2),sep="_"))),SeedSample[[1]])
    colnames(SeedSample[[1]])<-c("Library_Match_1","Library_Match_2","Library_Match_3",colnames(SeedSample[[1]])[4:length(colnames(SeedSample[[1]]))])
  }
  returnList<-list(FinalMatrix, SeedSample[[1]], MissingQMFrame)
  names(returnList)<-c("Alignment_Matrix","Peak_Info","Unmatched_Quant_Masses")
  return(returnList)
}
