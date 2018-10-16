###########################################################################
#' Aligning retention time using hydrophobicity index
#' @param dat1 A data frame of a spectrum library to be aligned to
#' @param dat2 A data frame containing the spectra to be predicted
#' @param datHydroIndex a data frame containing all the hydrophobicity
#' index of the peptides in both dat1 and dat2
#' @param includeLength a logic value indicating if including the peptide
#' length as a training feature (default is False)
#' @param method a string character indicating the method of RT alignment.
#' One of "HydroSequence" (default) and "Hydro".
#' @return a data frame of the dat2 with the newly aligned retention time
#' @examples
#' file1 <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' file2 <- paste(system.file("files",package="iSwathX"),"Lib3.txt",sep="/")
#' file3 <- paste(system.file("files",package="iSwathX"),"hydroIndex.txt",sep="/")
#' dat1 <- normalise(readLibFile(file1))
#' dat2 <- normalise(readLibFile(file2))
#' datHydroIndex <- read.delim2(file3,sep="\t",header=TRUE,as.is=TRUE)
#' dat2 <- alignRTbyHydro(dat1, dat2, datHydroIndex,includeLength=FALSE)
############################################################################

source("libraryFormat.R") 
# install.packages("e1071")
library(e1071)


alignRTbyHydro <- function(dat1,
                           dat2,
                           datHydroIndex,includeLength=FALSE,
                           method=c("hydrosequence", "hydro") )                      
  
{

	method <- match.arg(method)
  
  if(method=="hydro"){
    colnames(datHydroIndex)[grep("sequence",tolower(colnames(datHydroIndex)))]<-
      "sequence"
    
    colnames(datHydroIndex)[grep("hydro",tolower(colnames(datHydroIndex)))] <- "hydro"
    
    colnames(datHydroIndex)[grep("len",tolower(colnames(datHydroIndex)))] <- "length"
    
    dat1 <- dat1[!duplicated(dat1$stripped_sequence),]
    
    datHydroIndex <- datHydroIndex[!duplicated(datHydroIndex$sequence),]
    
    datHydroIndex$hydro <- as.numeric(datHydroIndex$hydro)
    
    dat.train <- merge(dat1[,c("stripped_sequence","RT_detected")],
                       datHydroIndex[,c("sequence","hydro","length")],
                       by.x="stripped_sequence",
                       by.y=tolower("sequence"),all=FALSE)
    
    dat.train$length <- log(dat.train$length)
    
    if(includeLength){
      f= as.formula("RT_detected~hydro+length")
      
      idx<-match(dat2$stripped_sequence,datHydroIndex$sequence, nomatch=0)
      
      if(0 %in% idx){
        warning("Some peptides have no hydro index found.")
      }
      
      
      dat2$hydro[idx!=0] <- datHydroIndex[idx,"hydro"]
      dat2$length[idx!=0] <- log(datHydroIndex[idx,"length"])
      
    } else {
      f= as.formula("RT_detected~hydro")
      
      idx<-match(dat2$stripped_sequence,datHydroIndex$sequence, nomatch=0)
      
      if(0 %in% idx){
        warning("Some peptides have no hydro index found.")
      }
      
      
      dat2$hydro[idx!=0] <- datHydroIndex[idx,"hydro"]  
      
      
    } 
    
    
    fit <- lm(f, data=dat.train)
    
    
    dat2$RT_detected <-predict(fit,dat2)
    
    
    dat2 <- libraryFormat(dat2)
    
    
  } else if(method=="hydrosequence"){
      
      colnames(datHydroIndex)[grep("sequence",tolower(colnames(datHydroIndex)))]<-
        "sequence"
      
      colnames(datHydroIndex)[grep("hydro",tolower(colnames(datHydroIndex)))] <- "hydro"
      
      colnames(datHydroIndex)[grep("len",tolower(colnames(datHydroIndex)))] <- "length"
      
      dat1 <- dat1[!duplicated(dat1$stripped_sequence),]
      
      datHydroIndex <- datHydroIndex[!duplicated(datHydroIndex$sequence),]
      
      datHydroIndex$hydro <- as.numeric(datHydroIndex$hydro)
      
      datRTTrain <- merge(dat1[,c("stripped_sequence","RT_detected")],
                          datHydroIndex[,c("sequence","hydro","length")],
                          by.x="stripped_sequence",
                          by.y=tolower("sequence"),all=FALSE)
      
      datRTTrain$length <- log(datRTTrain$length)
      
      datRTTrain.aa <-  do.call(rbind, lapply(unique(datRTTrain$stripped_sequence), 
                                              convertAACount))
      
      datRTTrain.aa <- merge(datRTTrain.aa, datRTTrain, by=1)    
      
      
      fit.hydroNseq.svm <- svm(RT_detected ~., datRTTrain.aa[,-1],scale=TRUE,
                               gamma=0.01,cost=10)    
      
      ## PREDICT FOR NEW    
      
      dat2.aa <- do.call(rbind,lapply(unique(dat2$stripped_sequence),
                                      convertAACount))
      idx<-match(dat2.aa$sequence,datHydroIndex$sequence, nomatch=0)
      
      if(0 %in% idx)
        warning("Some peptides have no hydro index found.")
      
      dat2.aa$hydro[idx!=0] <- datHydroIndex[idx,"hydro"]  
      dat2.aa$length[idx!=0] <- log(datHydroIndex[idx,"length"])    
      
      
      x.col <- grep("sequence",colnames(dat2.aa))
      
      dat2.aa$predicted.rt <- predict(fit.hydroNseq.svm, newdata=dat2.aa)
      
      m <- merge(dat2,dat2.aa[,c(x.col,ncol(dat2.aa))], 
                 by.x="stripped_sequence",by.y=x.col,all.x=TRUE)    
      
      
      dat2$RT_detected <-m$predicted.rt
      
    } else {
    stop("Unknow alignment method!")
  }
  
  dat2  
}