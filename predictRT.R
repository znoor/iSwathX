###############################################################################
#' Predicting retention time based on time
#' @param datBaseCommPart a data frame of common precursors in the base library
#' @param datExtLibCommPart a data frame of common precursors in the 
#' external library
#' @param datExtLibNewPart a data frame of new spectra in the external library
#' @param cutoff.size a number value specifying the mimum number of common 
#' peptides between two libraries requried for training
#' @param cutoff.r2 a number value specifying the mimum correlation coefficient
#' of retention time that this method requires.
#' @return a data frame which has the same spectra as datExtLibNewPart but with
#' new predicted retention time
#' @examples 
#' libfiles <- paste(system.file("files",package="SwathXtend"),
#'                              c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1], clean=TRUE, nomod=FALSE, nomc=FALSE)
#' datExtLib <- readLibFile(libfiles[2], clean=TRUE, nomod=FALSE, nomc=FALSE) 
#' datExtLibNewPart <- splitLib(datBaseLib, datExtLib, nomod=TRUE)[["ExtNew"]]
#' datExtLibNewPart <- 
#'        predictRT(datBaseLib,datBaseLib,datExtLibNewPart)
#'  
###############################################################################
source("libraryFormat.R")
source("selectModel.R")

predictRT <- function(datBaseCommPart,
                      datExtLibCommPart,  
                      datExtLibNewPart,
						cutoff.size = 50,
						cutoff.r2 = 0.8,
						cutoff.rse = 3.0)
{

  datBaseCommPart <- libraryFormat(datBaseCommPart)
  datExtLibCommPart <- libraryFormat(datExtLibCommPart)
  datExtLibNewPart <- libraryFormat(datExtLibNewPart)
  
  
  datRTrain <- getRTrain(datBaseCommPart, datExtLibCommPart) 


  
  if(nrow(datRTrain) < cutoff.size) 
    warnings(paste("Training set size is less than", cutoff.size))
  
  if(cor(datRTrain$RT_detected, datRTrain$RT_detected.base)^2 
     < cutoff.r2)	
    warnings(paste("R Squared of the RT training set is < ",
				cutoff.r2))
  if(summary(lm(RT_detected.base ~ RT_detected, data=datRTrain))$sigma > 
       cutoff.rse )
    warnings(paste("Residual SE is greater than  ",
                   cutoff.rse))
  
  bestModel <- try(selectModel(datRTrain)   )
  
  if(!is.null(bestModel) ) {
  
	datExtLibNewPart$RT_detected <- predict(bestModel, datExtLibNewPart[,c(1:3)])
  
  } else {
	warning("Retention time alignment failed due to lacking common peptides.")
  }
   
  datExtLibNewPart
}



