###########################################################################
#' Build a spectra library by integrating a pair of spectrum libraries
#' @param baseLib A base library data frame or file
#' @param extLib An external/addon library data frame or file
#' @param hydroIndex A data frame or file containing peptide hydrophobicity index
#' @param method A character string to specify the RT alignment method. One of "time" (default),
#' "hydro", "hydrosequence" can be selected
#' @param cutoff.size A number value specifying the mimum number of common
#' peptides between two libraries requried for training
#' @param cutoff.r2 A number value specifying the mimum correlation coefficient
#' of retention time that this method requires.
#' @param recalibrateRT A logic value indicating if retention times of external/addon library
#' adjust/shift in case of poor correlation
#' @param includeLength A logic value indicating if include the peptide
#' length as a feature for predicting retention time. Only applicable when method is "hydro"
#' @param labelBase A character string to specify the labels of proteins from the base library
#' @param labelAddon A character string to specify the labels of proteins from addon library
#' @param formatBase A character string denoting the file format of base library file. One of
#' "peakview" (default), "openswath", "skyline" and "spectronaut".
#' @param formatExt A character string denoting the file format of addon library file. One of
#' "peakview" (default), "openswath", "skyline" and "spectronaut".
#' @param outputFormat A character string denoting the file format of the output integrated library. One of
#' "peakview" (default), "openswath", "skyline" and "spectronaut".
#' @param outputFile A character string to specify the spectra library created.
#' @param plot a logic value representing if plots during processing will be plotted or not
#' @param clean A logic value representing if the input librarues will be cleaned before integration.
#' Default value is True
#' @param merge A logic value representing if the output will be the merged library (default),
#' or the adjested add-on library
#' @param consolidateAccession A logic value representing if the protein accessions will be consolidated to the
#' base library in the integrated library. Default value is True.
#' @return A data frame of integrated spectrum library, also saved in working directory. Library analysis plots saved in
#' working directory
#' @examples
#' libfiles <- paste(system.file("files",package="iSwathX"),
#'                              c("Lib2.txt","Lib3.txt"),sep="/")
#' Lib2_3 <- buildSpectraLibPair(libfiles[1], libfiles[2], outputFormat = "peakview", clean = T)
############################################################################


source("readLibFile.R")
source("normalise.R")
source("splitLib.R")
source("predictRT.R")
source("alignRTbyHydro.R")
source("libraryFormat.R")
source("consolidateAccession.R")
source("parseAccession.R")
source("plotStats.R")


buildSpectraLibPair<-function(baseLib,extLib, hydroIndex,
                           method = c("time", "hydro", "hydrosequence"),
                           cutoff.size = 50,
                           cutoff.r2 = 0.8,
                           recalibrateRT = TRUE,
                           includeLength = FALSE,
                           labelBase = NA,
                           labelAddon = NA,
                           formatBase = c("peakview", "openswath", "skyline", "spectronaut"),
                           formatExt = c("peakview", "openswath", "skyline", "spectronaut"),
                           outputFormat = c("peakview", "openswath", "skyline", "spectronaut"),
                           outputFile = "extendedLibrary.txt", 
                           plot = FALSE, 
                           clean = TRUE, 
                           merge = TRUE,
                           parseAcc = FALSE,
                           consolidateAccession=FALSE, ...)
{
  
  if(!is.data.frame(baseLib)){	
	  datBaseLib <- try(readLibFile(baseLib, format = formatBase, clean=clean, ...))
	  
	  if(inherits(datBaseLib, "try-error"))
	  stop(paste("Error with reading", baseLib))  
	} else {
		datBaseLib <- baseLib
	}
  if(!is.data.frame(extLib)){  
  datExtLib <- try(readLibFile(extLib, format = formatExt, clean=clean, ...))
  
  if(inherits(datExtLib,"try-error"))
    stop(paste("Error with reading", extLib))
  
  } else {
		datExtLib <- extLib
	}
  

  
  datBaseLib <- try(normalise(datBaseLib))
  
  if(inherits(datBaseLib,"try-error"))
	stop("Error with normalising datBaseLib")	

  datExtLib <- try(normalise(datExtLib))
  

  if(inherits(datExtLib,"try-error"))
	  stop("Error with normalising datExtLib")
	
  ## split datExtLib into common and new  
  list.datLibs <- try(splitLib(datBaseLib, datExtLib))
  

  if(inherits(list.datLibs,"try-error"))
	  stop("Error with splitting datExtLib")
	
  datExtLibCommPart <- list.datLibs[["ExtCommon"]]
  datExtLibNewPart <- list.datLibs[["ExtNew"]]
  datBaseCommPart <- list.datLibs[["BaseCommon"]]
 
  

  if(nrow(datBaseCommPart[!duplicated(datBaseCommPart$stripped_sequence),]) < 20)
  {  
    stop("Library merging not executed because common peptides are less than 20.") 
	
	res <- datBaseLib
  } else {
	  
	  if(plot) 
		
	    pall <-plotStats(datBaseLib, datExtLib)

	  datBaseCommPart <- libraryFormat(datBaseCommPart)
	  datExtLibCommPart <- libraryFormat(datExtLibCommPart)
	  datExtLibNewPart <- libraryFormat(datExtLibNewPart)
	  
	  
	  datRTrain <- getRTrain(datBaseCommPart, datExtLibCommPart) 
	  
	  
	  if(nrow(datRTrain) < cutoff.size) 
	    warnings(paste("Training set size is less than", cutoff.size))
	  
	  if(cor(datRTrain$RT_detected, datRTrain$RT_detected.base)^2 
	     < cutoff.r2){	
	    warnings(paste("R Squared of the RT training set is < ",
	                   cutoff.r2))
	    if(recalibrateRT){
	      
	      difference <- subtract(datRTrain$RT_detected.base, datRTrain$RT_detected)
	      
	      toadd <- mean(difference)
	      
	      if(sd(difference)<0.5)
	        datExtLib <- mutate(datExtLib, RT_detected = RT_detected + toadd)
	      
	      ## split datExtLib into common and new  
	      list.datLibs <- try(splitLib(datBaseLib, datExtLib))
	      
	      
	      if(inherits(list.datLibs,"try-error"))
	        stop("Error with splitting datExtLib")
	      
	      datExtLibCommPart <- list.datLibs[["ExtCommon"]]
	      datExtLibNewPart <- list.datLibs[["ExtNew"]]
	      datBaseCommPart <- list.datLibs[["BaseCommon"]]
	    }
	  }
	  
	 
		method <- match.arg(method)
	  if(method == "time"){
		
	  
		datExtLibNewPart <- try(predictRT(datBaseCommPart,
										 datExtLibCommPart, 
										  datExtLibNewPart,cutoff.size = cutoff.size, cutoff.r2 = cutoff.r2))
		
		if(inherits(datExtLibNewPart,"try-error"))
		  stop("Error with predictRT")
		
	  } else if (grepl("hydro",method)) {
		
	    if(!is.data.frame(hydroIndex)){
		datHydroIndex <- try(readLibFile(hydroIndex, type="Hydro"))
		
		if(inherits(datHydroIndex, "try-error"))
			stop("Error with reading hydroindex file")
		} else {
			datHydroIndex <- hydroIndex
		}
		
		datExtLibNewPart <- try(alignRTbyHydro(datBaseLib,
											   datExtLibNewPart, 
											   datHydroIndex,
											   method=method,
											   includeLength=includeLength))
		if(inherits(datExtLibNewPart,"try-error"))
		  stop("Error with alignRTbyHydro")
	  } else stop(paste("Unknow method ", method))
		  

	  
	  ## RE-GROUPING by common uniprot accession

	  datBaseCommPart <- libraryFormat(datBaseCommPart)
	  
	  datBaseLib <- libraryFormat(datBaseLib)
	  
	  
	  datExtLibNewPart <- libraryFormat(datExtLibNewPart)
	  
	  if(consolidateAccession)
	  datExtLibNewPart <- 
		consolidateAccession(datBaseLib, datExtLibNewPart)
	  

	  ## Merge base with aligned external new part
	  
	  datExtLibNewPart = libraryFormat(datExtLibNewPart)  
	  
	  if(!is.na(labelBase)){
		datBaseLib$uniprot_id <- paste(datBaseLib$uniprot_id, labelBase, sep="_")
	  }
	  
	  if(!is.na(labelAddon)){
		datExtLibNewPart$uniprot_id <- paste(datExtLibNewPart$uniprot_id, labelAddon, sep="_")
	  }
	  
	  res <- datExtLibNewPart
	  
	  if(merge) {
		# datlib.fin = rbind(datBaseLib,datExtLibNewPart)    
		# 
		# datlib.fin = libraryFormat(datlib.fin)
	    
	    res = rbind(datBaseLib,datExtLibNewPart)    
	    
	    reli <- reliabilityCheckLibrary(datBaseLib, res)
		
	  }
	  res = libraryFormat(res)
  }
  if(plot)
  return(list(res, pall))
  else
    return(list(res))
 
}  
  