########################################################################################
#' Get the training set form retention time prediction
#' @param dat1 a data frame of the base library
#' @param dat2 a data frame of the add-on library
#' @param nomod a logic value indicating if the modified peptides are included
#' when generating the training set, default value is T, i.e., no modified peptides
#' included
#' @return a data frame contained the retention time training set 
#' @examples 
#' libfiles <- paste(system.file("files",package="iSwathX"),
#'    c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1], clean=TRUE, nomod=FALSE, nomc=FALSE)
#' datExtLib <- readLibFile(libfiles[2], clean=TRUE, nomod=FALSE, nomc=FALSE)
#' datRTrain <- getRTrain(datBaseLib, datExtLib, nomod=TRUE)  
#########################################################################################

getRTrain <- function(dat1, dat2, nomod=TRUE)
{
  if(nomod) {
    
    
    dat1 <- dat1[!duplicated(dat1$stripped_sequence),]
    dat2 <- dat2[!duplicated(dat2$stripped_sequence),]
    
    datRTrain <- merge(dat1[,c("stripped_sequence","RT_detected")],
                       dat2[,c("stripped_sequence","RT_detected")],                       
                       by=1)
    
  } else { # include modified peptides    
    
    
    dat1$pepmod <- paste(dat1$stripped_sequence, dat1$modification_sequence)
    dat2$pepmod <- paste(dat2$stripped_sequence, dat2$modification_sequence)
    
    datRTrain <- merge(dat1[!duplicated(dat1$pepmod),c("pepmod", "RT_detected")],
                       dat2[!duplicated(dat2$pepmod),c("pepmod", "RT_detected")],
                       by=1)
  }
  
  colnames(datRTrain) <- c("sequence",
                           "RT_detected.base","RT_detected")
  
  
  datRTrain
}