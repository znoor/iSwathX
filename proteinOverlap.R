################################################################################################
#' Get the overlap proteins between two libraries
#' @param datBaseLib a data frame of a base library 
#' @param datExtLib a data frame of an addon library
#' @param parseAcc a logic value indicating if the protein accession to be parsed. Default value
#' is True.
#' @return A vector of protein accession that the two libraries overlap
#' @examples 
#' libfiles <- paste(system.file("files",package="SwathXtend"),
#' 					c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1])
#' datExtLib <- readLibFile(libfiles[2])
#' accs <- proteinOverlap(datBaseLib,datExtLib)
##################################################################################################               

proteinOverlap <- function(datBaseLib, datExtLib, parseAcc = TRUE){

  base.acc <- unique(datBaseLib$uniprot_id)
  ext.acc <- unique(datExtLib$uniprot_id)
  
  if (parseAcc){
	  base.acc <- sapply(base.acc, parseAccession)
	  ext.acc  <- sapply(ext.acc, parseAccession)
	}	
  # v <- c(rep(FALSE,length(base.acc)))
  
  # for(ii in 1:length(base.acc)){
    # e <- grepl(base.acc[ii], ext.acc)
    
    # if(TRUE %in% e){
      # v[ii] <- TRUE
    # }
    
  # }

  # res <- base.acc[v]
  
  # } else {
    res <- intersect(base.acc, ext.acc)
  # }
  
  res
}



