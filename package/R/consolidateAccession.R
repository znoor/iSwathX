###########################################################################
#' Consolidate the protein accessions between two libraries
#' @param datBaseLib a data frame for base library
#' @param datExtLib a data frame for external/addon library
#' @return a data frame of the external library with newly consolidated
#' protein accessions based on the base library.
#' @examples
#' file1 <- paste(system.file("files",package="iSwathX"),"Lib1.txt",sep="/")
#' file2 <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' datBaseLib <- readLibFile(file1)
#' datExtLib <- readLibFile(file2)
#' dat <- consolidateAccession(datBaseLib, datExtLib)
############################################################################

consolidateAccession <- function(datBaseLib, datExtLib)
{

  for(x in unique(datBaseLib$uniprot_id)){

      xx <- parseAccession(x)
      id.toreplace <- grepl(xx,datExtLib$uniprot_id)


    if(length(which(id.toreplace) > 0)) {

      datExtLib$uniprot_id[id.toreplace] <- x
    }

  }



  datExtLib
}

