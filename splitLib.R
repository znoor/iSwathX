###########################################################################
#' Split spectra libraries into common and new spectra
#' @param datBaseLib a data frame for base library
#' @param datExtLib a data frame for external/addon library
#' @param alignSpecies a character string or NA (default) representing the
#' species of the alignment base in the base library 
#' @param ignoreCharge a logic value indicating if using precursor charge
#' states as a splitting factor
#' @return a list of data frames composed of BaseCommon, ExtCommon and
#' ExtNew, corresponding to the common spectra of base library, common
#' spectra of external library and new spectra of external library 
#' @examples 
#' libfiles <- paste(system.file("files",package="SwathXtend"),
#'                              c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1], clean=TRUE, nomod=FALSE, nomc=FALSE)
#' datExtLib <- readLibFile(libfiles[2], clean=TRUE, nomod=FALSE, nomc=FALSE) 
#' list.datLibs <- splitLib(datBaseLib, datExtLib, nomod=FALSE)
############################################################################     

splitLib = function(datBaseLib, datExtLib, alignSpecies = NA)
{
  list.res = list()
  
  if(!is.na(alignSpecies))
    datBaseLib <- datBaseLib[grep(as.character(alignSpecies), datBaseLib$uniprot_id),]  
   
    datBaseLib$pepmod <- paste(datBaseLib$stripped_sequence, datBaseLib$modification_sequence)
    datExtLib$pepmod <- paste(datExtLib$stripped_sequence, datExtLib$modification_sequence)
    
  
    commpepmod <- intersect(datBaseLib$pepmod, datExtLib$pepmod)
    
    if(length(commpepmod) == 0)
      warning("The base and add-on libraries does NOT have any common peptides!\n")
    
    
    dat.extcomm <- datExtLib[datExtLib$pepmod %in% commpepmod,]
    
    dat.extnew <- datExtLib[!datExtLib$pepmod %in% commpepmod,]
    
    dat.basecomm <- datBaseLib[datBaseLib$pepmod %in% commpepmod,]    
  
    
    list.res[["ExtCommon"]] = dat.extcomm
    list.res[["ExtNew"]] = dat.extnew
    list.res[["BaseCommon"]] = dat.basecomm
    
    list.res
  

  
}