###########################################################################
#' Convert a data frame into Peak View library format
#' @param datlib a data frame of a canonical spectrum library
#' @param nodup A logic value, indicating if remove duplicated
#' sprectrum (default)
#' @return a data frame of a spectrum library in Peak View format
#' @examples
#' libfile <- paste(system.file("files",package="SwathXtend"),"Lib2.txt",sep="/")
#' dat <- readLibFile(libfile)
#' dat <- peakviewFormat(dat)
############################################################################

### LIBRARY FORMAT FILE (PEAK VIEW)

peakviewFormat <- function(datlib, nodup = TRUE)
{
  
  
  
  if(nodup) datlib <- datlib[!duplicated(datlib),]
  
  
  inclcols = c("Q1",
               "Q3",
               "RT_detected",
               "protein_name",
               "isotype",
               "relative_intensity",
               "stripped_sequence",
               "modification_sequence",
               "prec_z",
               "frg_type",
               "frg_z",
               "frg_nr",
               "iRT",
               "uniprot_id",
               "score",
               "decoy",
               "prec_y",
               "confidence",
               "shared",
               "N",
               "rank",
               "mods",
               "nterm",
               "cterm")
  
  if( all(inclcols %in% colnames(datlib)) )	{	
    res = datlib[,inclcols]
  }
  else {
    miss.col = inclcols[!inclcols %in% colnames(datlib)]
    
    toadd.col = matrix(NA, ncol=length(miss.col), nrow=nrow(datlib))
    colnames(toadd.col) = miss.col		
    
    res = cbind(datlib, toadd.col)
    
    # check if all columns are included
    
    stopifnot(all(inclcols %in% colnames(res)))
    
    # re-order
    res = res[,inclcols]
    
  }
  res
}