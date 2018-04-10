
###########################################################################
#' Normalisation of relative ion intensity
#' @param dat a data frame of a spectrum library
#' @param method a character string representing the normalisation method.
#' One of "top" (default) or "total".
#' @return a data frame with normalised ion intensity
#' @examples 
#' file <- paste(system.file("files",package="SwathXtend"),"Lib2.txt",sep="/")
#' dat <- normalise(readLibFile(file))
############################################################################     

normalise = function(dat, method=c("top", "total"))
{
  method <- match.arg(method)
  if(method == "top") {
    ## normalise to the top peak intensity
    dat.ref <- aggregate(dat$relative_intensity, 
                        by=list(dat$uniprot_id,dat$stripped_sequence),
                        max)
    dat$propep <- paste(dat$uniprot_id,dat$stripped_sequence)

    names(dat.ref) <- c("uniprot_id","stripped_sequence","max_intensity")
    
    dat.ref$propep <- paste(dat.ref$uniprot_id,dat.ref$stripped_sequence)    
    
    dat <- merge(dat,dat.ref[,c("propep","max_intensity")],by="propep")
    
    dat$relative_intensity <- dat$relative_intensity/dat$max_intensity*10000  

    selcols <- length(colnames(dat)) 
    
    dat <- dat[,c(2:(selcols-1))]      
    

  } else { # total
    dat.ref <- aggregate(dat$relative_intensity, 
                         by=list(dat$uniprot_id,dat$stripped_sequence),
                         sum)
    dat$propep <- paste(dat$uniprot_id,dat$stripped_sequence)
    
    names(dat.ref) <- c("uniprot_id","stripped_sequence","total_intensity")
    
    dat.ref$propep <- paste(dat.ref$uniprot_id,dat.ref$stripped_sequence)    
    
    dat <- merge(dat,dat.ref[,c("propep","total_intensity")],by="propep")
    
    dat$relative_intensity <- dat$relative_intensity/dat$total_intensity*10000  
    
    selcols <- length(colnames(dat)) 
    
    dat <- dat[,c(2:(selcols-1))]  
  }
  
  dat
}