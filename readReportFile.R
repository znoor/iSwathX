###########################################################################
#' Load a report file into the dataframe
#' @param file A file of a report file, in .txt or .csv format.
#' @return a data frame of the report file
#' @examples 
#' file <- paste(system.file("files",package="iSwathX"),"Report_file.txt",sep="/")
#' dat <- readReportFile(file)
############################################################################


readReportFile <- function(file)
{
  seps = c("\t", ",")
  
  if (grepl(".txt$",file) || grepl(".tsv$",file) ) 
  { 
    id <- 1
  } else if (grepl(".csv$",file)) {
    
    id <- 2 
  }  else {
    stop (paste(file, "have wrong file type !") )
  }
  
  sep <- seps[id]
  
  dat <- read.delim2(file,sep=sep,header=TRUE,
                     stringsAsFactors=FALSE)
  
  if(ncol(dat) < 2){
    stop(paste(file, "wrong library file format!") )
  }
  
  dat
}