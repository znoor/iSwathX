###########################################################################
#' Load a spectrum library into a data frame
#' @param file A file of a spectrum library, in .txt or .csv format, can be .gz files.
#' @param format A character string denoting the file format. One of "peakview" (default), "openswath",
#' "skyline" and "spectronaut"
#' @param type A character string denoting the file type. One of "spectrum" (default) and "hydro"
#' @param clean A logic value, representing if the library will be cleaned.
#' @param ... named arguments to be passed to the cleanLib function.
#' @return a data frame of the library with cleaning process
#' @examples
#' file <- paste(system.file("files",package="iSwathX"),"Lib1.txt",sep="/")
#' dat <- readLibFile(file)
############################################################################

# source("canonicalFormat.R")
# source("cleanLib.R")

readLibFile <- function(file, format = c("peakview", "openswath", "skyline", "spectronaut"),
                        type = c("spectrum", "hydro") ,clean = TRUE, ...)
{

  seps = c("\t",",")

  format <- match.arg(format)
  type <- match.arg(type)
  if(type == "spectrum"){
  if (grepl(".txt$",file) || grepl(".txt.gz$",file) || grepl(".tsv$",file))
  {
    id <- 1
  } else if (grepl(".csv$",file) || grepl(".csv.gz$",file) ) {

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

  dat <- try(canonicalFormat(dat, format = format) )

  if(inherits(dat,"try-error"))
    stop("Error with creating canonical format")

  dat <- try(cleanLib(dat, clean=clean, ...))

  if(inherits(dat,"try-error"))
    stop("Error with cleaning library")
  } else { ## Hydro
    dat <- read.delim2(file,sep=seps[1],header=TRUE,stringsAsFactor=FALSE)
  }
  dat

}
