###########################################################################
#' output a spectrum library into a PeakView format file
#' @param dat A data frame of a spectrum library
#' @param filename A character string for the name of the output.
#' @param nodup A logic value, indicating if remove duplicated
#' sprectrum (default)
#' @param format A character string representing the output format. One of
#' "peakview" (default), "openswath", "skyline" and "spectronaut".
#' @return a file with the specified file name (lib.txt as default) will be
#' saved under the current working directory
#' @examples
#' file <- paste(system.file("files",package="iSwathX"),"Lib1.txt",sep="/")
#' dat <- readLibFile(file)
#' outputLib(dat)
############################################################################

outputLib <- function(dat, filename = "NewLib.txt", format = c("peakview", "openswath", "skyline", "spectronaut"), nodup = TRUE)
{
  if(nodup) dat <- dat[!duplicated(dat),]

  dat <- libraryFormat(dat)

  format <- match.arg(format)
  if(format == "peakview"){
    dat.res <- peakviewFormat(dat)
  } else if (format == "openswath") {
    dat.res <- OswathFormat(dat)
  } else if (format == "skyline") {
  dat.res <- skylineFormat(dat)
  } else if(format == "skyline") {
  dat.res <- spectronautFormat(dat)
}   else
    stop("Wrong output format!")


  write.table(dat.res, file = filename,sep="\t",quote=FALSE,row.names=FALSE, na="")

  return(dat.res)
}
