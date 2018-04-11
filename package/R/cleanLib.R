
###########################################################################
#' Spectrum library cleanining
#' @param datLib a data frame for a spectrum library
#' @param nomod a logic value, representing if the modified peptides and its
#' fragment ions will be removed. True (default) means will be removed.
#' @param nomc a logic value, representing if peptides with miss cleavages
#' are removed. Default vaue is False (not to remove).
#' @param intensity.cutoff A number value to specify cut off for relative
#' intensity of fragment ions. Only ions with intensity higher than the
#' cut off value (default as 5) will be kept.
#' @param conf.cutoff A number value to specify cut off for precursor
#' confidence. Only ions with confidence higher than the cut off value
#' (default as 0.99) will be kept.
#' @param enz A character string representing the enzyme which can be one of
#' "trypsin" (defalut), "gluc", or "chymotrypsin"
#' @return a data frame of a cleaned spectrum library by the specified
#' criteria
#' @examples
#' file <- paste(system.file("files",package="iSwathX"),"Lib1.txt",sep="/")
#' dat <- read.delim2(file,sep="\t",header=TRUE,stringsAsFactors=FALSE)
#' dat <- canonicalFormat(dat)
#' dat <- cleanLib(dat)
############################################################################

#
# source("isMissCleavaged.R")
cleanLib <- function(datLib, clean=TRUE, intensity.cutoff=5,conf.cutoff=0.99,
                    nomod = TRUE, nomc = FALSE, enz=c("trypsin", "gluc", "chymotrypsin") )
{

  enz <- match.arg(enz)

  if(clean)
    dat.res <- datLib[as.numeric(datLib$relative_intensity)>
                       intensity.cutoff &
                       as.numeric(datLib$confidence) >
                       conf.cutoff,]

  else   dat.res <- datLib


  if(nomod)
    if(length(grep("\\[",dat.res$modification_sequence))>0)
      dat.res <- dat.res[-grep("\\[",dat.res$modification_sequence),]
  if(nomc)

    dat.res <- dat.res[-sapply(dat.res$stripped_sequence,
                               FUN=function(x){isMissCleavaged(x, enz)}),]

  dat.res
}
