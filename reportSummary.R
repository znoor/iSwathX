###########################################################################
#' Generate a report summary 
#' @param dat a data frame for a spectrum library
#' @return a list containing information related to proteins, peptides and transitions in the reference library
#' @examples
#' libfile <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' datlib <- readLibFile(libfile)
#' summary <- libSummary(datlib)
############################################################################

########## Library Summary function
# install.packages("plyr")
library("plyr")
source("canonicalFormat.R")


reportSummary <- function(dat)
{
  summarylist <- list()
  
  proteinCount <- dat[!duplicated(dat$Protein.Name),]
  unmodpeptidesCount <- dat[!duplicated(dat$Peptide),]
  modpeptidesCount <- dat[!duplicated(dat$Modified.Sequence),]
  
  summarylist[["proteins"]] = length(proteinCount$Protein.Name)
  summarylist[["unmodpeptides"]] = length(unmodpeptidesCount$Peptide)
  summarylist[["modpeptides"]] = length(modpeptidesCount$Modified.Sequence)
  
  summarylist
  
}
