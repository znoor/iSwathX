###########################################################################
#' Generate a library summary 
#' @param datlib a data frame for a spectrum library
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

libSummary <- function(datlib)
{
  summarylist <- list()
  
  # datlib <- canonicalFormat(datlib, format = format)
  proteinCount <- datlib[!duplicated(datlib$uniprot_id),]
  unmodpeptidesCount <- datlib[!duplicated(datlib$stripped_sequence),]
  modpeptidesCount <- datlib[!duplicated(datlib$modification_sequence),]
  #peptidesCount <- datlib[!duplicated(datlib$stripped_sequence),]
  #TransitionCount <- datlib[!duplicated(datlib$stripped_sequence),]
  
  if(all(is.na(proteinCount$protein_name))) {
    summarylist[["proteins"]] = paste("Proteins names are not available.")
    summarylist[["unmodpeptides"]] = length(unmodpeptidesCount$stripped_sequence)
    summarylist[["modpeptides"]] = length(modpeptidesCount$modification_sequence)
    summarylist[["transitions"]] = length(datlib$frg_type)
  } else {
  
  
  summarylist[["proteins"]] = length(proteinCount$uniprot_id)
  summarylist[["unmodpeptides"]] = length(unmodpeptidesCount$stripped_sequence)
  summarylist[["modpeptides"]] = length(modpeptidesCount$stripped_sequence)
  summarylist[["transitions"]] = length(datlib$frg_type)
  }
  
  summarylist
}