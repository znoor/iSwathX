###########################################################################
#' Convert a data frame into Skyline library format
#' @param dat a data frame of a canonical spectrum library
#' @param nodup A logic value, indicating if remove duplicated
#' sprectrum (default)
#' @return a data frame of a spectrum library in Skyline format
#' @examples
#' libfile <- paste(system.file("files",package="SwathXtend"),"Lib2.txt",sep="/")
#' dat <- readLibFile(libfile)
#' dat <- skylineFormat(dat)
############################################################################

### LIBRARY FORMAT FILE (SKYLINE)

source("libraryFormat.R")

skylineFormat <- function(dat, nodup = TRUE)
{
  if(nodup) dat <- dat[!duplicated(dat),]
  
  dat <- libraryFormat(dat)
  
  ProteinName  <- dat$protein_name
  PeptideSequence <- dat$stripped_sequence
  ModificationSequence <- dat$modification_sequence
  UniprotID  <- dat$uniprot_id
  Tr_recalibrated  <- as.numeric(as.character(dat$RT_detected))
  PrecursorMz  <- as.numeric(as.character(dat$Q1))
  ProductMz  <- as.numeric(as.character(dat$Q3))
  LibraryIntensity  <- as.numeric(as.character(dat$relative_intensity))
  PrecursorCharge  <- as.numeric(as.character(dat$prec_z))
  FragmentType  <- dat$frg_type
  FragmentCharge  <- as.numeric(as.character(dat$frg_z))
  FragmentSeriesNumber  <- as.numeric(as.character(dat$frg_nr))
  decoy  <- as.factor(dat$decoy)
  
  
  dat.res <- data.frame(PrecursorMz,ProductMz,Tr_recalibrated,
                        ProteinName,LibraryIntensity, PeptideSequence,
                        ModificationSequence,UniprotID,
                        PrecursorCharge,FragmentType,FragmentCharge,FragmentSeriesNumber, decoy)
  
  dat.res
  
}  