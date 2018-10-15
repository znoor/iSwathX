###########################################################################
#' Convert a data frame into Spectronaut library format
#' @param dat a data frame of a canonical spectrum library
#' @param nodup A logic value, indicating if remove duplicated
#' sprectrum (default)
#' @return a data frame of a spectrum library in Spectronaut format
#' @examples
#' libfile <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' dat <- readLibFile(libfile)
#' dat <- spectronautFormat(dat)
############################################################################

### LIBRARY FORMAT FILE (SPECTRONAUT)

source("libraryFormat.R")

spectronautFormat <- function(dat, nodup = TRUE)
{
  if(nodup) dat <- dat[!duplicated(dat),]
  
  dat <- libraryFormat(dat)
  
  PrecursorMz  <- as.numeric(as.character(dat$Q1))
  FragmentMz  <- as.numeric(as.character(dat$Q3))
  RetentionTime  <- as.numeric(as.character(dat$RT_detected))
  iRT  <- as.numeric(as.character(dat$iRT))
  RelativeIntensity  <- as.numeric(as.character(dat$relative_intensity))
  StrippedPeptide  <- dat$stripped_sequence
  PrecursorCharge  <- as.numeric(as.character(dat$prec_z))
  FragmentType  <- dat$frg_type
  FragmentNumber  <- as.numeric(as.character(dat$frg_nr))
  FragmentCharge  <- as.numeric(as.character(dat$frg_z))
  FragmentLossType <- c(rep("", nrow(dat)))
  ExcludeFromAssay <- c(rep(FALSE, nrow(dat)))
  ModifiedPeptide  <- dat$modification_sequence
  LabeledSequence <- dat$modification_sequence
  ProteinGroups  <- dat$protein_name
  UniProtIds  <- dat$uniprot_id
  UserGroup <- c(rep("", nrow(dat)))
  
  dat.res <- data.frame(PrecursorMz,FragmentMz,RetentionTime, iRT,
                        RelativeIntensity, StrippedPeptide, PrecursorCharge,
                        FragmentType, FragmentNumber, FragmentCharge, FragmentLossType,
                        ExcludeFromAssay, ModifiedPeptide, LabeledSequence,
                        ProteinGroups, UniProtIds, UserGroup)
  
  dat.res
  
}  