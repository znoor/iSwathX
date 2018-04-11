###########################################################################
#' Convert a data frame into OpenSWATH library format
#' @param dat a data frame of a canonical spectrum library
#' @param nodup A logic value, indicating if remove duplicated
#' sprectrum (default)
#' @return a data frame of a spectrum library in OpenSWATH format
#' @examples
#' libfile <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' dat <- readLibFile(libfile)
#' dat <- OswathFormat(dat)
############################################################################

### LIBRARY FORMAT FILE (OPENSWATH)
#
# source("libraryFormat.R")

OswathFormat <- function(dat, nodup = TRUE)
{
  if(nodup) dat <- dat[!duplicated(dat),]

  dat <- libraryFormat(dat)

  PrecursorMz  <- as.numeric(as.character(dat$Q1))
  ProductMz  <- as.numeric(as.character(dat$Q3))
  Tr_recalibrated  <- as.numeric(as.character(dat$RT_detected))
  ProteinName  <- dat$protein_name
  GroupLabel  <- as.factor(dat$isotype)
  LibraryIntensity  <- as.numeric(as.character(dat$relative_intensity))
  PeptideSequence  <- dat$stripped_sequence
  FullUniModPeptideName  <- dat$modification_sequence
  UniprotID  <- dat$uniprot_id
  decoy  <- as.factor(dat$decoy)
  PrecursorCharge  <- as.numeric(as.character(dat$prec_z))
  FragmentType  <- dat$frg_type
  FragmentCharge  <- as.numeric(as.character(dat$frg_z))
  FragmentSeriesNumber  <- as.numeric(as.character(dat$frg_nr))


  dat.res <- data.frame(PrecursorMz,ProductMz,Tr_recalibrated,
                        ProteinName,GroupLabel,LibraryIntensity,
                        PeptideSequence,FullUniModPeptideName,
                        UniprotID,decoy,PrecursorCharge,
                        FragmentType,FragmentCharge,FragmentSeriesNumber)


  dat.res
}
