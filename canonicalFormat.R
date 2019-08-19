
###########################################################################
#' Standardise a sprectrum library data frame 
#' @param dat a data frame of a spectrum library
#' @param format a character string, representing the format of the input
#' spectrum library. One of "peakview" (default), "openswath", "skyline" and "spectronaut"
#' #' @return a data frame of the library in canonical format
#' @examples
#' file <- paste(system.file("files",package="iSwathX"),"Lib1.txt",sep="/")
#' dat <- read.delim2(file,sep="\t",stringsAsFactor = FALSE, header=TRUE)
#' dat <- try(canonicalFormat(dat, format = "peakview"))
############################################################################ 

canonicalFormat = function(dat, format=c("peakview", "openswath", "skyline", "spectronaut") )
{
  format <- match.arg(format) 
  if(format == "peakview"){
  
    dat.res <- dat
    dat.res$Q1 <- as.numeric(as.character(dat.res$Q1))
    dat.res$Q3 <- as.numeric(as.character(dat.res$Q3))    
    dat.res$RT_detected <- as.numeric(as.character(dat.res$RT_detected))
#     dat.res$score <- as.numeric(as.character(dat.res$score))
    dat.res$isotype <- as.factor(dat.res$isotype) 
    dat.res$iRT <- as.numeric(as.character(dat.res$iRT))
    dat.res$relative_intensity <- as.numeric(as.character(dat.res$relative_intensity))
    dat.res$decoy <- as.factor(dat.res$decoy)
    dat.res$shared <- as.factor(dat.res$shared)
    dat.res$confidence <- as.numeric(as.character(dat.res$confidence))

  } else if(format == "openswath"){ 
    
    Q1 <- as.numeric(as.character(dat$PrecursorMz))
    Q3 <- as.numeric(as.character(dat$ProductMz))
    RT_detected <- as.numeric(as.character(dat$Tr_recalibrated))
    protein_name <- dat$ProteinName
    isotype <- as.factor(dat$LabelType)
    if(length(isotype) == 0) {
      isotype <- as.factor(dat$GroupLabel)
    }
    iRT <- as.numeric(as.character(dat$Tr_recalibrated))
    relative_intensity <- as.numeric(as.character(dat$LibraryIntensity))
    stripped_sequence <- dat$PeptideSequence
    modification_sequence <- dat$FullUniModPeptideName
    uniprot_id <- dat$UniprotID
    decoy <- as.factor(dat$decoy)
    shared <- c(rep(FALSE, nrow(dat)))
    confidence <- c(rep(1, nrow(dat)))
    prec_z <- as.numeric(as.character(dat$PrecursorCharge))
    frg_type <- dat$FragmentType
    frg_z <- as.numeric(as.character(dat$FragmentCharge))
    frg_nr <- as.numeric(as.character(dat$FragmentSeriesNumber))
    N <- c(rep(1, nrow(dat)))
    
    
    dat.res <- data.frame(Q1,Q3,RT_detected,protein_name,isotype,
                          iRT,relative_intensity,stripped_sequence,
                          modification_sequence,uniprot_id,decoy,shared,
                          confidence,prec_z,frg_type,frg_z,
                          frg_nr,N)
    
    if(is.factor(dat.res$protein_name))
      dat.res$protein_name <- as.character(dat.res$protein_name)
    if(is.factor(dat.res$stripped_sequence))
      dat.res$stripped_sequence <- as.character(dat.res$stripped_sequence)
    if(is.factor(dat.res$modification_sequence))
      dat.res$modification_sequence <- as.character(dat.res$modification_sequence)
    if(is.factor(dat.res$uniprot_id))
      dat.res$uniprot_id <- as.character(dat.res$uniprot_id)
    if(is.factor(dat.res$frg_type))
      dat.res$frg_type <- as.character(dat.res$frg_type)
  
  } else if (format == "skyline") {
	
	Q1 <- as.numeric(as.character(dat$PrecursorMz))
    Q3 <- as.numeric(as.character(dat$ProductMz))
    RT_detected <- as.numeric(as.character(dat$Tr_recalibrated))
    protein_name <- dat$ProteinName
    iRT <- as.numeric(as.character(dat$Tr_recalibrated))
    relative_intensity <- as.numeric(as.character(dat$LibraryIntensity))
    stripped_sequence <- dat$PeptideSequence
    modification_sequence <- dat$ModificationSequence
    uniprot_id <- as.character(dat$UniprotID)
    prec_z <- as.numeric(as.character(dat$PrecursorCharge))
    frg_type <- as.character(dat$FragmentType)
    frg_z <- as.numeric(as.character(dat$FragmentCharge))
    frg_nr <- as.numeric(as.character(dat$FragmentSeriesNumber))
    decoy <- as.factor(dat$decoy)
	shared <- c(rep(FALSE, nrow(dat)))
    confidence <- c(rep(1, nrow(dat)))
	N <- c(rep(1, nrow(dat)))
	
	dat.res <- data.frame(Q1,Q3,RT_detected,protein_name,
                          iRT,relative_intensity,stripped_sequence,
                          modification_sequence,uniprot_id, prec_z, frg_type,
	                        frg_z, frg_nr, decoy, shared, confidence, N)
						  
	if(is.factor(dat.res$protein_name))
      dat.res$protein_name <- as.character(dat.res$protein_name)
    if(is.factor(dat.res$stripped_sequence))
      dat.res$stripped_sequence <- as.character(dat.res$stripped_sequence)
    if(is.factor(dat.res$modification_sequence))
      dat.res$modification_sequence <- as.character(dat.res$modification_sequence)
	if(is.factor(dat.res$uniprot_id))
	  dat.res$uniprot_id <- as.character(dat.res$uniprot_id)
	if(is.factor(dat.res$frg_type))
	  dat.res$frg_type <- as.character(dat.res$frg_type)
  }
  
  else { # spectronaut
  
	Q1 <- as.numeric(as.character(dat$PrecursorMz))
    Q3 <- as.numeric(as.character(dat$FragmentMz))
    RT_detected <- as.numeric(as.character(dat$iRT))
    protein_name <- dat$ProteinGroups
    iRT <- as.numeric(as.character(dat$iRT))
    relative_intensity <- as.numeric(as.character(dat$RelativeIntensity))
    stripped_sequence <- dat$StrippedPeptide
    modification_sequence <- dat$ModifiedPeptide
    uniprot_id <- as.character(dat$UniProtIds)
    shared <- c(rep(FALSE, nrow(dat)))
    confidence <- c(rep(1, nrow(dat)))
    prec_z <- as.numeric(as.character(dat$PrecursorCharge))
    frg_type <- as.character(dat$FragmentType)
    frg_z <- as.numeric(as.character(dat$FragmentCharge))
    frg_nr <- as.numeric(as.character(dat$FragmentNumber))
    N <- c(rep(1, nrow(dat)))
	
	dat.res <- data.frame(Q1,Q3,RT_detected,protein_name,
                          iRT,relative_intensity,stripped_sequence,
                          modification_sequence,uniprot_id, shared,
                          confidence,prec_z,frg_type,frg_z,
                          frg_nr,N)
						  
	if(is.factor(dat.res$protein_name))
      dat.res$protein_name <- as.character(dat.res$protein_name)
    if(is.factor(dat.res$stripped_sequence))
      dat.res$stripped_sequence <- as.character(dat.res$stripped_sequence)
    if(is.factor(dat.res$modification_sequence))
      dat.res$modification_sequence <- as.character(dat.res$modification_sequence)
    if(is.factor(dat.res$frg_type))
      dat.res$frg_type <- as.character(dat.res$frg_type)
	if(is.factor(dat.res$uniprot_id))
	  dat.res$uniprot_id <- as.character(dat.res$uniprot_id)
	if(is.factor(dat.res$frg_type))
	  dat.res$frg_type <- as.character(dat.res$frg_type)
  
  }
  
dat.res


}