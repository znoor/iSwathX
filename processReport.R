
###########################################################################
#' Processing Report File 
#' @param dat a data frame for a report file
############################################################################                      

library("dplyr")
library("stringr")

processReport <- function(dat, spiked = "sp|", avgDotP = TRUE,
                          avgMass = TRUE, avgRT = TRUE, avgQval = TRUE
                          # , allRep = 3
                          # , bioRep = 3
                          )
  
{
   if(spiked != ""){
    dat <- dat %>%
      filter(!(dat$Protein.Name %in% str_subset(string = dat$Protein.Name, pattern = fixed(spiked))))
  }
  
  if(avgDotP) {
    
    dat[, str_which(colnames(dat), pattern = "Library.Dot.Product")] <- apply(m <- (dat[, str_which(colnames(dat), pattern = "Library.Dot.Product")]), 2, as.numeric)
    
    dat <- dat %>%
      mutate(Average.Library.Dot.Product = apply(m <- (dat[, str_which(colnames(dat), pattern = "Library.Dot.Product")]), 1, mean, na.rm = T))
    
    dat$Average.Library.Dot.Product <- round(dat$Average.Library.Dot.Product, digits = 2)
    # dat[dat == 0] <- NA
  }
   
  if(avgMass) { 
    
    dat[, str_which(colnames(dat), pattern = "Average.Mass.Error.PPM")] <- apply(m <- (dat[, str_which(colnames(dat), pattern = "Average.Mass.Error.PPM")]), 2, as.numeric)
    
    dat <- dat %>%
      mutate(Average.Mass.Error.PPM = apply(m <- (dat[, str_which(colnames(dat), pattern = "Average.Mass.Error.PPM")]), 1, mean, na.rm = T))
  }
  
  if(avgRT) {
    
    dat[, str_which(colnames(dat), pattern = "Peptide.Retention.Time")] <- apply(m <- (dat[, str_which(colnames(dat), pattern = "Peptide.Retention.Time")]), 2, as.numeric)
    
    dat <- dat %>%
      mutate(Average.Measured.Retention.Time = apply(m <- (dat[, str_which(colnames(dat), pattern = "Peptide.Retention.Time")]), 1, mean, na.rm = T))
  }
  
  if(avgQval) {
    
    dat[, str_which(colnames(dat), pattern = "Detection.Q.Value")] <- apply(m <- (dat[, str_which(colnames(dat), pattern = "Detection.Q.Value")]), 2, as.numeric)
    
    dat <- dat %>%
      mutate(Average.Annotation_QValue = apply(m <- (dat[, str_which(colnames(dat), pattern = "Detection.Q.Value")]), 1, mean, na.rm = T))
  }
  
  dat
}




