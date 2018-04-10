###########################################################################
#' Extract a (list of)  Uniprot accession(s) from a string
#' @param proteinID a character string of protein ID which contains Uniprot
#' accession(s). Two formats are accepted:  1/P0CG40 or sp|O75369|FLNB_HUMAN
#' @param method a character string representing the method, one of "keep"
#' or "replace"
#' @return all Uniprot accessions separated by space
#' @examples 
#' parseAccession("sp|O75369|FLNB_HUMAN")  ## O75369
#' parseAccession("5/Q9ULI2/Q7KZI7/Q9P0L2/P27448/DECOY_Q6ZMI3") #Q9ULI2 Q7KZI7 Q9P0L2 P27448 Q6ZMI3
############################################################################ 
## extract all accession in the uniprot-id
## requires the uniprot_id is in format of 1/P0CG40 or sp|O75369|FLNB_HUMAN

parseAccession <- function(proteinID, method=c("keep","replace"))
{
  ##[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}  
  if(is.na(proteinID) | !grepl("/", proteinID) &
       !grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",proteinID)) {
    proteinID
  } else {
  list.s<-unlist(strsplit(proteinID,split="/"))
  
  method <- match.arg(method)
  paste(na.omit(sapply(list.s, 
         FUN=function(x)
           {
           if(grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",x))
               gsub(".*?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*",
       "\\1",x,perl=TRUE)
          else if (method=="keep") proteinID
           else NA}
         )),collapse=" ")
  
  }
            
}