###########################################################################
#' Check if a sequence has miss cleavages
#' @param seq A character string representing the peptide sequence
#' @param enz A character string representing the enzyme which can be one of
#' "trypsin" (defalut), "gluc", or "chymotrypsin"
#' @return a boolean value
#' @examples 
#' isMissCleavaged("VSQGSKDPAEGDGAQPEETPR")
############################################################################ 

isMissCleavaged <- function(seq, enz = c("trypsin", "gluc", "chymotrypsin") )
{
  
  enz <- match.arg(enz)
  if(enz == "trypsin") ## Trypsin rule: K or R not before P  
	  grepl("^\\w*[KR][^P]\\w*",seq)
  else if (enz == "gluC") ## GluC rule: D or E not before P
	  grepl("^\\w*[DE][^P]\\w*",seq)
  else if (enz == "chymotrypsin") ## Chymotrypsin rule: F/Y/W/M/L, not before P
	  grepl("^\\w*[FYWML][^P]\\w*",seq)
  else
    warning("Unknow enzyme type. Must be one of trypsin, gluc, or chymotrypsin")
}