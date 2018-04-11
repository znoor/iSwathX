###########################################################################
#' Standardise a data frame into library format
#' @param datlib a data frame for a spectrum library
#' @param colnums a number representing the number of columns to keep
#' (default as 18)
#' @return a data frame with specified number of columns
#' @examples
#' libfile <- paste(system.file("files",package="iSwathX"),"Lib2.txt",sep="/")
#' datlib <- readLibFile(libfile)
#' dat <- libraryFormat(datlib)
############################################################################

libraryFormat <- function(datlib)
{

  inclcols = c("Q1",
			"Q3",
			"RT_detected",
			"protein_name",
			"isotype",
			"relative_intensity",
			"stripped_sequence",
			"modification_sequence",
			"prec_z",
			"frg_type",
			"frg_z",
			"frg_nr",
			"iRT",
			"uniprot_id",
			"score",
			"decoy",
			"prec_y",
			"confidence",
			"shared",
			"N",
			"rank",
			"mods",
			"nterm",
			"cterm")

	if( all(inclcols %in% colnames(datlib)) )	{
		res = datlib[,inclcols]
	}
	else {
		miss.col = inclcols[!inclcols %in% colnames(datlib)]

		toadd.col = matrix(NA, ncol=length(miss.col), nrow=nrow(datlib))
		colnames(toadd.col) = miss.col

		res = cbind(datlib, toadd.col)

		# check if all columns are included

		stopifnot(all(inclcols %in% colnames(res)))

		# re-order
		res = res[,inclcols]

	}
	res
}
