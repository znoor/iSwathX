source("readLibFile.R")
# source("coverage.R")
# source("fdrFunctions.R")



# library("actuar", lib.loc="~/R/win-library/3.4")


reliabilityCheckLibrary <- function(seedlib.file, extlib.file)
{

dlib = list()

if(!is.data.frame(seedlib.file)){
  dlib[[1]] = try(readLibFile(seedlib.file))

  if(inherits(dlib[[1]], "try-error"))
    stop(paste("Error with reading", seedlib.file))
} else {
  dlib[[1]] = seedlib.file
}


if(!is.data.frame(extlib.file)){
  dlib[[2]] = try(read.delim(extlib.file, sep='\t', stringsAsFactors=FALSE))

  if(inherits(dlib[[2]],"try-error"))
    stop(paste("Error with reading", extLib.file))

} else {
  dlib[[2]] = extlib.file
}



# pep.col = c('Q1','stripped_sequence','modification_sequence', 'uniprot_id')
lib.labels = c('seed.library', 'extended.library')


lib.stats = matrix(NA, nrow=length(dlib), ncol=2)
rownames(lib.stats) = lib.labels
colnames(lib.stats) = c('proteins','peptides')


for( ii in 1:length(dlib)) {						
	
	x = dlib[[ii]]

	lib.stats[ii,] = c(length(unique(x$uniprot_id)), length(unique(x$modification_sequence)) )
					

}

#####/////////////////////////////////////////////////////////////

# barplots 

filepath <- getwd()
filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
png(filepath2)

# png('plots_of_library_info.png', 2500, 3500, res=300)

sl = ceiling(max(lib.stats[,1])/1000)/ceiling(max(lib.stats[,2])/1000)

pepnum.scaled = lib.stats[,2]*sl

ymax = 1.1*max(c(lib.stats[,1], pepnum.scaled))

bp1 = barplot(lib.stats[,1], space=0.5, names.arg=lib.labels, ylim=c(0, ymax),
	ylab='Number of proteins', axes=FALSE, col=gray(0.8), main='Libraries')

at1 = axis(side=2)

points(bp1, lib.stats[,2]*sl, type='b', lwd=3, col=gray(0.4) )
at2 = round(at1/(sl*1000))*1000
axis(side=4, at = at1, labels=at2)
mtext('Number of peptides', side=4, line=3)

par(xpd=TRUE)
legend(bp1[1,1], ceiling(1.05*ymax) , bty='n', fill=gray(0.8), legend='protein')
legend(bp1[2,1], ceiling(1.05*ymax), bty='n', col=gray(0.4),  lty=1, legend='peptide', lwd=2)

dev.off()

####?//////////////////////////////////////////////////////////////////////////////////////

# lib.stats = rbind(lib.stats, c(coverage(dlib[[1]]$uniprot_id, dlib[[2]]$uniprot_id),
# 					coverage(dlib[[1]]$modification_sequence, dlib[[2]]$modification_sequence) ) )
# 
# rownames(lib.stats)[3] = 'coverage.perc'

# filepath <- getwd()
# filepath2 <- paste0(filepath,"/output_spectral_libraries")

# write.csv(lib.stats, 'table_library_coverage.csv')

lib.stats
}












