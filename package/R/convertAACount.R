
convertAACount<-function(sequence,with.charge=FALSE,with.time=TRUE){
  aminoAcids=c("A","R","N","D","C","E","Q","G","H","I","L",
               "K","M","F","P","S","T","W","Y","V")
  
  aa.count=sapply(aminoAcids, FUN=function(x)
			{
			 if(sum(gregexpr(x,sequence)[[1]]) > 0){
				length(gregexpr(x,sequence)[[1]])
			} else {0} 
			} )
										
  
  data.frame(sequence,t(aa.count))
}



