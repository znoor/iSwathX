###########################################################################
#' Plot for ion intensity ranking correlation of two libraries
#' @param dat1 A data frame containing the first spectrum library
#' @param dat2 A data frame containing the second spectrum library
#' @param method a character string indicating the method for calculating
#'  correlation coefficient. Either of "spearman" (default) or "kendall".
#' @return a list of boxplot and correlation dataframe for the ion intensity ranking correlation
#' between the two libraries
#' @examples
#' file1 <- paste(system.file("files",package="SwathXtend"),"Lib2.txt",sep="/")
#' file2 <- paste(system.file("files",package="SwathXtend"),"Lib3.txt",sep="/")
#' dat1 <- normalise(readLibFile(file1))
#' dat2 <- normalise(readLibFile(file2))
#' list.intensitycor <- computeIntensityCor(dat1, dat2)
############################################################################

#### different from swathxtend
#install.packages("pryr")
library("pryr")
library("graphics")

library("ggplot2")
library("ggthemes")
library("grid")

source("getRTrain.R")
source("normalise.R")


computeIntensityCor <- function(dat1,dat2, method = "spearman", wb=NULL, sheet=NULL, nomod=FALSE)
  
{
  
  datRTrain <- getRTrain(dat1, dat2, nomod)
  if(nrow(datRTrain) < 20) {
    stop("Size of RT training set is too small (<20) to get a good model!")
  } else {
  
  comm.peps <- unique(datRTrain$sequence)
  
  #   comm.peps <- unique(intersect(dat1$stripped_sequence,
  #                                 dat2$stripped_sequence))
  
  cols <- c("stripped_sequence", "relative_intensity", "frg_type", "frg_nr",            
            "frg_z" )
  
  if(nomod) { 
    
    
    datlibs.sub = lapply(list(dat1, dat2),
                         function(x){x[x$stripped_sequence %in% 
                                         comm.peps,cols]})
    
    datlibs.sub = lapply(datlibs.sub, 
                         function(x) {x$frg_type_nr=
                           paste(x$stripped_sequence,x$frg_type,x$frg_nr,x$frg_z)
                         x})
  } else { # with modification
    
    
    datlibs.sub = lapply(list(dat1, dat2),
                         function(x){
                           x$stripped_sequence <- paste(x$stripped_sequence,x$modification_sequence)
                           x[x$stripped_sequence %in% comm.peps,cols]})
    
    datlibs.sub = lapply(datlibs.sub, 
                         function(x) {x$frg_type_nr=
                           paste(x$stripped_sequence,x$frg_type,x$frg_nr,x$frg_z)
                         x})
  }
  
  
  mm = merge(datlibs.sub[[1]][,c(1,2,6)],datlibs.sub[[2]][,c(1,2,6)],by=3,all=FALSE)
  if(nrow(mm) > 5){
    mm = mm[!duplicated(mm$frg_type_nr),c(1,2,3,5)]
    
    colnames(mm) <- c("frg_type_nr","stripped_sequence","relative_intensity.x",
                      "relative_intensity.y")
    ## correlation statistics
    
    list.ms2 <- split(mm,as.factor(mm$stripped_sequence))
    
    
    if(length(list.ms2) > 1){
      allCor<-unlist(sapply(list.ms2, FUN=function(x)
      {if(nrow(x)>1 & length(unique(x[,3]))>1 & length(unique(x[,4]))>1) {cor(x[,3],x[,4],method=method)}}) ) 
      
      ionCorGS<-NULL; rm(ionCorGS);
      data(ionCorGS)
      
      
      ioncordata <- data.frame(ionCorGS)
       allCordata <- data.frame(allCor)
    
       # png("IonIntensityCorrelationggplot.png")
       plotg <- ggplot(data = allCordata, aes(x= "", y = allCordata)) +
         geom_boxplot(outlier.colour = "red", alpha = 0.7, outlier.shape = 8, outlier.size = 2, fill = "#146d14", colour = "darkgreen")+
         labs(x = "In Study", y = "")+
         ggtitle ("Relative ion intensity correlation between libraries")+
         scale_y_continuous(limits=c(min(allCordata),1), breaks = seq(format(min(allCordata), digits = 2), 1, 0.2))+
         geom_hline(yintercept = 0.8, color = 'red')+
       theme_classic()+
         theme(plot.title = element_text(size = 12, face = "bold"),
               text = element_text(size = 12),
               axis.title = element_text(size = 11, face = "bold"),
               axis.text.x = element_text(face = "bold", size = 10),
               axis.text.y = element_text(face = "bold", size = 10))
       filepath <- getwd()
       filepath2 <- paste0(filepath,"/graphs")
       ggsave("Ion_Intensity_Correlation.png", width = 6, height = 4, path = filepath2)
       # dev.off()
    }
    
   
  }
  }
  
  list(allCor, mm, plotg)
}