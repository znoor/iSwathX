###########################################################################
#' Plot for hydrophobicity based retention time correlation of two libraries
#' @param dat1 A data frame containing the first spectrum library
#' @param dat2 A data frame containing the second spectrum library
#' @param label1 a character string representing the x axis label for plotting
#' @param label2 a character string representing the y axis label for plotting
#' @param nomod a logic value, representing if the modified peptides and its
#' fragment ions will be removed. FALSE (default) means not removing.
#' @return two .png files of RT correlation and RT residual plots will be saved
#' under current working directory. a list of training set dataframe, RTcor plot, RT residual plot, and RT correlation value.
#' @examples
#' libfiles <- paste(system.file("files",package="iSwathX"),
#'    c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1])
#' datExtLib <- readLibFile(libfiles[2])
#' list.timecor <- computeRTTime(datBaseLib, datExtLib)
############################################################################
#
#      Retention time correlation function (separate from SwathXtend computeRTCor)

source("getRTrain.R")
source("selectModel.R")
source("convertAACount.R")

library("ggplot2")
library("ggthemes")
library("grid")


computeRTTime <- function(dat1, dat2, label1="Base Lib", label2 = "External Lib", nomod = FALSE)
{
    datRTrain <- getRTrain(dat1, dat2, nomod)
    bestModel <- selectModel(datRTrain)
    r2 <- NULL
    rmse <- NULL
      
      
    if(!is.null(bestModel)){
      
      r2 <- cor(datRTrain[,2],datRTrain[,3])^2
      
      
      predicted <- predict(bestModel, newdata=datRTrain)
      
      rmse <- sqrt(mean((datRTrain[,2]-predicted)^2))      
      
      r2label <- bquote(italic(R)^2 ==. (format(r2, digits = 2)))
      
      resids <- predicted-datRTrain[,2]
      
      rmselabel <- bquote(italic(RMSError) ==. (format(rmse, digits=2)))
      
      # png(file = "RT_correlation_ggplot.png")
      
      rtcor <- ggplot(data = datRTrain, aes(x = datRTrain[,2], y = datRTrain[,3])) +
        geom_point(shape=21, fill="#56B4E9", color = "#56B4E9", size = 2, alpha = 0.5) + 
        
        annotate(geom="text", x=min(datRTrain[,2])+4, y=max(datRTrain[,3])-4, label=deparse(r2label), parse = T,
                 color="red")+
        geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +
        
        scale_x_continuous(name=label1, limits=c(min(datRTrain[,2]), max(datRTrain[,2])))+
        scale_y_continuous(name=label2, limits=c(min(datRTrain[,3]), max(datRTrain[,3])) )+
        scale_color_brewer(palette="Accent")+
        ggtitle("Retention time correlation between base and external library")+
        theme_classic() +
        theme(plot.title = element_text(size = 12, face = "bold"),
              text = element_text(size = 12),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10))
     
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs") 
      ggsave("RT_Correlation.png", width = 6, height = 4, path = filepath2)
      
      # dev.off()
      
      
      # png("RT_residual_ggplot.png")
      
      x= seq(1,length(resids))
      data <- data.frame(x, resids)
      
      rtres <- ggplot(data = data, aes(x = x, y = resids)) +
        geom_point(shape=21, fill="#56B4E9", color = "#56B4E9", size = 2, alpha = 0.5)+
        annotate(geom="text", x=min(data$x)+500, y=max(resids)-4, label=deparse(rmselabel), parse = T,
                 color="red")+
        
        ggtitle(paste("Residual plot for",
                           attr(bestModel,"class")[length(attr(bestModel,"class"))]))+
        scale_x_continuous(name="Index", limits=c(0, length(data$x)))+
        scale_y_continuous(name="Residuals", limits=c(min(data$resids), max(data$resids)) )+
        
        geom_hline(yintercept = floor(rmse), color = "red") +
        geom_hline(yintercept = -floor(rmse), color = "red")+
        theme_classic()+
        theme(plot.title = element_text(size = 12, face = "bold"),
              text = element_text(size = 12),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10))
      
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs")
      ggsave("RT_residual.png", width = 6, height = 4, path = filepath2)
      
      # dev.off()
    
    r2 <- cor(datRTrain[,2],datRTrain[,3])^2
   
      list(rtcor,rtres, datRTrain, r2)

    
  } 
  
}