



library(plyr)
library(dplyr)
library(e1071)

library("ggpubr")
library("grid")
# install.packages("VennDiagram")
library("VennDiagram")


source("readLibFile.R")
source("normalise.R")
source("splitLib.R")
source("predictRT.R")
# source("alignRTbyHydro.R")
source("libraryFormat.R")
source("consolidateAccession.R")
source("parseAccession.R")
# source("plotStats.R")
source("getRTrain.R")
source("libSummary.R")


multiCorrLibThree <- function(baseLib,extLib1, extLib2,
                              method = "time",
                              cutoff.size = 50,
                              cutoff.r2 = 0.8,
                              formatBase = c("peakview", "openswath", "skyline", "spectronaut"),
                              formatExt1 = c("peakview", "openswath", "skyline", "spectronaut"),
                              formatExt2 = c("peakview", "openswath", "skyline", "spectronaut"),
                              clean = FALSE, 
                              nomod = FALSE,
                              label1 = NA,
                              label2 = NA, 
                              label3 = NA)
{
  #### Reading the library files
  
  if(!is.data.frame(baseLib)){
    datBaseLib <- try(readLibFile(baseLib, format = formatBase, clean=clean, ...))
    
    if(inherits(datBaseLib, "try-error"))
      stop(paste("Error with reading", baseLib))
  } else {
    datBaseLib <- baseLib
  }
  if(!is.data.frame(extLib1)){
    datExtLib1 <- try(readLibFile(extLib1, format = formatExt, clean=clean, ...))
    
    if(inherits(datExtLib1,"try-error"))
      stop(paste("Error with reading", extLib1))
    
  } else {
    datExtLib1 <- extLib1
  }
  if(!is.data.frame(extLib2)){
    datExtLib2 <- try(readLibFile(extLib2, format = formatExt, clean=clean, ...))
    
    if(inherits(datExtLib2,"try-error"))
      stop(paste("Error with reading", extLib2))
    
  } else {
    datExtLib2 <- extLib2
  }
  
  ##### Normalising Relative Intensity Values
  
  datBaseLib <- try(normalise(datBaseLib))
  if(inherits(datBaseLib,"try-error"))
    stop("Error with normalising datBaseLib")	
  
  datExtLib1 <- try(normalise(datExtLib1))
  if(inherits(datExtLib1,"try-error"))
    stop("Error with normalising datExtLib")
  
  datExtLib2 <- try(normalise(datExtLib2))
  if(inherits(datExtLib2,"try-error"))
    stop("Error with normalising datExtLib")
  
  
  #####  CHECKING WHICH EXTERNAL LIBRARY IS BETTER CORRELATED TO THE BASE LIBRARY
  ##### splitLib.R
  ##***
  list.res = list()
  datBaseLib$pepmod <- paste(datBaseLib$stripped_sequence, datBaseLib$modification_sequence)
  datExtLib1$pepmod <- paste(datExtLib1$stripped_sequence, datExtLib1$modification_sequence)
  datExtLib2$pepmod <- paste(datExtLib2$stripped_sequence, datExtLib2$modification_sequence)
  
  
  commpepmod1 <- intersect(datBaseLib$pepmod, datExtLib1$pepmod)
  commpepmod2 <- intersect(datBaseLib$pepmod, datExtLib2$pepmod)
  
  if(length(commpepmod1) == 0)
    warning("The base and first add-on libraries does NOT have any common peptides!\n")
  if(length(commpepmod2) == 0)
    warning("The base and second add-on libraries does NOT have any common peptides!\n")
  
  
  dat.extcomm1 <- datExtLib1[datExtLib1$pepmod %in% commpepmod1,]
  
  dat.extnew1 <- datExtLib1[!datExtLib1$pepmod %in% commpepmod1,]
  
  dat.basecomm1 <- datBaseLib[datBaseLib$pepmod %in% commpepmod1,]    
  
  
  dat.extcomm2 <- datExtLib2[datExtLib2$pepmod %in% commpepmod2,]
  
  dat.extnew2 <- datExtLib2[!datExtLib2$pepmod %in% commpepmod2,]
  
  dat.basecomm2 <- datBaseLib[datBaseLib$pepmod %in% commpepmod2,]
  
  
  list.res[["ExtCommon1"]] = dat.extcomm1
  list.res[["ExtNew1"]] = dat.extnew1
  list.res[["BaseCommon1"]] = dat.basecomm1
  
  list.res[["ExtCommon2"]] = dat.extcomm2
  list.res[["ExtNew2"]] = dat.extnew2
  list.res[["BaseCommon2"]] = dat.basecomm2
  
  list.datLibs <- list.res
  
  ##***
  ##### splitLib.R
  
  datExtLibCommPart1 <- list.datLibs[["ExtCommon1"]]
  datExtLibNewPart1 <- list.datLibs[["ExtNew1"]]
  datBaseCommPart1 <- list.datLibs[["BaseCommon1"]]
  
  datExtLibCommPart2 <- list.datLibs[["ExtCommon2"]]
  datExtLibNewPart2 <- list.datLibs[["ExtNew2"]]
  datBaseCommPart2 <- list.datLibs[["BaseCommon2"]]
  
  
  
  ## RT correlation between base library and external library 1 is 
  datBaseCommPart1 <- libraryFormat(datBaseCommPart1)
  datExtLibCommPart1 <- libraryFormat(datExtLibCommPart1)
  datExtLibNewPart1 <- libraryFormat(datExtLibNewPart1)
  #### Training dataset for RT prediction == training dataset is the common peptides in both libraries for making a regression model
  datRTrain1 <- getRTrain(datBaseCommPart1, datExtLibCommPart1, nomod = nomod) 
  
  
  # label1="Base Lib"
  # label2 = "External Lib"
  # 
  corRT1 <- cor(datRTrain1$RT_detected, datRTrain1$RT_detected.base)^2
  r2label <- bquote(italic(R)^2 ==. (format(corRT1, digits = 2)))
  comm_pep1 <- nrow(datRTrain1)
  corRI1 <- median(computeIntensityCor(datBaseLib, datExtLib1)[[1]]$allCordata)
  
  
  ###### RT plot for seed lib and external library 1
  
  rtcor1 <- ggplot(data = datRTrain1, aes(x = datRTrain1[,2], y = datRTrain1[,3])) +
    geom_point(shape=21, fill="black", color = "gray", size = 2, alpha = 0.5) + 
    
    annotate(geom="text", x=min(datRTrain1[,2])+4, y=max(datRTrain1[,3])-4, label=deparse(r2label), parse = T,
             color="red")+
    geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +
    
    scale_x_continuous(name=label1, limits=c(min(datRTrain1[,2]), max(datRTrain1[,2])))+
    scale_y_continuous(name=label2, limits=c(min(datRTrain1[,3]), max(datRTrain1[,3])) )+
    scale_color_brewer(palette="Accent")+
    ggtitle("Retention time correlation between base and external library")+
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10))
  
  
  ## RT correlation between base library and external library 2 is 
  datBaseCommPart2 <- libraryFormat(datBaseCommPart2)
  datExtLibCommPart2 <- libraryFormat(datExtLibCommPart2)
  datExtLibNewPart2 <- libraryFormat(datExtLibNewPart2)
  #### Training dataset for RT prediction == training dataset is the common peptides in both libraries for making a regression model
  datRTrain2 <- getRTrain(datBaseCommPart2, datExtLibCommPart2, nomod = nomod) 
  
  corRT2 <- cor(datRTrain2$RT_detected, datRTrain2$RT_detected.base)^2
  r2label <- bquote(italic(R)^2 ==. (format(corRT2, digits = 2)))
  comm_pep2 <- nrow(datRTrain2)
  corRI2 <- median(computeIntensityCor(datBaseLib, datExtLib2)[[1]]$allCordata)
  
  rtcor2 <- ggplot(data = datRTrain2, aes(x = datRTrain2[,2], y = datRTrain2[,3])) +
    geom_point(shape=21, fill="black", color = "gray", size = 2, alpha = 0.5) + 
    
    annotate(geom="text", x=min(datRTrain2[,2])+4, y=max(datRTrain2[,3])-4, label=deparse(r2label), parse = T,
             color="red")+
    geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +
    
    scale_x_continuous(name=label1, limits=c(min(datRTrain2[,2]), max(datRTrain2[,2])))+
    scale_y_continuous(name=label3, limits=c(min(datRTrain2[,3]), max(datRTrain2[,3])) )+
    scale_color_brewer(palette="Accent")+
    ggtitle("Retention time correlation between base and external library")+
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10))
  
  corr_stats <- data.frame("Libraries" = c(label2, label3),
                           "Common.Peptides" = c(comm_pep1, comm_pep2),
                           "RT.Correlation" = c(corRT1, corRT2),
                           "RI.Correlation" = c(corRI1, corRI2))
  
  ### order on the basis of RT cor
  corr_stats<- corr_stats[order(-corr_stats$RT.Correlation), ]
  corr_stats$Rank <- c(1,2)
  
  rtcor <- ggarrange(rtcor1, rtcor2, nrow = 1, ncol = 2)
  
  # ### order on the basis of common peptides
  # corr_stats<- corr_stats[order(-corr_stats$Common.Peptides), ]
  # corr_stats$index <- c(1,2)
  
  list(corr_stats, rtcor)
}