


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
source("buildSpectraLibPair.R")


buildSpectralLibThree <- function(baseLib,extLib1, extLib2,
                                  method = "time",
                                  cutoff.size = 50,
                                  cutoff.r2 = 0.8,
                                  recalibrate = FALSE,
                                  label1 = NA,
                                  label2 = NA,
                                  label3 = NA,
                                  formatBase = c("peakview", "openswath", "skyline", "spectronaut"),
                                  formatExt1 = c("peakview", "openswath", "skyline", "spectronaut"),
                                  formatExt2 = c("peakview", "openswath", "skyline", "spectronaut"),
                                  outputFormat = c("peakview", "openswath", "skyline", "spectronaut"),
                                  plot = FALSE,
                                  clean = FALSE, 
                                  nomod = FALSE,
                                  merge = FALSE,
                                  parseAcc = FALSE,
                                  consolidateAccession = FALSE, ...)
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
  
  # 
  corRT1 <- cor(datRTrain1$RT_detected, datRTrain1$RT_detected.base)^2
  r2label <- bquote(italic(R)^2 ==. (format(corRT1, digits = 2)))
  comm_pep1 <- nrow(datRTrain1)
  corRI1 <- median(computeIntensityCor(datBaseLib, datExtLib1)[[1]]$allCordata)
  
  
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
  
  ####
  
  corr_stats <- data.frame("Libraries" = c("External Lib 1", "External Lib 2"),
                           "Common.Peptides" = c(comm_pep1, comm_pep2),
                           "RT.Correlation" = c(corRT1, corRT2),
                           "RI.Correlation" = c(corRI1, corRI2))
  
  ### order on the basis of RT cor
  corr_stats<- corr_stats[order(-corr_stats$RT.Correlation), ]
  corr_stats$Rank <- c(1,2)
  
  #####  CHECKING CLOSE
  
  list.ext = list()
  list.ext[[corr_stats[(corr_stats$Libraries == "External Lib 1"), "Rank"]]] = datExtLib1
  list.ext[[corr_stats[(corr_stats$Libraries == "External Lib 2"), "Rank"]]] = datExtLib2
  
  list.venn = list()
  
  for(i in 1:2)
  {
    datExtLib <- list.ext[[i]]
    
    
    ##### splitLib.R
    ##***
    # list.res = list()
    # datBaseLib$pepmod <- paste(datBaseLib$stripped_sequence, datBaseLib$modification_sequence)
    # datExtLib$pepmod <- paste(datExtLib$stripped_sequence, datExtLib$modification_sequence)
    # 
    
    
    prot_base <- datBaseLib[!duplicated(datBaseLib$uniprot_id), 9]
    pep_base <- datBaseLib[!duplicated(datBaseLib$stripped_sequence), 7]

    prot_ext <- datExtLib[!duplicated(datExtLib$uniprot_id), 9]
    pep_ext <- datExtLib[!duplicated(datExtLib$stripped_sequence), 7]

    pprotvenn <- venn.diagram(list(prot_base, prot_ext),
                              category.names = c("Base Library Proteins", "Add-on Library Proteins"),
                              resolution = 500,
                              cat.default.pos = "outer",
                              fill = c("#fb1f77", "#0080ff"),
                              cat.pos = c(-170, 170),
                              alpha = c(0.5, 0.5),
                              cex = 2,
                              cat.fontface = 2,
                              lty =2,
                              # fontfamily ="sans",
                              # filename = filepath2
                              filename = NULL
    )
    pprot <- gTree(children=pprotvenn)
    pprotgg <- as_ggplot(pprot)

    ppepvenn <- venn.diagram(list(pep_base, pep_ext),
                             category.names = c("Base Library Peptides", "Add-on Library Peptides"),
                             resolution = 500,
                             cat.default.pos = "outer",
                             fill = c("#fb1f77", "#0080ff"),
                             cat.pos = c(-170, 170),
                             alpha = c(0.5, 0.5),
                             cex = 2,
                             cat.fontface = 2,
                             lty =2,
                             # fontfamily ="sans",
                             # filename = filepath2
                             filename = NULL
    )
    ppep <- gTree(children=ppepvenn)
    ppepgg <- as_ggplot(ppep)

    venngg <- ggarrange(pprotgg, ppepgg)

    list.venn[[i]] <- venngg
    
    res <- buildSpectraLibPair(datBaseLib, datExtLib, method = "time", cutoff.size = cutoff.size, cutoff.r2 = cutoff.r2, recalibrateRT = recalibrateRT,
                               clean = FALSE, merge = TRUE, consolidateAccession = consolidateAccession, plot = FALSE, parseAcc = FALSE)
    datBaseLib <- res[[1]]
    
  }
  if(plot)
  list(res[[1]], list.venn)
  else
    list(res[[1]])
  
}