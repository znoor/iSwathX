###########################################################################
#' Plot for hydrophobicity based retention time correlation of two libraries
#' @param dat1 A data frame containing the first spectrum library
#' @param dat2 A data frame containing the second spectrum library
#' @param label1 a character string representing the x axis label for plotting
#' @param label2 a character string representing the y axis label for plotting
#' @param nomod a logic value, representing if the modified peptides and its
#' fragment ions will be removed. FALSE (default) means not removing.
#' @return Plots of RT correlation and RT residual will be saved
#' under current working directory. a list of training set dataframe, hydro-based plot and hydro-based residual plot.
#' @examples
#' libfiles <- paste(system.file("files",package="iSwathX"),
#'    c("Lib2.txt","Lib3.txt"),sep="/")
#' hydroFile <- paste(system.file("files", package = "iSwathX"),
#'                                       "hydroIndex.txt", sep = "/")
#' datBaseLib <- readLibFile(libfiles[1])
#' datExtLib <- readLibFile(libfiles[2])
#' hydro <- readLibFile(hydroFile, type = "hydro")
#'
#' list.hydrocor <- computeRTHydro(datBaseLib, datExtLib, hydro)
#' list.hydrocor[[2]]
#' list.hydrocor[[3]]
############################################################################

computeRTHydro <- function(dat1, dat2, datHydroIndex, nomod=FALSE, label1="Base Lib", label2 = "External Lib")

{

  if(missing(datHydroIndex)) stop("Missing data frame of hydrophibicity index.")

  colnames(datHydroIndex)[grep("sequence",tolower(colnames(datHydroIndex)))]<-
    "sequence"

  colnames(datHydroIndex)[grep("hydro",tolower(colnames(datHydroIndex)))] <- "hydro"

  colnames(datHydroIndex)[grep("len",tolower(colnames(datHydroIndex)))] <- "length"

  dat1 <- dat1[!duplicated(dat1$stripped_sequence),]

  datHydroIndex <- datHydroIndex[!duplicated(datHydroIndex$sequence),]

  datHydroIndex$hydro <- as.numeric(datHydroIndex$hydro)

  datRTrain <- merge(dat1[,c("stripped_sequence","RT_detected")],
                     datHydroIndex[,c("sequence","hydro","length")],
                     by.x="stripped_sequence",
                     by.y=tolower("sequence"),all=FALSE)


  datRTrain$length <- log(datRTrain$length)

  f.h= as.formula("RT_detected~hydro")

  f.hl <- as.formula("RT_detected~hydro+length")


  if(nrow(datRTrain) > 10){
    fit.h <- lm(f.h, data=datRTrain)
    fit.hl <- lm(f.hl, data=datRTrain)

    #### RT correlation hydro plot ///////////////////////////////////////////////

    # png("RT_correlation_hydroggplot.png")

    r2 <- cor(datRTrain[,3],datRTrain[,2])^2
    r2label <- bquote(italic(R)^2 ==. (format(r2, digits = 2)))

    rtcorhydro <- ggplot(data = datRTrain, aes(x = datRTrain[,2], y = datRTrain[,3])) +
      geom_point(shape=21, fill='#4c4c4c', color = '#191919', size = 2, alpha = 0.5) +

      annotate(geom="text", x=min(datRTrain[,2])+4, y=max(datRTrain[,3])-4, label=deparse(r2label), parse = T,
               color="red")+
      geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +

      scale_x_continuous(name="Hydrophibicity index", limits=c(min(datRTrain[,2]), max(datRTrain[,2])))+
      scale_y_continuous(name="Base library RT", limits=c(min(datRTrain[,3]), max(datRTrain[,3])) )+
      scale_color_brewer(palette="Accent")+
      ggtitle("Correlation between Retention time and Hydrophobicity index")+
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face = "bold", size = 10))


    ggsave("RT_Correlation_Hydro.png", width = 6, height = 5)
    # dev.off()


    #### RT correlation hydro length plot ///////////////////////////////////////////////


    # png("RT_correlation_hydroNlengthggplot.png")


    r2 <- cor(datRTrain[,3],datRTrain[,2])^2
    r2label <- bquote(italic(R)^2 ==. (format(r2, digits = 2)))

    rtcorhydrolen <- ggplot(data = datRTrain, aes(x = datRTrain[,2], y = datRTrain[,3])) +
      geom_point(shape=21, fill='#4c4c4c', color = '#191919', size = 2, alpha = 0.5) +

      annotate(geom="text", x=min(datRTrain[,2])+4, y=max(datRTrain[,3])-4, label=deparse(r2label), parse = T,
               color="red")+
      geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +

      scale_x_continuous(name="Base Lib Hydro", limits=c(min(datRTrain[,2]), max(datRTrain[,2])))+
      scale_y_continuous(name="Base Lib RT", limits=c(min(datRTrain[,3]), max(datRTrain[,3])) )+
      scale_color_brewer(palette="Accent")+
      ggtitle("Correlation between Retention time and Hydrophobicity index")+
      theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10))


    r2 <- cor(datRTrain[,4],datRTrain[,2])^2
    r2label <- bquote(italic(R)^2 ==. (format(r2, digits = 2)))

    rtcorrtlen <- ggplot(data = datRTrain, aes(x = datRTrain[,4], y= datRTrain[,2])) +
      geom_point(shape=21, fill='#4c4c4c', color = '#191919', size = 2, alpha = 0.5) +

      annotate(geom="text", x=min(datRTrain[,4])+4, y=max(datRTrain[,2])-4, label=deparse(r2label), parse = T,
               color="red")+
      geom_smooth(method="lm", color="#ff0000", se = F, fullrange = T) +
      scale_x_continuous(name="ln(Length)", limits=c(min(datRTrain[,4]), max(datRTrain[,4])))+
      scale_y_continuous(name="RT", limits=c(min(datRTrain[,2]), max(datRTrain[,2])) )+

      scale_color_brewer(palette="Accent")+
      ggtitle("Correlation between ln(Length) and RT")+
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            text = element_text(size = 12),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face = "bold", size = 10))


    rthydro <- ggarrange(rtcorhydrolen, rtcorrtlen, nrow = 2)

    ggsave("RT_Correlation_HydroNLength.png", width = 6, height = 5)

    # dev.off()


    rmse.h <- sqrt(mean(resid(fit.h)^2))
    rmse.hl <- sqrt(mean(resid(fit.hl)^2))


    ###### RT residual hydro plot ///////////////////////////////////////////////

    # png("RT_residual_hydroggplot.png")
    rmselabel <- bquote(italic(Rmse) ==. (format(rmse.h, digits=2)))

    resids <- residuals(fit.h)
    x= seq(1,length(resids))
    data <- data.frame(x, resids)

    reshydro <- ggplot(data = data, aes(x = x, y = resids)) +
      geom_point(shape=21, fill='#4c4c4c', color = '#191919', size = 2, alpha = 0.5)+
      annotate(geom="text", x=min(data$x)+500, y=max(resids)-4, label=deparse(rmselabel), parse = T,
               color="red")+

      scale_x_continuous(name="Index", limits=c(0, length(data$x)))+
      scale_y_continuous(name="Resids", limits=c(min(data$resids), max(data$resids)) )+

      geom_hline(yintercept = floor(rmse.h), color = "red") +
      geom_hline(yintercept = -floor(rmse.h), color = "red")+
      ggtitle("Residual plot for lm (time and hydro)")+
      theme_classic() +
        theme(plot.title = element_text(size = 12, face = "bold"),
              text = element_text(size = 12),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10))

    ggsave("RT_residual_Hydro.png", width = 6, height = 5)

    # dev.off()

    #### RT residual hydro length plot ///////////////////////////////////////////////

    # png("RT_residual_hydrolengthggplot.png")
    rmselabel <- bquote(italic(Rmse) ==. (format(rmse.hl, digits=2)))

    resids <- residuals(fit.hl)
    x= seq(1,length(resids))
    data <- data.frame(x, resids)

    reshydrolen <- ggplot(data = data, aes(x = x, y = resids)) +
      geom_point(shape=21, fill='#4c4c4c', color = '#191919', size = 2, alpha = 0.5)+
      annotate(geom="text", x=min(data$x)+500, y=max(resids)-4, label=deparse(rmselabel), parse = T,
               color="red")+

      scale_x_continuous(name="Index", limits=c(0, length(data$x)))+
      scale_y_continuous(name="Resids", limits=c(min(data$resids), max(data$resids)) )+

      geom_hline(yintercept = floor(rmse.hl), color = "red") +
      geom_hline(yintercept = -floor(rmse.hl), color = "red")+
      ggtitle("Residual plot for lm (time and hydro/length)")+
      theme_classic() +
        theme(plot.title = element_text(size = 12, face = "bold"),
              text = element_text(size = 12),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 10),
              axis.text.y = element_text(face = "bold", size = 10))

    ggsave("RT_residual_HydroNLength.png", width = 6, height = 5)

    # dev.off()



    r2hydro <- cor(datRTrain[,3],datRTrain[,2])^2

    r2len <- cor(datRTrain[,4],datRTrain[,2])^2

    list(r2hydro, rtcorhydro, reshydro, datRTrain)

  }

  }
