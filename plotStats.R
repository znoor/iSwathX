###########################################################################
#' Plot the statistical summary of the base and external libraries
#' @param datBaseLib a data frame for base library
#' @param datExtLib a data frame for external/addon library
#' @param ... arguments to pass onto cleanLib function
#' @return a list of ggplot objects of statistical plots of the two libraries
#' @examples
#' libfiles <- paste(system.file("files",package="iSwathX"),c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1])
#' datExtLib <- readLibFile(libfiles[2])
#' list.statplots <- plotStats(datBaseLib, datExtLib)
############################################################################

# install.packages("ggpubr")
library("ggpubr")
# install.packages("VennDiagram")
library("VennDiagram")

source("readLibFile.R")
source("proteinOverlap.R")


plotStats <- function(datBaseLib, datExtLib, wb=NULL, sheet=NULL, ...)
{
  
  list.plots = list()
 
  if(!missing(datExtLib)){ ## 2 libs
    
    ## protein peptide numbers
    nums.pro <- sapply(list(datBaseLib$uniprot_id,datExtLib$uniprot_id),
                       FUN=function(x){length(unique(x))})
   
    nums.pep <- sapply(list(datBaseLib$stripped_sequence,datExtLib$stripped_sequence),
                       FUN=function(x){length(unique(x))})
  
    names <- c("Base", "Addon")
    provalues <- c(Base = nums.pro[[1]], Addon = nums.pro[[2]])
    pepvalues <- c(Base = nums.pep[[1]], Addon = nums.pep[[2]])
    
    npro <- data.frame(names, provalues)
    npep <- data.frame(names, pepvalues)
    
  # png("protein_peptide_numbers.png")
    
  
  p1 <- ggplot(data = npro, aes(names, provalues)) +
    geom_bar(stat = "identity", fill =c("#0080ff", "#fb1f77") , position = "dodge", alpha = 0.8) +
    ggtitle("Proteins") + labs(x="", y="") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face= "bold", size = 10)) 
  
  p2 <- ggplot(data = npep, aes(names, pepvalues)) +
    geom_bar(stat = "identity", fill =c("#0080ff", "#fb1f77") , position = "dodge", alpha = 0.8) +
    ggtitle("Peptides") + labs(x="", y="") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face= "bold", size = 10)) 
  
  ppnum <- ggarrange(p1, p2, ncol = 2)
  filepath <- getwd()
  filepath2 <- paste0(filepath,"/graphs")
  
  ggsave("Protein_Peptide_Numbers.png", width = 6, height = 5, path = filepath2)
  
  # dev.off()
    
    #////////////////////////////////////////////////////////////////////////////////////////////////////////
    
      ## VENN DIAGRAMS

    #grid.newpage()


    if(nums.pro[1]>0 && nums.pro[2]>0){

  #     filepath <- getwd()
  # filepath2 <- paste0(filepath,"/graphs/proteinVennDigram.png")
  # # png(filepath2)
      # png("proteinVennDigram.png")

  prot_base <- datBaseLib[!duplicated(datBaseLib$uniprot_id), "uniprot_id"]
  pep_base <- datBaseLib[!duplicated(datBaseLib$stripped_sequence), "stripped_sequence"]
  
  prot_ext <- datExtLib[!duplicated(datExtLib$uniprot_id), "uniprot_id"]
  pep_ext <- datExtLib[!duplicated(datExtLib$stripped_sequence), "stripped_sequence"]
  
 pprotvenn <- venn.diagram(list(prot_base, prot_ext),
               category.names = c("Base Library Proteins", "External Library Proteins"),
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
 filepath <- getwd()
 filepath2 <- paste0(filepath,"/graphs")
 ggsave("ProteinVennDiagram.png", width = 6, height = 5, path = filepath2)
      # draw.pairwise.venn(nums.pro[1], nums.pro[2],
      #                    length(proteinOverlap(datBaseLib,datExtLib, ...)),
      #                    category = c("Base Library Proteins", "External Library Proteins"),
      #                    lty = rep("blank", 2),
      #                    fill = c("light blue", "pink"),
      #                    alpha = rep(0.5, 2), cat.pos = c(0,20),
      #                    cat.dist = rep(0.025, 2))

 
 # library(grDevices)
 # 
 # png(file=filepath2)
 # grid.draw(temp1)
 #      dev.off()

    }

    if(nums.pep[1]>0 && nums.pep[2]>0){

      
      # filepath <- getwd()
      # filepath2 <- paste0(filepath,"/graphs/peptideVennDigram.png")
      # # png(filepath2)
      # png("peptideVennDigram.png")
      
      ppepvenn <- venn.diagram(list(pep_base, pep_ext),
                   category.names = c("Base Library Peptides", "External Library Peptides"),
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
    }
  
  
  
  ppepgg <- as_ggplot(ppep)
  filepath <- getwd()
  filepath2 <- paste0(filepath,"/graphs")
  ggsave("PeptideVennDiagram.png", width = 6, height = 5, path = filepath2)
  
      # draw.pairwise.venn(nums.pep[1], nums.pep[2],
      #                    length(unique(intersect(datBaseLib$stripped_sequence,
      #                                            datExtLib$stripped_sequence))),
      #                    category = c("Base Library Peptides", "External Library Peptides"),
      #                    lty = rep("blank", 2),
      #                    fill = c("light blue", "pink"),
      #                    alpha = rep(0.5, 2), cat.pos = c(0,20),
      #                    cat.dist = rep(0.025, 2))

      # png(file=filepath2)
      # grid.draw(temp2)
      # dev.off()
      
      ####/////////////////////////////////////////////////////////////////////////////////
 
    # Density plots
    
    
    # png("densityPlots.png")
    
    p1 <- ggplot(data = datBaseLib, aes(datBaseLib$confidence)) +
      stat_density(position="identity",geom="line", color = "#0080ff") +
      ggtitle("Base lib conf") + labs(x= "Confidence", y = "Density") +
       geom_vline(xintercept=0.99,col="red") +
      scale_x_continuous(breaks = seq(format(min(datBaseLib$confidence), digits = 2), 1, 0.4))+
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p2 <- ggplot(data = datExtLib, aes(datExtLib$confidence)) +
      stat_density(position="identity",geom="line", color = "#0080ff") +
      ggtitle("Addon lib conf") + labs(x= "Confidence", y = "Density") +
      geom_vline(xintercept=0.99,col="red") +
      scale_x_continuous(breaks = seq(format(min(datExtLib$confidence), digits = 2), 1, 0.4))+
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p3 <- ggplot(data = datBaseLib, aes(log(datBaseLib$relative_intensity))) +
      stat_density(position="identity",geom="line", color = "#fb1f77") +
      ggtitle("Base lib ion intensity") + labs(x = "log ion intensity", y = "Density") +
      geom_vline(xintercept=1.6,col="red") +
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p4 <- ggplot(data = datExtLib, aes(log(datExtLib$relative_intensity))) +
      stat_density(position="identity",geom="line", color = "#fb1f77") +
      ggtitle("Addon lib ion intensity") + labs(x = "log ion intensity", y = "Density") +
      geom_vline(xintercept=1.6,col="red") +
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    pdens <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    
    filepath <- getwd()
    filepath2 <- paste0(filepath,"/graphs")
    ggsave("Density_Plots.png", width = 6, height = 5, path = filepath2)
    # dev.off()
  ###/////////////////////////////////////////////////////////////////////////////  
    
    # histgrams of peptides/protein and ions/peptide
    
    # png("histogram_peptide_ion.png")
    
    # x<-aggregate(datBaseLib$uniprot_id, by=list(Protein=datBaseLib$uniprot_id, Peptide=datBaseLib$stripped_sequence), 
    #              FUN=function(x){length(x)})
    # y<-aggregate(x$Peptide, by=list(Protein=x$Protein), FUN=function(x){length(x)})
    
    datBaseLib2 <- datBaseLib[!duplicated(datBaseLib$modification_sequence),]
    x <- aggregate(x = datBaseLib2$modification_sequence, by = list(Protein = datBaseLib2$uniprot_id, charge = datBaseLib2$prec_z),
                   FUN = function(x) {length(x)})
    x$charge <- as.character(x$charge)
    breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
    bwidth1 <- breaks1[3] - breaks1[1]
    
    x2 <- aggregate(x = datBaseLib$modification_sequence, by = list(Peptide = datBaseLib$modification_sequence, 
                                                                  frg_type = datBaseLib$frg_type), FUN = function(x2) {length(x2)})
    breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
    bwidth2 <- breaks2[2] - breaks2[1]
    
    p1 <- ggplot(data = x, aes(x = x, fill = charge)) +
      geom_histogram(binwidth = bwidth1) +
      ggtitle("Base lib") +
      labs(x ="Number of Peptides", y = "Number of Proteins") +
      # scale_y_continuous(name = "Frequency") +
      scale_fill_brewer("Precursor Charge", palette = "YlOrRd") +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10))  
    
    p2 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
      geom_histogram(binwidth = bwidth2) +
      labs(x ="Number of Ions", y = "Number of Peptides") +
      ggtitle("Base lib") +
      scale_fill_brewer("Fragment Type", palette = "PuRd") +
      # scale_fill_manual(values = c("#eeba30", "#ae0001")) +
      # scale_y_continuous(name = "Frequency") +
      # scale_fill_gradient("Frequency", low = "blue", high = "red")+
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10))  
    
    #################################################################################################
    
    # x<-aggregate(datExtLib$uniprot_id, by=list(Protein=datExtLib$uniprot_id, Peptide=datExtLib$stripped_sequence), 
    #              FUN=function(x){length(x)})
    # y<-aggregate(x$Peptide, by=list(Protein=x$Protein), FUN=function(x){length(x)})
    
    datExtLib2 <- datExtLib[!duplicated(datExtLib$modification_sequence),]
    x <- aggregate(x = datExtLib2$modification_sequence, by = list(Protein = datExtLib2$uniprot_id, charge = datExtLib2$prec_z),
                   FUN = function(x) {length(x)})
    x$charge <- as.character(x$charge)
    breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
    bwidth1 <- breaks1[3] - breaks1[1]
    
    x2 <- aggregate(x = datExtLib$modification_sequence, by = list(Peptide = datExtLib$modification_sequence, 
                                                                  frg_type = datExtLib$frg_type), FUN = function(x2) {length(x2)})
    breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
    bwidth2 <- breaks2[2] - breaks2[1]
    
    p3 <- ggplot(data = x, aes(x = x, fill = charge)) +
      geom_histogram(binwidth = bwidth1) +
      labs(x ="Number of Peptides", y = "Number of Proteins") +
      ggtitle("Addon lib") +
      # scale_y_continuous(name = "Frequency") +
      scale_fill_brewer("Precursor Charge", palette = "YlOrRd") +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10))  
    
    p4 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
      geom_histogram(binwidth = bwidth2) +
      labs(x ="Number of Ions", y = "Number of Peptides") +
      ggtitle("Addon lib") +
      scale_fill_brewer("Fragment type", palette = "PuRd") +
      # scale_y_continuous(name = "Frequency") +
      # scale_fill_gradient("", low = "blue", high = "red")+
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10))
    
    phist <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    
    filepath <- getwd()
    filepath2 <- paste0(filepath,"/graphs")
    ggsave("Histogram_Peptide_Ion.png", width = 8, height = 5, path = filepath2)
    
    # dev.off()
    
 ####///////////////////////////////////////////////////////////////////////////////////
    
  } else { # single lib
    nums.pro <- length(unique(datBaseLib$uniprot_id)) 
    
    nums.pep <- length(unique(datBaseLib$stripped_sequence))
    
    names <- "Base"
    provalues <- c(Base = nums.pro[[1]])
    pepvalues <- c(Base = nums.pep[[1]])
    
    npro <- data.frame(names, provalues)
    npep <- data.frame(names, pepvalues)
    
    
    # png("protein_peptide_numbersggplot.png")
    
    p1 <- ggplot(data = npro, aes(names, provalues)) +
      geom_bar(stat = "identity", fill ="#0080ff" , position = "dodge", alpha = 0.8) +
      ggtitle("Proteins") + labs(x="", y="") +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p2 <- ggplot(data = npep, aes(names, pepvalues)) +
      geom_bar(stat = "identity", fill ="lightseagreen" , position = "dodge", alpha = 0.8) +
      ggtitle("Peptides") + labs(x="", y="") +
      theme_classic() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    p <- ggarrange(p1, p2, ncol = 2)
   
    filepath <- getwd()
    filepath2 <- paste0(filepath,"/graphs")
     ggsave("Protein_Peptide_Numbers.png", width = 6, height = 5, path = filepath2)
    # dev.off()
    
###//////////////////////////////////////////////////
    
    # png("densityplotsggplot.png")
    
    p1 <- ggplot(data = datBaseLib, aes(datBaseLib$confidence)) +
      stat_density(position="identity",geom="line", color = "#0080ff") +
      ggtitle("Base lib conf") + labs(x= "Confidence", y = "Density") +
      geom_vline(xintercept=0.99,col="red") +
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p2 <- ggplot(data = datBaseLib, aes(log(datBaseLib$relative_intensity))) +
      stat_density(position="identity",geom="line", color = "#fb1f77") +
      ggtitle("Ion intensity") + labs(x = "log ion intensity", y = "Density") +
      geom_vline(xintercept=1.6,col="red") +
      theme_classic2() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(face = "bold", size = 10),
            axis.text.y = element_text(face= "bold", size = 10)) 
    
    p <- ggarrange(p1, p2, ncol = 2)
    
    filepath <- getwd()
    filepath2 <- paste0(filepath,"/graphs")
    ggsave("Density_Plots.png", width = 6, height = 5, path = filepath2)
    # dev.off()
    
    
    
    
  }
  
  list.plots[["ppnum"]] = ppnum
  list.plots[["pdens"]] = pdens
  list.plots[["phist"]] = phist
  list.plots[["pvennprot"]] = pprotgg
  list.plots[["pvennpep"]] = ppepgg
  
  return(list.plots)
}
