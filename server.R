#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=10000*1024^2)
library(shiny)
library(shinyBS)

#source("outputLib.R")
#source("openSwathFormat.R")

source("readLibFile.R")
source("cleanLib.R")
source("OswathFormat.R")
source("peakviewFormat.R")
source("skylineFormat.R")
source("spectronautFormat.R")
source("computeRTTime.R")
source("computeRTHydro.R")
source("computeIntensityCor.R")
source("buildSpectraLibPair.R")
source("reliabilityCheckLibrary.R")
source("libSummary.R")
#source("reliabilityCheckSwath.R")


shinyServer(function(input, output, session) {

  
  shinyjs::disable("apply")
  shinyjs::disable("inclength")
  shinyjs::disable("downloadseedinputLibs")
  shinyjs::disable("downloadextinputLibs")
  # shinyjs::disable("wizdownloadoutputLib")
  # shinyjs::disable("autowizrun")
#### Apply button which is activating the tabs

  observeEvent(input$apply, {
      
      shinyjs::show(id = "libreadpanel")
      showTab(inputId = "libraries", target = "Seed Library", select = TRUE)
    })
  
  
  ##### auto wizard start over button try
  
  ##### Parameters and Log Summary for auto wizard run
  
  observeEvent(input$autowizrun, 
                {
                  
                  
                output$autocombparameters <- renderText({
                  paste0("Intensity cut-off = 5", "\n",
                         "Confidence cut-off = 0.99", "\n",
                         "Training set cut-off = 50", "\n",
                         "Intensity correlation cut-off = 0.8", "\n",
                         "Residual standard error cut-off = 3.0", "\n",
                         "Method = Time", "\n",
                         "Consolidate accessions = TRUE", "\n",
                         "Recalibration = TRUE, (in case of poor correlation)", "\n",
                         "Generate plots = TRUE")
                })
                    
                
                output$autocomblog <- renderText({
                  
                  dat1 <- autolib1read()
                  dat2 <- autolib2read()
                  newdat <- autocombLib()
                  
                  withProgress(message = "Generating summary", style = "notification", value = 0.1, {
                    for(i in 1:1) {
                      
                      
                      lib1sumdata <- libSummary(dat1)
                      lib2sumdata <- libSummary(dat2)
                      comblibdata2 <- libSummary(newdat[[1]])
                      rtcorrelation <- autortcor()
                      ricorrelation <- autoricor()
                      
                      incProgress(0.5)
                      Sys.sleep(1)
                    }
                  })
                  
                  # if(is.numeric(lib1sumdata[["proteins"]])) {
                  paste0("Seed assay library contains \n", " Proteins = ", formatC(lib1sumdata[["proteins"]], format = "d", big.mark = ","),
                          "\n", " Peptides = ", formatC(lib1sumdata[["peptides"]], format = "d", big.mark = ","), "\n",
                         " Transitions = ",formatC(lib1sumdata[["transitions"]], format = "d", big.mark = ",") ,  "\n", "\n",
                
                         "External assay library contains \n",  " Proteins = ", formatC(lib2sumdata[["proteins"]], format = "d", big.mark = ","),
                         "\n", " Peptides = ", formatC(lib2sumdata[["peptides"]], format = "d", big.mark = ","),  "\n",
                         " Transitions = ", formatC(lib2sumdata[["transitions"]], format = "d", big.mark = ",") , "\n", "\n",
                         
                         "Retention time correlation between seed and external library = ", format(rtcorrelation[[4]], digits = 2) , "\n","\n",
                         
                         "Combined assay library contains \n", " Proteins = ", formatC(comblibdata2[["proteins"]], format = "d", big.mark = ","),
                          "\n", " Peptides = ", formatC(comblibdata2[["peptides"]], format = "d", big.mark = ","),  "\n",
                         " Transitions = " , formatC(comblibdata2[["transitions"]], format = "d", big.mark = ","))
                  # }
                        # else paste(lib1sumdata[["proteins"]])
                }) ########## have to set all lib summary codes
                             
                             })
                
                observeEvent(input$autowizrun, 
                             {
                               output$dvricor <- renderPlot({
                  
                  ricor <- autoricor()
                  
                  # ionCorGS<-NULL; rm(ionCorGS);
                  # data(ionCorGS)
                  # 
                  if(!is.null(ricor)){ 
                    print(ricor[[3]])}
                  else return()
                  
                })
                })
                
#////////////////////////////////////////////////////////////////////////////////////////////////
  
  #### Data table for combined library in auto wizard
                
                observeEvent(input$autowizrun, 
                             {
                               output$autocombineLib <- DT::renderDataTable({
                    newdat = autocombLib()
                   if(!is.null(newdat)) {
                     DT::datatable(data = newdat[[1]], 
                                   options = list(pageLength = 50,
                                                  scrollX = T, 
                                                  scrollY = "500px", 
                                                  scrollCollapse = T, 
                                                  autoWidth = T,
                                                  scroller.loadingIndicator = T),
                                   caption = paste(input$autoseedlib$name,input$autoextlib$name, sep = "_"))
                   } else return()
                   
                    
                 })
                             })
               
  #//////////////////////////////////////////////////////////////////////////             
                 
      ####  Output / autowiz Combined Libraries Download 
                 output$wizdownloadoutputLib <- downloadHandler(filename = function () {
                   paste("Combined", input$autoseedlib$name,input$autoextlib$name , Sys.Date(), ".csv", sep = "_")
                   },
                                                             content = function (file) {
                                                               
                                                                data = autocombLib()
                                                                data = data[[1]]
                                                                
                                                                # write.csv(x = peakviewFormat(data), file = paste0(filename, "peakview" , ".csv"), sep = ",")
                                                               
                                                               # filename <-  paste("Combined", input$autoseedlib$name,input$autoextlib$name , sep = "_")
                                                               
                                                               if(input$autooutputlibformat == "PeakView"){
                                                                 write.csv(x = peakviewFormat(data), file = file)
                                                               } else if(input$autooutputlibformat == "OpenSwath")
                                                               {
                                                                 write.csv(x = OswathFormat(data), file = file)
                                                               } else if(input$autooutputlibformat == "Skyline")
                                                               {
                                                                 write.csv(x = skylineFormat(data), file = file)
                                                               } else # spectronaut
                                                               {
                                                                 write.csv(x = spectronautFormat(data), file = file)
                                                               }
                                                               
                                                               }
                   )
         
  #/////////////////////////////////////////////////////////////////////////////////
                 
                 
#### auto wizard library combination
  
  autocombLib <- eventReactive(input$autowizrun, {

    dat1 <- autolib1read()
    dat2 <- autolib2read()

    if(!is.null(dat1) && !is.null(dat2)) {
      withProgress(message = "Combining libraries", style = "notification", value = 0.1, {
        for(i in 1:1) {

          newdat <- buildSpectraLibPair(baseLib = dat1, extLib = dat2, method = "time",
                                        plot = TRUE, parseAcc = FALSE, consolidateAccession = TRUE)

          incProgress(0.5)
          Sys.sleep(1)
        }
      })

       return(newdat)
      }
    else if(is.null(dat1) && is.null(dat2))
      showNotification("No libraries available", duration = 5, closeButton = TRUE, type = "error")
    else if(is.null(dat1) || is.null(dat2)) showNotification("Library to be combine is missing", duration = 5, closeButton = TRUE, type = "error")
    
    })
  
  
  autortcor <- eventReactive(input$autowizrun, {
    dat1 <- autolib1read()
    dat2 <- autolib2read()
    
    dat <- computeRTTime(dat1, dat2, label1 = input$autoseedlib$name, label2 = input$autoextlib$name)
    
    return(dat)
  })
  autoricor <- eventReactive(input$autowizrun, {
    dat1 <- autolib1read()
    dat2 <- autolib2read()
    
    dat <- computeIntensityCor(dat1, dat2)
    
    return(dat)
  })

#////////////////////////////////////////////////////////////////////////////////////////    
  
  seedlibdata <- eventReactive(input$apply, {
    lib1clean()
  })
  extlibdata <- eventReactive(input$apply,{
    lib2clean()
  })

#//////////////////////////////////////////////////////////////////////////////////
  
#### Clean button for cleaning update
  
  # seedcleanlibdata <- eventReactive(input$cleanupdate, {
  #   lib1clean()
  # })
  # extcleanlibdata <- eventReactive(input$cleanupdate,{
  #   lib2clean()
  # })
  
  #//////////////////////////////////////////////////////////////////////////////////
  
  
#### auto Seed library reading reactive
  
  autolib1read <- eventReactive(input$autowizrun,
                                {
  inFile <- input$autoseedlib
  if(is.null(inFile))
    return(NULL)
    else{
  withProgress(message = "Reading seed library", style = "notification", value = 0.1, {
    for(i in 1:1) {
      
      lib <- readLibFile(inFile$datapath, input$autoformlib1, "spectrum", clean = TRUE)
      
      incProgress(0.75)
      Sys.sleep(1)
    }
  }
  )}
  return(lib)
  })

  ###////////////////////////////////////////////////////////////////
  
  #### auto Ext library reading reactive
  
  autolib2read <- eventReactive(input$autowizrun,
                                {
    inFile <- input$autoextlib
    if(is.null(inFile))
      return(NULL)
    else{
      withProgress(message = "Reading External library", style = "notification", value = 0.1, {
        for(i in 1:1) {
          
          lib <- readLibFile(inFile$datapath, input$autoformlib2, "spectrum", clean = TRUE)
          
          incProgress(0.75)
          Sys.sleep(1)
        }
      }
      )}
    return(lib)
  })
  
  ###////////////////////////////////////////////////////////////////
  
#### Seed Library reading reactive
  
  lib1read <- reactive({
    
    
    inFile <- input$seedlib
  
  if(is.null(inFile)) {
    
    return(NULL)
  } else
  
  readLibFile(inFile$datapath, input$libformats11, "spectrum", clean = FALSE)
  })

#/////////////////////////////////////////////////////////////////////////////////  
  
#### Seed Library cleaning reactive

  lib1clean <- reactive(
  {
    lib1 <- lib1read()
    cleanLib(lib1, TRUE, intensity.cutoff = input$intensity,
             conf.cutoff = input$confidence, nomod = input$modified,
             nomc = input$misscleavage, enz = input$enzyme)
  })
  

#/////////////////////////////////////////////////////////////////////////////////
  
#### External Library Reading reactive
  
  lib2read <- reactive({
    inFile <- input$extlib
  
  if(is.null(inFile)){
    return(NULL)
    } else
  
  readLibFile(inFile$datapath, input$libformats12, "spectrum", clean = FALSE)
  })
  
#/////////////////////////////////////////////////////////////////////////////////  
  
#### External Library Claening reactive

  lib2clean <- reactive(
  {
    lib2 <- lib2read()
    cleanLib(lib2, TRUE, intensity.cutoff = input$intensity,
             conf.cutoff = input$confidence, nomod = input$modified,
             nomc = input$misscleavage, enz = input$enzyme)
  })

#/////////////////////////////////////////////////////////////////////////////////

#### Enable apply button    
  output$readapply <- renderUI({
    if (is.null(lib1clean()) && is.null(lib2clean())) return()
    shinyjs::enable("apply")
  })  
  
#/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable downloadseed button    
  output$downloadseed <- renderUI({
    if (is.null(lib1clean()) && is.null(lib2clean())) return()
    else if(!is.null(lib1clean()))
    shinyjs::enable("downloadseedinputLibs")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable downloadseed button    
  output$downloadext <- renderUI({
    if (is.null(lib1clean()) && is.null(lib2clean())) return()
    else if(!is.null(lib2clean()))
      shinyjs::enable("downloadextinputLibs")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
#### Enable auto run button    
  # output$autoreadapply <- renderUI({
  #   if (is.null(autolib1read()) && is.null(autolib2read())) return()
  #   else shinyjs::enable("autowizrun")
  # })  
  
#/////////////////////////////////////////////////////////////////////////////////
  
#### Enable auto download button    
  # output$autoreadapply <- renderUI({
  #   if (is.null(autolib1read()) && is.null(autolib2read())) return()
  #   else shinyjs::enable("downloadseedinputLibs")
  # })  
  
#/////////////////////////////////////////////////////////////////////////////////
  
#### Enable include length option    
  output$includelength <- renderUI({
    if (input$combmethod == "hydro" || input$combmethod == "hydrosequence")
    shinyjs::enable("inclength")
  })  
  
#/////////////////////////////////////////////////////////////////////////////////
  
#### Seed Library Data Table
  
  output$seedlibcontents <- DT::renderDataTable({
    
    if (!is.null(lib1clean())){
    withProgress(message = "Generating seed library table", style = "notification", value = 0.1, {
      for(i in 1:1) {
        data = seedlibdata()
        incProgress(0.1, detail = "reading library")
        Sys.sleep(0.25)
      }
    })
    
    DT::datatable(data = data,
                  options = list(pageLength = 10, 
                                 scrollX = T, 
                                 scrollY = "500px", 
                                 scrollCollapse = T, 
                                 autoWidth = T,
                                 scroller.loadingIndicator = T),
                  caption = input$seedlib$name
                  ) 
    }
    else return()
})

  
#/////////////////////////////////////////////////////////////////////////////////  
  
#### External Library Data Table 
  
  output$extlibcontents <- DT::renderDataTable({
   
    if (!is.null(lib2clean())){ 
    withProgress(message = "Generating external library table", style = "notification", value = 0.1, {
      for(i in 1:1) {
        data = extlibdata()
        incProgress(0.1, detail = "reading library")
        Sys.sleep(0.25)
      }
    })
    
    DT::datatable(data = data,
                  options = list(pageLength = 10,
                                 scrollX = T, 
                                 scrollY = "500px", 
                                 scrollCollapse = T, 
                                 autoWidth = T,
                                 scroller.loadingIndicator = T),
                  caption = input$extlib$name
                  )
    } else return()
  })

#/////////////////////////////////////////////////////////////////////////////////  
  
#### Seed Lib Plot

   output$lib1plot <- renderPlot({
     
     if (!is.null(lib1clean())){
     withProgress(message = "Generating seed library plot", style = "notification", value = 0.1, {
       for(i in 1:1) {
         lib1data <- seedlibdata()
         x<-aggregate(lib1data$uniprot_id, by=list(Protein=lib1data$uniprot_id, Peptide=lib1data$stripped_sequence), 
                      FUN=function(x){length(x)})
         y<-aggregate(x$Peptide, by=list(Protein=x$Protein), FUN=function(x){length(x)})
         incProgress(0.1, detail = "plotting")
         Sys.sleep(0.25)
       }
     })
       
       p1 <- ggplot(data = x, aes(x = x)) +
         geom_histogram(aes(fill = ..count..), binwidth = 5) +
         labs(x ="Number of ions per peptide") +
         scale_y_continuous(name = "Frequency") +
         scale_fill_gradient("Frequency", low = "blue", high = "red")+
         theme_classic() +
         theme(plot.title = element_text(size = 14, face = "bold"),
               axis.title = element_text(size = 11, face = "bold"),
               axis.text.x = element_text(face = "bold", size = 12),
               axis.text.y = element_text(face= "bold", size = 12)) 
       
       p2 <- ggplot(data = y, aes(x = x)) +
         geom_histogram(aes(fill = ..count..), binwidth = 5) +
         labs(x ="Number of peptides per protein") +
         scale_y_continuous(name = "Frequency") +
         scale_fill_gradient("Frequency", low = "blue", high = "red")+
         theme_classic() +
         theme(plot.title = element_text(size = 14, face = "bold"),
               axis.title = element_text(size = 11, face = "bold"),
               axis.text.x = element_text(face = "bold", size = 12),
               axis.text.y = element_text(face= "bold", size = 12)) 
       
       p <- ggarrange(p1, p2, ncol = 2)
       print(p)
       
       filepath <- getwd()
       filepath2 <- paste0(filepath,"/graphs")
       ggsave("Seed_library_plot.png", width = 8, height = 5, path = filepath2)
       
       
       
       
    # Relative_Intensity <- lib1data$relative_intensity
    # Retention_Time <- lib1data$RT_detected
    # mass_to_charge <- lib1data$frg_z
    # Fragment_Type <- factor(lib1data$frg_type)
    # libplot <- ggplot(data = lib1data, aes(x = Retention_Time, 
    #                             y = Relative_Intensity, 
    #                             #z = mass_to_charge 
    #                             col = Fragment_Type)) + 
    #   geom_line()
    # libplot + theme(panel.background = element_rect(fill = 'white', colour = 'red'))
     }else return()
  })

  
#/////////////////////////////////////////////////////////////////////////////////  
  
#### External Library Plot  
  
 output$lib2plot <- renderPlot({
   
   if (!is.null(lib2clean)){
   withProgress(message = "Generating external library plot", style = "notification", value = 0.1, {
     for(i in 1:1) {
       lib2data <- extlibdata()
       x<-aggregate(lib2data$uniprot_id, by=list(Protein=lib2data$uniprot_id, Peptide=lib2data$stripped_sequence), 
                    FUN=function(x){length(x)})
       y<-aggregate(x$Peptide, by=list(Protein=x$Protein), FUN=function(x){length(x)})
       incProgress(0.1, detail = "plotting")
       Sys.sleep(0.25)
     }
   })
     
     p1 <- ggplot(data = x, aes(x = x)) +
       geom_histogram(aes(fill = ..count..), binwidth = 5) +
       labs(x ="Number of ions per peptide") +
       scale_y_continuous(name = "Frequency") +
       scale_fill_gradient("", low = "blue", high = "red")+
       theme_classic() +
       theme(plot.title = element_text(size = 14, face = "bold"),
             axis.title = element_text(size = 11, face = "bold"),
             axis.text.x = element_text(face = "bold", size = 12),
             axis.text.y = element_text(face= "bold", size = 12)) 
     
     p2 <- ggplot(data = y, aes(x = x)) +
       geom_histogram(aes(fill = ..count..), binwidth = 5) +
       labs(x ="Number of peptides per protein") +
       scale_y_continuous(name = "Frequency") +
       scale_fill_gradient("", low = "blue", high = "red")+
       theme_classic() +
       theme(plot.title = element_text(size = 14, face = "bold"),
             axis.title = element_text(size = 11, face = "bold"),
             axis.text.x = element_text(face = "bold", size = 12),
             axis.text.y = element_text(face= "bold", size = 12)) 
     p <- ggarrange(p1, p2, ncol = 2)
     print(p)
     
     filepath <- getwd()
     filepath2 <- paste0(filepath,"/graphs")
     ggsave("External_library_plot.png", width = 8, height = 5, path = filepath2)
     
    # Relative_Intensity <- lib2data$relative_intensity
    # Retention_Time <- lib2data$RT_detected
    # mass_to_charge <- lib2data$frg_z
    # Fragment_Type <- factor(lib2data$frg_type)
    # libplot <- ggplot(data = lib2data, aes(x = Retention_Time, 
    #                                    y = Relative_Intensity, 
    #                                    #z = mass_to_charge
    #                                    col = Fragment_Type)) + 
    #          geom_line()
    # 
    # libplot + theme(panel.background = element_rect(fill = 'white', colour = 'red'))
   } else return()
  })
 
#/////////////////////////////////////////////////////////////////////////////////
 
 
# #### Hover and click function (seed lib)
#  
#  output$seed_info <- renderText({
#    xy_str <- function(e) {
#      if(is.null(e)) return("NULL\n")
#      paste0("Retention time =", round(e$x, 1), " Relative intensity =", round(e$y, 1), "\n")
#    }
#    xy_range_str <- function(e) {
#      if(is.null(e)) return("NULL\n")
#      paste0("Min_retention time =", round(e$xmin, 1), " Max_retention time =", round(e$xmax, 1), 
#             " Min_relative intensity =", round(e$ymin, 1), " Max_relative intensity =", round(e$ymax, 1))
#    }
#    
#    paste0(
#      "click: ", xy_str(input$seed_click),
#      "dblclick: ", xy_str(input$seed_dblclick),
#      "hover: ", xy_str(input$seed_hover),
#      "brush: ", xy_range_str(input$seed_brush)
#    )
#  })
#  
# #///////////////////////////////////////////////////////////////////////////////// 
#  
# #### Hover and click function  (seed lib)
#  
#  output$ext_info <- renderText({
#    xy_str <- function(e) {
#      if(is.null(e)) return("NULL\n")
#      paste0("Retention time =", round(e$x, 1), " Relative intensity =", round(e$y, 1), "\n")
#    }
#    xy_range_str <- function(e) {
#      if(is.null(e)) return("NULL\n")
#      paste0("Min_retention time =", round(e$xmin, 1), " Max_retention time =", round(e$xmax, 1), 
#             " Min_relative intensity =", round(e$ymin, 1), " Max_relative intensity =", round(e$ymax, 1))
#    }
#    
#    paste0(
#      "click: ", xy_str(input$ext_click),
#      "dblclick: ", xy_str(input$ext_dblclick),
#      "hover: ", xy_str(input$ext_hover),
#      "brush: ", xy_range_str(input$ext_brush)
#    )
#  })
#  
# #///////////////////////////////////////////////////////////////////////////////// 
 

#### All other Hovers and clicks function  (seed lib)
 
 addTooltip(session = session, id = "rtcorplot", 
           title = "Graphs can be downloaded in 'Data Visualization' menu",
           placement = "top", trigger = "hover")
 
 addTooltip(session = session, id = "rtresdplot", 
            title = "Graphs can be downloaded in 'Data Visualization' menu",
            placement = "top", trigger = "hover")
 
 addTooltip(session = session, id = "intensitycorplot", 
            title = "Graphs can be downloaded in 'Data Visualization' menu",
            placement = "top", trigger = "hover")
 
 addTooltip(session = session, id = "inclength", 
            title = "Include length of peptides in predicting retention times",
            placement = "right", trigger = "hover")
 
 # addTooltip(session = session, id = "parseacc", 
 #            title = "Consolidate the protein accessions between two libraries",
 #            placement = "top", trigger = "hover")
 # 
 addTooltip(session = session, id = "consacc", 
            title = "Consolidate the protein accessions between two libraries",
            placement = "right", trigger = "hover")
 
 addTooltip(session = session, id = "gplots", 
            title = "Generate plots during processing. Can be visualized in 'Data Visualization' menu",
            placement = "right", trigger = "hover")
 
 addTooltip(session = session, id = "mergelib", 
            title = "Merge two libraries after predicting retention times",
            placement = "right", trigger = "hover")
 
 addTooltip(session = session, id = "recalrt", 
            title = "Recalibrate retention times of external library before combining in case of low correlation",
            placement = "right", trigger = "hover")
 
#///////////////////////////////////////////////////////////////////////////////// 
 
  
#### PeakView and OpenSwath Library Conversion and Writing/Saving reactives  
  
  # peaklib1 <- reactive({peakviewFormat(lib1clean())})
  # 
  # openswathlib1 <- reactive({OswathFormat(lib1clean())})
  # 
  # skylinelib1 <- reactive({skylineFormat(lib1clean())})
  # 
  # spectlib1 <- reactive({spectronautFormat(lib1clean())})
  # 
  # peaklib2 <- reactive({peakviewFormat(lib2clean())})
  # 
  # openswathlib2 <- reactive({OswathFormat(lib2clean())})
  # 
  # skylinelib2 <- reactive({skylineFormat(lib2clean())})
  # 
  # spectlib2 <- reactive({spectronautFormat(lib2clean())})
  
  
#///////////////////////////////////////////////////////////////////////////////// 
  
####  Input Libraries Download 
 output$downloadseedinputLibs <- downloadHandler(filename = function()
 {paste("Library_seed",input$seedlib$name , Sys.Date(), ".csv", sep = "_")},
 
 content= function(file)
 {
   data = lib1clean()
   
   if(input$inputseedlibformat == "PeakView"){
     write.csv(x = peakviewFormat(data), file = file)
     
   } else if(input$inputseedlibformat == "OpenSwath")
   {
     write.csv(x = OswathFormat(data), file = file)
     
   } else if(input$inputseedlibformat == "Skyline")
   {
     write.csv(x = skylineFormat(data), file = file)
     
   } else # spectronaut
   {
     write.csv(x = spectronautFormat(data), file = file)
   }
   
 }
 )
  
  output$downloadextinputLibs <- downloadHandler(filename = function()
    {paste("Library_ext",input$extlib$name , Sys.Date(), ".csv", sep = "_")},
    
    content= function(file)
    {
      data = lib2clean()
      
      if(input$inputextlibformat == "PeakView"){
        write.csv(x = peakviewFormat(data), file = file)
        
      } else if(input$inputextlibformat == "OpenSwath")
      {
        write.csv(x = OswathFormat(data), file = file)
        
      } else if(input$inputextlibformat == "Skyline")
      {
        write.csv(x = skylineFormat(data), file = file)
        
      } else # spectronaut
      {
        write.csv(x = spectronautFormat(data), file = file)
      }
      
    }
    )
  
#/////////////////////////////////////////////////////////////////////////////////
  
  
####  Output / Combined Libraries Download 
  output$downloadoutputLib <- downloadHandler(filename = function () 
    {paste("Combined", input$seedlib$name,input$extlib$name , Sys.Date(), ".csv", sep = "_")
    },
                                              content = function (file) {
                                                
                                                data = combdatalib()
                                                
                                                if(input$outputlibformat == "PeakView"){
                                                  write.csv(x = peakviewFormat(data), file = file)
                                                } else if(input$outputlibformat == "OpenSwath")
                                                {
                                                  write.csv(x = OswathFormat(data), file = file)
                                                } else if(input$outputlibformat == "Skyline")
                                                {
                                                  write.csv(x = skylineFormat(data), file = file)
                                                } else # spectronaut
                                                {
                                                  write.csv(x = spectronautFormat(data), file = file)
                                                }
                                                
                                              }
    )
  
#/////////////////////////////////////////////////////////////////////////////////
  
  
#### Available Libraries
  
  output$libseedfilesrt <- renderUI({
    if (is.null(lib1clean())) return()
    else
    paste(input$seedlib$name)
  })
  
  output$libextfilesrt <- renderUI({
    if (is.null(lib2clean())) return()
    else
      paste(input$extlib$name)
  })

#//////////////////////////////////////////////////////////////////////////////////
  
#### Available Libraries
  
  output$libseedfilesint <- renderUI({
    if (is.null(lib1clean())) return()
    else
    paste(input$seedlib$name)
  })
  
  output$libextfilesint <- renderUI({
    if (is.null(lib2clean())) return()
    else
      paste(input$extlib$name)
  })
  
#//////////////////////////////////////////////////////////////////////////////////
  
#### Available Libraries for stats
  
  output$statlibnames <- renderUI({
    
    if(input$inputfiles == "current"){
    if (is.null(input$seedlib$name) && is.null(input$extlib$name)) return()
      else{
        extendedlib <- paste(input$seedlib$name,input$extlib$name, sep = "_")
        paste(input$seedfile$name, "," ,extendedlib )}
    } 
    
    if(input$inputfiles == "new"){
      if (is.null(input$seedfile$name) && is.null(input$extendedfile$name)) return()
      else
        paste(input$seedfile$name, "," , input$extendedfile$name)
    } 
  })
  
#//////////////////////////////////////////////////////////////////////////////////
  
#### Seed Library Summary
  
  output$seedlibsummary <- renderText({
    seeddata <- lib1clean()
    summarydata <- libSummary(seeddata)
    
   paste0("Seed assay library contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
         "\n", " Peptides = ", formatC(summarydata[["peptides"]], format = "d", big.mark = ","),  "\n",
         " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
   
  })
  
#//////////////////////////////////////////////////////////////////////////////////

#### External Library Summary
  
  output$extlibsummary <- renderText({
    extdata <- lib2clean()
    summarydata <- libSummary(extdata)
    paste0("External assay library contains \n", " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","), 
           "\n", " Peptides = ", formatC(summarydata[["peptides"]], format = "d", big.mark = ","), "\n",
           " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
    
  })
  
#//////////////////////////////////////////////////////////////////////////////////  
  
#### Compute RT Time reactive function
  
  ### by time tables/ graphs
  
  rtcordt <- eventReactive(input$calRT, {
    dat1 <- lib1clean()
    dat2 <- lib2clean()

    if(!is.null(dat1) && !is.null(dat2)){                                   
      withProgress(message = "Retention time correlation", style = "notification", value = 0.1, {
        for(i in 1:1) {
          
          if(input$rttime == "time") {
            rtcordata <- computeRTTime(dat1, dat2, label1 = input$seedlib$name, label2 = input$extlib$name)
          } 
          incProgress(0.1, detail = "computing")
          # Sys.sleep(0.25)
        }
      })
      
      if(is.null(rtcordata)) return ()
      else{
        return(rtcordata)
      }
      
    } else if(is.null(dat1) && is.null(dat2))
      showNotification("No libraries available", duration = 5, closeButton = TRUE, type = "error")
    else if(is.null(dat1) || is.null(dat2)) showNotification("Library to be compared is missing", duration = 5, closeButton = TRUE, type = "error") 
  })

#//////////////////////////////////////////////////////////////////////////////////
  
  #### Compute RT Time Hydro reactive function
  
  ###  by hydro tables/ graphs
  
  hydrocordt <- eventReactive(input$calRT,{
    dat1 <- lib1clean()
    dat2 <- lib2clean()
    
    if(!is.null(dat1) && !is.null(dat2) && !is.null(hydrofile())){                                   
      withProgress(message = "Retention time correlation", style = "notification", value = 0.1, {
        for(i in 1:1) {
          
          if(input$rttime == "hydro") {
            rtcordata <- computeRTHydro(dat1, dat2, datHydroIndex = hydrofile(), label1 = input$seedlib$name, label2 = input$extlib$name)
          }
          incProgress(0.1, detail = "computing")
          # Sys.sleep(0.25)
        }
      })
      
      if(is.null(rtcordata)) return ()
      else{
        return(rtcordata)
      }
      
    } else if(is.null(dat1) && is.null(dat2))
      showNotification("No libraries available", duration = 5, closeButton = TRUE, type = "error")
    else if(is.null(dat1) || is.null(dat2)) showNotification("Library to be compared is missing", duration = 5, closeButton = TRUE, type = "error") 
  })
  
  #//////////////////////////////////////////////////////////////////////////////////
  
  
  
#### hydro file reactive function
  
  hydrofile <- reactive({
  
    dat1 <- lib1clean()
    dat2 <- lib2clean()
    
    inFile <- input$hydrofile
    
    if(is.null(inFile))
      return(NULL)
    
    hydrodat <- read.delim(inFile$datapath, stringsAsFactors = F)
    
    hydrodat
    
  })
  
#//////////////////////////////////////////////////////////////////////////////////
  
  
#### Retention Time Correlation graph
  
  observeEvent(input$calRT,{
    output$rtcorplot <- renderPlot({
    
    
    rr <- rtcordt()
    if (!is.null(rr)){
    
      print(rr[[1]]) 
  } else return()
  })
})
  

   output$dvrtcor <- renderPlot({
     dat1 <- lib1clean()
     dat2 <- lib2clean()
     
     if(!is.null(dat1) && !is.null(dat2))
     rr <- rtcordt()
     else rr <- autortcor()
     
     if (!is.null(rr)){
       
       print(rr[[1]]) 
     } else return()
   })
 
   
  

#/////////////////////////////////////////////////////////////////////////////////
  
  #### Retention Time Hydro Correlation graph
  
  observeEvent(input$calRT,{
    output$hydrocorplot <- renderPlot({
    
    
    rr <- hydrocordt()
    if (!is.null(rr)){
      
        print(rr[[2]])
    } else return()
  })
  })
  
 
    output$dvhydrocor <- renderPlot({

        rr <- hydrocordt()
      if (!is.null(rr)){
        
        print(rr[[2]])
      } else return()
    })
  
  
  
  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
#### Retention Time Residual Plot
  
  observeEvent(input$calRT,{
    output$rtresdplot <- renderPlot({
    
    rr <- rtcordt()
    if (!is.null(rr)){
   
     print(rr[[2]])
 
    } else return()
  })
  })
  
    output$dvrtresd <- renderPlot({
      
      dat1 <- lib1clean()
      dat2 <- lib2clean()
      
      if(!is.null(dat1) && !is.null(dat2))
        rr <- rtcordt()
      else rr <- autortcor()

      if (!is.null(rr)){
        
        print(rr[[2]])
        
      } else return()
    })
  
  
  
#/////////////////////////////////////////////////////////////////////////////////  
  
  #### Retention Time Hydro Residual Plot
  
  observeEvent(input$calRT,{
               output$hydroresdplot <- renderPlot({
    
    rr <- hydrocordt()
    if (!is.null(rr)){
      
        print(rr[[3]])
    } else return()
  })
  })
  
  
    output$dvhydroresd <- renderPlot({
      
        rr <- hydrocordt()
      if (!is.null(rr)){
        
        print(rr[[3]])
      } else return()
    })
  
  
  
  
  #/////////////////////////////////////////////////////////////////////////////////  
  
#### Retention Time Correlation datatable

  observeEvent(input$calRT, {
    output$rtcordata <- DT::renderDataTable({
      
      if(is.null(rtcordt())) 
        return ()
      
      else {
          
        data <- rtcordt() 
          data <- data[[3]]
        caption <- "Retention times of common peptides in Seed and External libraries (Training Set)"
      
        
        if(is.null(data)) {return ()
          }else {
    DT::datatable(data = data,
                  options = list(pageLength = 50,
                                 scrollX = T, 
                                 scrollY = "500px", 
                                 scrollCollapse = T, 
                                 autoWidth = T,
                                 scroller.loadingIndicator = T),
                  caption = caption)
          }
      }
  })
  
})
  
#/////////////////////////////////////////////////////////////////////////////////
  
  
  #### Retention Time Hydro Correlation datatable
  
  observeEvent(input$calRT, {
    output$hydrocordata <- DT::renderDataTable({
      
      if(is.null(hydrocordt())) 
        return ()
      
      else {
       
          data <- hydrocordt()
          data <- data[[1]]
          caption <- "Retention times and hydrophobicity values of common peptides in Seed and External libraries (Training Set)"
        
        if(is.null(data)) {return ()
        }else {
          DT::datatable(data = data,
                        options = list(pageLength = 50,
                                       scrollX = T, 
                                       scrollY = "500px", 
                                       scrollCollapse = T, 
                                       autoWidth = T,
                                       scroller.loadingIndicator = T),
                        caption = caption)
        }
      }
    })
    
  })
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Compute RI reactive function
  
  computeRIbyintensity <- reactive({
    
    dat1 <- lib1clean()
    dat2 <- lib2clean() 
    
    if(!is.null(dat1) && !is.null(dat2)){
      withProgress(message = "Relative intensity correlation", style = "notification", value = 0.1, {
        for(i in 1:1) {
          if(input$rintensity == "spearman") {
            
            ri <- computeIntensityCor(dat1, dat2, method = "spearman")
            
          } else {
            ri <- computeIntensityCor(dat1, dat2, method = "kendall")
          }
          
          incProgress(0.1, detail = "computing")
          Sys.sleep(0.25)
        }
      })
    }
    else if(is.null(dat1) && is.null(dat2))
      showNotification("No libraries available", duration = 5, closeButton = TRUE, type = "error")
    else if(is.null(dat1) || is.null(dat2)) showNotification("Library to be compared is missing", duration = 5, closeButton = TRUE, type = "error")
    
    
    return(ri)
  })
  
  #//////////////////////////////////////////////////////////////////////////////////
  
#### Relative Intensity Correlation Graph
  
  observeEvent(input$calRI,{
    output$intensitycorplot <- renderPlot({
 
    
   
          ricor <- computeRIbyintensity()
          
        ionCorGS<-NULL; rm(ionCorGS);
        data(ionCorGS)
      # 
        if(!is.null(ricor)){ 
       # 
       # dataplot
        
        # allCordata <- data.frame(dataplot)


        # p <- ggplot(data = allCordata, aes(x= "", y = allCordata)) +
        #   geom_boxplot()
        
        print(ricor[[3]])}
        else return()
       
         # boxplot(dataplot, names="in study",
        #         col="darkgreen", main="Relative ion intensity correlation between libraries",outline=FALSE)
        # abline(h=0.8, col="red")
      
      
  })
})
  
  output$dvricor <- renderPlot({

    ricor <- computeRIbyintensity()

    ionCorGS<-NULL; rm(ionCorGS);
    data(ionCorGS)
    #
    if(!is.null(ricor)){
      print(ricor[[3]])}
    else return()

  })
  
#/////////////////////////////////////////////////////////////////////////////////// 
  
#### Relative Intensity Correlation datatable
  
  observeEvent(input$calRI,{
    output$intensitycordata <- DT::renderDataTable({
    dat1 <- lib1clean()
    dat2 <- lib2clean()
    
    if(!is.null(dat1) && !is.null(dat2)){
      ricor <- computeRIbyintensity()
      
    DT::datatable(data = ricor[[2]],
                  options = list(pageLength = 50,
                                 scrollX = T, 
                                 scrollY = "500px", 
                                 scrollCollapse = T, 
                                 autoWidth = T,
                                 scroller.loadingIndicator = T),
                  caption = "Relative Intensities of common peptides in Seed and External libraries")
    } else return()
  })
})
  
#////////////////////////////////////////////////////////////////////////////////// 

#### Building Combined Library reactive
  
  combdatalib <- eventReactive(input$combrun, { 
    dat1 <- lib1clean()
    dat2 <- lib2clean()
    
    if(!is.null(dat1) && !is.null(dat2)){
    
    withProgress(message = "Combining libraries", style = "notification", value = 0.1, {
      for(i in 1:1) {
        if(input$combmethod == "time"){
          
          newdat <- buildSpectraLibPair(baseLib = dat1, extLib = dat2, method = input$combmethod, merge = input$mergelib,
                                        includeLength = input$inclength, recalibrateRT = input$recalrt, cutoff.size = input$cutofftsize, cutoff.r2 = input$cutoffr2, 
                                        plot = input$gplots, consolidateAccession = input$consacc)}
        
        else { 
          inFile <- input$combhydrofile
          
          if(is.null(inFile))
            return(NULL)
          
          datHydroIndex <- read.delim(inFile$datapath, stringsAsFactors = F)
          newdat <- buildSpectraLibPair(baseLib = dat1, extLib = dat2, method = input$combmethod, 
                                        includeLength = input$inclength, hydroIndex = datHydroIndex, recalibrateRT = input$recalrt,
                                        cutoff.size = input$cutofftsize, cutoff.r2 = input$cutoffr2,
                                        plot = input$gplots, parseAcc = input$parseacc, consolidateAccession = input$consacc)}
        
        incProgress(0.5)
        Sys.sleep(1)
      }
    })
      
      
    return(newdat)
    } else if(is.null(dat1) && is.null(dat2))
      showNotification("No libraries available to combine", duration = 5, closeButton = TRUE, type = "error")
    else if(is.null(dat1) || is.null(dat2)) showNotification("Library to be combined is missing", duration = 5, closeButton = TRUE, type = "error")
  })
    
    output$combineLib <- DT::renderDataTable({
      if(!is.null(lib1clean()) && !is.null(lib2clean())){
        data = combdatalib()
      DT::datatable(data = data[[1]], 
                    options = list(pageLength = 50,
                                   scrollX = T, 
                                   scrollY = "500px", 
                                   scrollCollapse = T, 
                                   autoWidth = T,
                                   scroller.loadingIndicator = T),
                    caption = paste(input$seedlib$name,input$extlib$name, sep = "_"))
      } else return()
      
    })
    
    observeEvent(input$autowizrun, {
      
      output$dvhisto <- renderPlot({
      data = autocombLib()
      data2 <- data[[2]]
      if(is.null(data2))
        return()
      else print(data2$phist)
      })
    })
    
    observeEvent(input$combrun, {
      
      output$dvhisto <- renderPlot({
        data = combdatalib()
        data2 <- data[[2]]
        if(is.null(data2))
          return()
        else print(data2$phist)
      })
    })
    
 
    
    observeEvent(input$autowizrun, {
      
      output$dvdensity <- renderPlot({
        data = autocombLib()
        data2 <- data[[2]]
        if(is.null(data2))
          return()
        else print(data2$pdens)
      })
    })
    
    observeEvent(input$combrun, {
      
      output$dvdensity <- renderPlot({
        data = combdatalib()
        data2 <- data[[2]]
        if(is.null(data2))
          return()
        else print(data2$pdens)
      })
    })
    
    observeEvent(input$autowizrun, {
      
      output$dvpepprotnum <- renderPlot({
        data = autocombLib()
        data2 <- data[[2]]
        if(is.null(data2))
          return()
        else print(data2$ppnum)
      })
    })
    
    observeEvent(input$combrun, {
      
      output$dvpepprotnum <- renderPlot({
        data = combdatalib()
        data2 <- data[[2]]
        if(is.null(data2))
          return()
        else print(data2$ppnum)
      })
    })
    
    observeEvent(input$autowizrun, {

      output$dvpepvenn <- renderImage({
        data = autocombLib()
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/peptideVennDigram.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })

    observeEvent(input$combrun, {

      output$dvpepvenn <- renderImage({
        data = combdatalib()
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/peptideVennDigram.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })

    observeEvent(input$autowizrun, {

      output$dvprotvenn <- renderImage({
        data = autocombLib()
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/proteinVennDigram.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })

    observeEvent(input$combrun, {

      output$dvprotvenn <- renderImage({
        data = combdatalib()
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/proteinVennDigram.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })

    observeEvent(input$autowizrun, {

      output$dvlibinfo <- renderImage({
        # data = autocombLib()
        # dat1 = autolib1read()
        # reliabilityCheckLibrary(dat1, data)
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })

    observeEvent(input$combrun, {

      output$dvlibinfo <- renderImage({
        # data = combdatalib()
        # dat1 = lib1clean()
        # reliabilityCheckLibrary(dat1, data)
        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
        list(src = filepath2, width = "100%", height = "100%")
      })
    })
 
  
# ///////////////////////////////////////////////////////////////////////////////////
    
    
#### Graphs Download
    
    # output$downloadextinputLibs <- downloadHandler(filename = function()
    # {paste("Library_ext",input$extlib$name , Sys.Date(), ".csv", sep = "_")},
    # 
    # content= function(file)
    # {
    #   data = lib2clean()
    #   
    #   if(input$inputextlibformat == "PeakView"){
    #     write.csv(x = peakviewFormat(data), file = file)
    #     
    #   } else if(input$inputextlibformat == "OpenSwath")
    #   {
    #     write.csv(x = OswathFormat(data), file = file)
    #     
    #   } else if(input$inputextlibformat == "Skyline")
    #   {
    #     write.csv(x = skylineFormat(data), file = file)
    #     
    #   } else # spectronaut
    #   {
    #     write.csv(x = spectronautFormat(data), file = file)
    #   }
    #   
    # }
    # )
    
    # output$downloadlibplots <- downloadHandler(filename = function()
    # {paste("Plot", Sys.Date(), ".tiff", sep = "_")},
    # 
    # content= function(file)
    # {
    #   observeEvent(input$combrun, {
    #     data = combdatalib()
    #     data2 <- data[[2]]
    #   
    #  if(input$libanalysisplots == "Histograms")
    #   {
    #    tiff(file, width = 480, height = 480, units = "px", bg = "white")
    #    plot(data2$phist)
    #    dev.off()
    #     
    #   } else if(input$libanalysisplots == "Density plots")
    #   {
    #     tiff(file, data2$pdens, width = '100%', height = '100%', bg = "white")
    #     
    #   } else if(input$libanalysisplots == "Peptide protein numbers")
    #   {
    #     tiff(file, data2$ppnum, width = '100%', height = '100%', bg = "white")
    #     
    #   } else if(input$libanalysisplots == "Peptide venn diagram")
    #   {
    #     filepath <- getwd()
    #     filepath2 <- paste0(filepath,"/graphs/peptideVennDigram.png")
    #     tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #     
    #   }else if(input$libanalysisplots == "Protein venn diagram")
    #   {
    #     filepath <- getwd()
    #     filepath2 <- paste0(filepath,"/graphs/proteinVennDigram.png")
    #     tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #     
    #   } else # Libraries coverage
    #   {
    #     filepath <- getwd()
    #     filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
    #     tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #   }
    #   
    #   })
    #   
    #   observeEvent(input$autowizrun, {
    #     data = autocombLib()
    #     data2 <- data[[2]]
    #     
    #     if(input$libanalysisplots == "Histograms")
    #     {
    #       tiff(file, width = 480, height = 480, units = "px", bg = "white")
    #       plot(data2$phist)
    #       dev.off()
    #       
    #     } else if(input$libanalysisplots == "Density plots")
    #     {
    #       tiff(file, data2$pdens, width = '100%', height = '100%', bg = "white")
    #       
    #     } else if(input$libanalysisplots == "Peptide protein numbers")
    #     {
    #       tiff(file, data2$ppnum, width = '100%', height = '100%', bg = "white")
    #       
    #     } else if(input$libanalysisplots == "Peptide venn diagram")
    #     {
    #       filepath <- getwd()
    #       filepath2 <- paste0(filepath,"/graphs/peptideVennDigram.png")
    #       tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #       
    #     }else if(input$libanalysisplots == "Protein venn diagram")
    #     {
    #       filepath <- getwd()
    #       filepath2 <- paste0(filepath,"/graphs/proteinVennDigram.png")
    #       tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #       
    #     } else # Libraries coverage
    #     {
    #       filepath <- getwd()
    #       filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
    #       tiff(file, filepath2, width = '100%', height = '100%', bg = "white")
    #     }
    #     
    #   })
    #   }, contentType = 'image/tiff'
    # )
  

## //////////////////////////////////////////////////////////////////////////////////        
    
  
# #### Library Reliability Check
#   
  # librarychecking <- eventReactive(
  #   input$librelcheck, {
  #     if(input$inputfiles == "current") {
  # 
  #        dat1 <- lib1clean()
  #        if(!is.null(dat1))
  #     dat2 <- combdatalib()
  #        else {
  #          dat1<- autolib1read()
  #          dat2 <- autocombineLib()
  #        }
  # 
  #     } else
  #     {
  #       inFile1 <- input$seedfile
  # 
  #       if(is.null(inFile1))
  #         return(NULL)
  # 
  #       dat1 <- inFile1$datapath
  # 
  #       inFile2 <- input$extendedfile
  # 
  #         if(is.null(inFile2))
  #       return(NULL)
  # 
  #       dat2 <- inFile2$datapath
  #     }
  # 
  #   res<- reliabilityCheckLibrary(dat1, dat2)
  # })
#     
#     # output$rellibcheck <- renderUI({
#     # 
#     #  paste("plot has been saved in directory")
#     #   # lib.stats <- librarychecking()
#     #   # paste(lib.stats)
#     #   
#     #   
#     #   
#     # })
#   
#   output$libcheckplot <- renderImage({
#      
#    #  
#    #  lib.labels = c('seed.library', 'extended.library')
#    # 
#      lib.stats <- librarychecking()
#    # 
#    #  sl = ceiling(max(lib.stats[,1])/1000)/ceiling(max(lib.stats[,2])/1000)
#    # 
#    #  pepnum.scaled = lib.stats[,2]*sl
#    # 
#    #  ymax = 1.1*max(c(lib.stats[,1], pepnum.scaled))
#    # 
#    # bp1 = barplot(lib.stats[,1], space=0.5, names.arg=lib.labels, ylim=c(0, ymax),
#    #                ylab='Number of proteins', axes=FALSE, col=gray(0.8), main='Libraries')
#    # 
#    #  at1 = axis(side=2)
#    # 
#    #  points(bp1, lib.stats[,2]*sl, type='b', lwd=3, col=gray(0.4) )
#    #  at2 = round(at1/(sl*1000))*1000
#    #  axis(side=4, at = at1, labels=at2)
#    #  mtext('Number of peptides', side=4, line=3)
#    # 
#    #  par(xpd=TRUE)
#    #  legend(bp1[1,1], ceiling(1.05*ymax) , bty='n', fill=gray(0.8), legend='protein')
#    #  legend(bp1[2,1], ceiling(1.05*ymax), bty='n', col=gray(0.4),  lty=1, legend='peptide', lwd=2)
#     width = 500
#     height = 700
#      
#      filepath <- getwd()
#     filepath2 <- paste0(filepath,"/plots_of_library_info.png")
#   list(src = filepath2, width = width, height = height)
#     
#     
#   })
#   
#   
#   
# #////////////////////////////////////////////////////////////////////////////////// 
# 
# #### Swath Reliability Check
#   
#   swathresults <- eventReactive(input$swathr,
#                                 {
#                                   swath1 <- input$seedswathfile
#                                   swath2 <- input$extendedswathfile
#                                   fdrpass <- input$fdrpass
#                                   fdrpeptide <- input$maxfdrpep
#                                   
#                                   if(!is.null(swath1) && !is.null(swath2)){
#                                     withProgress(message = "Performing FDR calculations", style = "notification", value = 0.1, {
#                                       for(i in 1:1) {
#                                         reliabilityCheckSwath(swath1, swath2, fdrpass, fdrpeptide)
#                                         incProgress(1)
#                                         Sys.sleep(0.5)
#                                       }
#                                     }
#                                     )
#                                   }
#                                   
#                                   # reliabilityCheckSwath(swath1, swath2, fdrpass, fdrpeptide)
#                                 })
#   
#   observeEvent(input$swathr, {
#     
#     shinyjs::show(id = "swathtables")
#     showTab(inputId = "swathrestables", target = "FDR Bins", select = TRUE)
#   })
#   
#   
#   output$fdrbinsseed <- DT::renderDataTable({
#     
#     swath <- swathresults()
# 
#       DT::datatable(data = swath[[1]], 
#                     options = list(pageLength = 50,
#                                    scrollX = T, 
#                                    scrollY = "500px", 
#                                    scrollCollapse = T, 
#                                    autoWidth = T,
#                                    scroller.loadingIndicator = T),
#                     caption = "FDR Bins")
#     
#   })
#   
#   
#/////////////////////////////////////////////////////////////////////////////////  
  
  ####### Help page 
  
  output$helppage <- renderText({
    
    filepath <- getwd()
    filepath2 <- paste0(filepath,"/www/help-app.html")
    helptext <- read_file(filepath2)
    HTML(helptext)
    
  })
  
  ####//////////////////////////////////////////////////////////////////////////////////
  
})   ##### Server Function Ends Here
