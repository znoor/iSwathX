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
library(RColorBrewer)
library(grid)
library(gridExtra)
library(pheatmap)

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
source("readReportFile.R")
source("processReport.R")
source("reportSummary.R")
source("multiCorrLibThree.R")
source("multiCorrLibFour.R")
source("buildSpectralLibThree.R")
source("buildSpectralLibFour.R")
# source("reliabilityCheckSwath.R")


shinyServer(function(input, output, session) {
  
  value <- reactiveVal(0)
  value1 <- reactiveVal(0)
  
  
  shinyjs::disable("apply")
  shinyjs::disable("tmultiapply")
  shinyjs::disable("fmultiapply")
  shinyjs::disable("cleanupdate")
  shinyjs::disable("inclength")
  shinyjs::disable("downloadseedinputLibs")
  shinyjs::disable("downloadextinputLibs")
  shinyjs::disable("report")
  shinyjs::disable("processFile")
  shinyjs::disable("downloadreport")
  # shinyjs::disable("wizdownloadoutputLib")
  # shinyjs::disable("autowizrun")
  #### Apply button which is activating the tabs
  
  observeEvent(input$apply, {
    
    shinyjs::show(id = "libreadpanel")
    showTab(inputId = "libraries", target = "Seed Library", select = TRUE)
  })
  
  observeEvent(input$report,{
    shinyjs::show(id = "reportreadpanel")
  })
  
  
  observeEvent(input$multireport,{
    shinyjs::show(id = "multireportpanel")
  })
  
  observeEvent(input$tmultiapply, {
    shinyjs::show(id = "multilibsummary")
  })
  
  observeEvent(input$fmultiapply, {
    shinyjs::show(id = "multilibsummary")
  })
  
  observeEvent(input$tmultiapply, {
    shinyjs::show(id = "multilibcompare")
  })
  
  observeEvent(input$fmultiapply, {
    shinyjs::show(id = "multilibcompare")
  })
  
  observeEvent(input$corcompute, {
    shinyjs::show(id = "corrstats")
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
                          "\n", " Unmodified Peptides = ", formatC(lib1sumdata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
                          " All Peptides = ", formatC(lib1sumdata[["modpeptides"]], format = "d", big.mark = ","), "\n",
                          " Transitions = ",formatC(lib1sumdata[["transitions"]], format = "d", big.mark = ",") ,  "\n", "\n",
                          
                          "External assay library contains \n",  " Proteins = ", formatC(lib2sumdata[["proteins"]], format = "d", big.mark = ","),
                          "\n", " Unmodified Peptides = ", formatC(lib2sumdata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
                          " All Peptides = ", formatC(lib2sumdata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
                          " Transitions = ", formatC(lib2sumdata[["transitions"]], format = "d", big.mark = ",") , "\n", "\n",
                          
                          "Retention time correlation between seed and external library = ", format(rtcorrelation[[4]], digits = 2) , "\n","\n",
                          
                          "Combined assay library contains \n", " Proteins = ", formatC(comblibdata2[["proteins"]], format = "d", big.mark = ","),
                          "\n", " Unmodified Peptides = ", formatC(comblibdata2[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
                          " All Peptides = ", formatC(comblibdata2[["modpeptides"]], format = "d", big.mark = ","),  "\n",
                          " Transitions = " , formatC(comblibdata2[["transitions"]], format = "d", big.mark = ","))
                   # }
                   # else paste(lib1sumdata[["proteins"]])
                 }) ########## have to set all lib summary codes
                 
               })
  
  observeEvent(input$autowizrun, 
               {
                 output$dvricor <- renderPlot({
                   
                   ricor <- autoricor()
                   
                   ionCorGS<-NULL; rm(ionCorGS);
                   data(ionCorGS)
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
    paste("Combined", input$autoseedlib$name,input$autoextlib$name , Sys.Date(), ".txt", sep = "_")
  },
  content = function (file) {
    
    data = autocombLib()
    data = data[[1]]
    
    # write.table(x = peakviewFormat(data), file = paste0(filename, "peakview" , ".csv"), sep = ",")
    
    # filename <-  paste("Combined", input$autoseedlib$name,input$autoextlib$name , sep = "_")
    
    if(input$autooutputlibformat == "PeakView"){
      write.table(x = peakviewFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$autooutputlibformat == "OpenSwath")
    {
      write.table(x = OswathFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$autooutputlibformat == "Skyline")
    {
      write.table(x = skylineFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else # spectronaut
    {
      write.table(x = spectronautFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
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
  ### Seed Library Read
  
  seedlibdata <- eventReactive(input$apply, {
    inFile <- input$seedlib
    
    if(is.null(inFile)) {
      
      return(NULL)
    } else
      
      lib1 <- readLibFile(inFile$datapath, input$libformats11, "spectrum", clean = FALSE)
    return(lib1)
  })
  
  ### External Library Read
  extlibdata <- eventReactive(input$apply,{
    inFile <- input$extlib
    
    if(is.null(inFile)){
      return(NULL)
    } else
      
      lib2 <- readLibFile(inFile$datapath, input$libformats12, "spectrum", clean = FALSE)
    return(lib2)
  })
  
  ### Report File Read
  
  reportfiledata <- eventReactive(input$report, {
    inFile <- input$inputreport
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  ### Multi Report File Read 1
  
  reportfiledata1 <- eventReactive(input$report, {
    inFile <- input$inputreport1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  ### Multi Report File Read 2
  
  reportfiledata2 <- eventReactive(input$report, {
    inFile <- input$inputreport2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  ### Multi Report File Read 3
  
  reportfiledata3 <- eventReactive(input$report, {
    inFile <- input$inputreport3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  #//////////////////////////////////////////////////////////////////////////////////
  
  #### Clean button for cleaning update
  
  seedcleanlibdata <- eventReactive(input$cleanupdate, {
    clib1 <- lib1clean()
    # lib1 <- seedlibdata()
    # clib1 <- cleanLib(lib1, TRUE, intensity.cutoff = input$intensity,
    #          conf.cutoff = input$confidence, prec.charge = input$precursorcharge, prod.charge = input$productcharge,
    #          nomod = input$modified, frag.number = input$fragmentnumber, nomc = input$misscleavage, enz = input$enzyme)
    return(clib1)
    # return(clib1)
  })
  extcleanlibdata <- eventReactive(input$cleanupdate,{
    clib2 <- lib2clean()
    # lib2 <- extlibdata()
    # clib2 <- cleanLib(lib2, TRUE, intensity.cutoff = input$intensity,
    #          conf.cutoff = input$confidence, prec.charge = input$precursorcharge, prod.charge = input$productcharge,
    #          nomod = input$modified, frag.number = input$fragmentnumber, nomc = input$misscleavage, enz = input$enzyme)
    return(clib2)
  })
  
  reportprocessdata <- eventReactive(input$processFile, {
    reportf <- reportprocess()
    return(reportf)
  })
  
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
               conf.cutoff = input$confidence, prec.charge = input$precursorcharge, prod.charge = input$productcharge,
               nomod = input$modified, frag.number = input$fragmentnumber, nomc = input$misscleavage, enz = input$enzyme)
    })
  
  
  observeEvent(input$cleanupdate, {
    newValue <- value() + 1     
    value(newValue)             
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
  
  
  #### Report Reading Reactive
  
  reportread <- reactive({
    inFile <- input$inputreport
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  
  #### Multi Report Reading Reactive 1
  
  reportread1 <- reactive({
    inFile <- input$inputreport1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  #### Multi Report Reading Reactive 2
  
  reportread2 <- reactive({
    inFile <- input$inputreport2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  
  #### Multi Report Reading Reactive 3
  
  reportread3 <- reactive({
    inFile <- input$inputreport3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      
      report <- readReportFile(inFile$datapath)
    return(report)
  })
  #/////////////////////////////////////////////////////////////////////////////////  
  
  #### External Library Cleaning reactive
  
  lib2clean <- reactive(
    {
      lib2 <- lib2read()
      cleanLib(lib2, TRUE, intensity.cutoff = input$intensity,
               conf.cutoff = input$confidence, prec.charge = input$precursorcharge, prod.charge = input$productcharge,
               nomod = input$modified, frag.number = input$fragmentnumber, nomc = input$misscleavage, enz = input$enzyme)
    })
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Report File processing reactive
  
  reportprocess <- reactive(
    {
      repp <- reportread()
      processReport(repp, spiked = input$refineText, avgDotP = input$averagedotp,
                    avgMass = input$averagemasserr, avgRT = input$averagert, avgQval = input$averageqval)
    })
  
  
  observeEvent(input$processFile, {
    newValue1 <- value1() + 1     
    value1(newValue1)             
  })
  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable apply button    
  output$readapply <- renderUI({
    
    inFile1 <- input$seedlib
    inFile2 <- input$extlib
    if(is.null(inFile1) && is.null(inFile2)) return()
    
    # if (is.null(lib1read()) && is.null(lib2read())) return()
    shinyjs::enable("apply")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable report apply button    
  output$reportapply <- renderUI({
    
    inFile<- input$inputreport
    if(is.null(inFile)) return()
    
    # if (is.null(lib1read()) && is.null(lib2read())) return()
    shinyjs::enable("report")
  })  
  
  #### Enable multireport apply button    
  output$multireportapply <- renderUI({
    
    inFile<- input$inputreport
    if(is.null(inFile)) return()
    
    # if (is.null(lib1read()) && is.null(lib2read())) return()
    shinyjs::enable("report")
  }) 
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable clean button    
  output$clean.lib <- renderUI({
    if (is.null(seedlibdata()) && is.null(extlibdata())) return()
    shinyjs::enable("cleanupdate")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable report process button    
  output$process.rep <- renderUI({
    if (is.null(reportread())) return()
    shinyjs::enable("processFile")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable downloadseed button    
  output$downloadseed <- renderUI({
    if (is.null(seedlibdata()) && is.null(extlibdata())) return()
    else if(!is.null(seedlibdata()))
      shinyjs::enable("downloadseedinputLibs")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable downloadext button    
  output$downloadext <- renderUI({
    if (is.null(seedlibdata()) && is.null(seedlibdata())) return()
    else if(!is.null(extlibdata()))
      shinyjs::enable("downloadextinputLibs")
  })  
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  #### Enable download report button    
  output$downloadrep <- renderUI({
    if (is.null(reportread())) return()
    else if(!is.null(reportread()))
      shinyjs::enable("downloadreport")
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
  
  observeEvent(input$apply, {
    output$seedlibcontents <- DT::renderDataTable({
      withProgress(message = "Generating seed library table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          data = lib1read()
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
      
    })
  })
  
  observeEvent(input$cleanupdate,{
    output$seedlibcontents <- DT::renderDataTable({
      withProgress(message = "Generating seed library table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          data = lib1clean()
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
      
    })
  })
  
  ############# Report Data Table
  
  observeEvent(input$report, {
    output$reportContents <- DT::renderDataTable({
      withProgress(message = "Generating report file table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          data <- reportread()
          incProgress(0.1, detail = "reading report file")
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
                    caption = input$inputreport$name
      )
      
    })
  }) 
  
  
  observeEvent(input$processFile, {
    output$reportContents <- DT::renderDataTable({
      withProgress(message = "Generating report file table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          data <- reportprocess()
          incProgress(0.1, detail = "reading report file")
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
                    caption = input$inputreport$name
      )
      
    })
  }) 
  
  #/////////////////////////////////////////////////////////////////////////////////  
  
  #### External Library Data Table 
  
  
  observeEvent(input$apply, {
    output$extlibcontents <- DT::renderDataTable({
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
      
    })
  })
  
  observeEvent(input$cleanupdate,{
    output$extlibcontents <- DT::renderDataTable({
      if (!is.null(extcleanlibdata())){
        withProgress(message = "Generating external library table", style = "notification", value = 0.1, {
          for(i in 1:1) {
            data = extcleanlibdata()
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
      }
    })
  })
  
  #///////////////////////////////////////////////////////////////////////////////// 
  
  #### Comb Lib Plot
  
  output$combplot <- renderPlot({
    
    if(!is.null(seedlibdata()) && !is.null(extlibdata())) {
      withProgress(message = "Generating combined library plot", style = "notification", value = 0.1, {
        for(i in 1:1) {
          data = combdatalib()
          data <- data[[1]]
          # x<-aggregate(lib1data$uniprot_id, by=list(Protein=lib1data$uniprot_id, Peptide=lib1data$stripped_sequence), 
          #              FUN=function(x){length(x)})
          # y<-aggregate(x$Peptide, by=list(Protein=x$Protein), FUN=function(x){length(x)})
          
          data2 <- data[!duplicated(data$modification_sequence),]
          x <- aggregate(x = data2$modification_sequence, by = list(Protein = data2$uniprot_id, charge = data2$prec_z),
                         FUN = function(x) {length(x)})
          x$charge <- as.character(x$charge)
          breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
          bwidth1 <- breaks1[3] - breaks1[1]
          
          x2 <- aggregate(x = data$modification_sequence, by = list(Peptide = data$modification_sequence, 
                                                                    frg_type = data$frg_type), FUN = function(x2) {length(x2)})
          breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
          bwidth2 <- breaks2[3] - breaks2[1]
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
      
      p1 <- ggplot(data = x, aes(x = x, fill = charge)) +
        geom_histogram(binwidth = bwidth1) +
        labs(x ="Number of Peptides", y = "Number of Proteins") +
        # scale_y_continuous(name = "Frequency") +
        scale_fill_brewer("Precursor Charge", palette = "YlOrRd")+
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      
      p2 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
        geom_histogram(binwidth = bwidth2) +
        labs(x ="Number of Ions", y = "Number of Peptides") +
        scale_fill_brewer("Fragment Type", palette = "PuRd") +
        # scale_fill_manual(values = c("#eeba30", "#ae0001")) +
        # scale_y_continuous(name = "Frequency") +
        # scale_fill_gradient("Frequency", low = "blue", high = "red")+
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      
      p <- annotate_figure(ggarrange(p1, p2, ncol = 2), top = "Frequency distribution of Peptides per Protein and Ions per Peptide")
      print(p)
      
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs")
      ggsave("Combined_library_plot.png", width = 8, height = 5, path = filepath2)
      
    } else return()
    
  })
  
  #/////////////////////////////////////////////////////////////////////////////////  
  
  #### Seed Lib Plot
  
  output$lib1plot <- renderPlot({
    
    if(value() == 1)
      lib1data = seedcleanlibdata()
    else
      lib1data = seedlibdata()
    
    if (!is.null(lib1data)){
      withProgress(message = "Generating seed library plot", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib1data2 <- lib1data[!duplicated(lib1data$modification_sequence),]
          x <- aggregate(x = lib1data2$modification_sequence, by = list(Protein = lib1data2$uniprot_id, 
                                                                        charge = lib1data2$prec_z), FUN = function(x) {length(x)})
          x$charge <- as.character(x$charge)
          breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
          bwidth1 <- breaks1[3] - breaks1[1]
          
          x2 <- aggregate(x = lib1data$modification_sequence, by = list(Peptide = lib1data$modification_sequence, 
                                                                        frg_type = lib1data$frg_type), FUN = function(x2) {length(x2)})
          breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
          bwidth2 <- breaks2[2] - breaks2[1]
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
      p1 <- ggplot(data = x, aes(x = x, fill = charge)) +
        geom_histogram(binwidth = bwidth1) +
        labs(x ="Number of Peptides", y = "Number of Proteins") +
        # scale_y_continuous(name = "Frequency") +
        scale_fill_brewer("Precursor Charge", palette = "YlOrRd")+
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      p2 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
        geom_histogram(binwidth = bwidth2) +
        labs(x ="Number of Ions", y = "Number of Peptides") +
        scale_fill_brewer("Fragment Type", palette = "PuRd") +
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      p <- annotate_figure(ggarrange(p1, p2, ncol = 2), top = "Frequency distribution of Peptides per Protein and Ions per Peptide", fig.lab.face = "bold")
      print(p)
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs")
      ggsave("Seed_library_plot.png", width = 8, height = 5, path = filepath2)
    }
  })
  
  
  #/////////////////////////////////////////////////////////////////////////////////  
  
  #### External Library Plot  
  
  output$lib2plot <- renderPlot({
    
    if(value() == 1)
      lib2data = extcleanlibdata()
    else
      lib2data = extlibdata() 
    
    if (!is.null(lib2data)){
      withProgress(message = "Generating external library plot", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib2data2 <- lib2data[!duplicated(lib2data$modification_sequence),]
          x <- aggregate(x = lib2data2$modification_sequence, by = list(Protein = lib2data2$uniprot_id, 
                                                                        charge = lib2data2$prec_z), FUN = function(x) {length(x)})
          x$charge <- as.character(x$charge)
          breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
          bwidth1 <- breaks1[3] - breaks1[1]
          x2 <- aggregate(x = lib2data$modification_sequence, by = list(Peptide = lib2data$modification_sequence, 
                                                                        frg_type = lib2data$frg_type), FUN = function(x2) {length(x2)})
          breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
          bwidth2 <- breaks2[2] - breaks2[1]
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
      p1 <- ggplot(data = x, aes(x = x, fill = charge)) +
        geom_histogram(binwidth = bwidth1) +
        labs(x ="Number of Peptides", y = "Number of Proteins") +
        scale_fill_brewer("Precursor Charge", palette = "YlOrRd")+
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      
      p2 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
        geom_histogram(binwidth = bwidth2) +
        labs(x ="Number of Ions", y = "Number of Peptides") +
        scale_fill_brewer("Fragment type", palette = "PuRd") +
        theme_classic() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face= "bold", size = 12)) 
      
      p <- annotate_figure(ggarrange(p1, p2, ncol = 2), "Frequency distribution of Peptides per Protein and Ions per Peptide", fig.lab.face = "bold")
      print(p)
      
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs")
      ggsave("External_library_plot.png", width = 8, height = 5, path = filepath2)
    }
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
  {paste("Library_seed",input$seedlib$name , Sys.Date(), ".txt", sep = "_")},
  
  content= function(file)
  {
    if(value() == 1)
      data = seedcleanlibdata()
    else
      data = seedlibdata()
    
    if(input$inputseedlibformat == "PeakView"){
      write.table(x = peakviewFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else if(input$inputseedlibformat == "OpenSwath")
    {
      write.table(x = OswathFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else if(input$inputseedlibformat == "Skyline")
    {
      write.table(x = skylineFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else # spectronaut
    {
      write.table(x = spectronautFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    }
    
  }
  )
  
  output$downloadextinputLibs <- downloadHandler(filename = function()
  {paste("Library_ext",input$extlib$name , Sys.Date(), ".txt", sep = "_")},
  
  content= function(file)
  {
    if(value() == 1)
      data = extcleanlibdata()
    else
      data = extlibdata()
    
    if(input$inputextlibformat == "PeakView"){
      write.table(x = peakviewFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else if(input$inputextlibformat == "OpenSwath")
    {
      write.table(x = OswathFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else if(input$inputextlibformat == "Skyline")
    {
      write.table(x = skylineFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
      
    } else # spectronaut
    {
      write.table(x = spectronautFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    }
    
  }
  )
  
  
  ### download report file
  
  output$downloadreport <- downloadHandler(filename = function()
  {paste("Report",input$inputreport$name , Sys.Date(), ".csv", sep = "_")},
  
  content= function(file)
  {
    if(value1() == 1)
      data = reportprocess()
    else
      data = reportread()
    
    write.csv(x = data, file = file, row.names = F, na = " ")
    
  }
  )
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  
  ####  Output / Combined Libraries Download 
  output$downloadoutputLib <- downloadHandler(filename = function () 
  {paste("Combined", input$seedlib$name,input$extlib$name , Sys.Date(), ".txt", sep = "_")
  },
  content = function (file) {
    
    data = combdatalib()
    data = data[[1]]
    
    if(input$outputlibformat == "PeakView"){
      write.table(x = peakviewFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$outputlibformat == "OpenSwath")
    {
      write.table(x = OswathFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$outputlibformat == "Skyline")
    {
      write.table(x = skylineFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else # spectronaut
    {
      write.table(x = spectronautFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    }
    
  }
  )
  
  #/////////////////////////////////////////////////////////////////////////////////
  
  
  #### Available Report File and Reps
  
  output$reportfilename <- renderUI({
    if (is.null(reportfiledata())) return()
    else
      paste(input$inputreport$name)
  })
  
  output$reps <- renderUI({
    if (is.null(reportfiledata())) return()
    else
    { 
      dat <- reportread()
      reps <- length(colnames(dat[, str_which(colnames(dat), 
                                              pattern = "Library.Dot.Product")]))
      if(is.null(reps)) {
        reps <- length(colnames(dat[, str_which(colnames(dat),
                                                pattern = "Peptide.Retention.Time")]))
        
        if(is.null(reps))
          reps <- length(colnames(dat[, str_which(colnames(dat), 
                                                  pattern = "Average.Mass.Error.PPM")]))
        if(is.null(reps))
          reps <- "null_paste"
      }
      
      paste(reps)
    }
  })
  
  #//////////////////////////////////////////////////////////////////////////////////    
  
  #### Available Libraries
  
  output$libseedfilesrt <- renderUI({
    if (is.null(seedlibdata())) return()
    else
      paste(input$seedlib$name)
  })
  
  output$libextfilesrt <- renderUI({
    if (is.null(extlibdata())) return()
    else
      paste(input$extlib$name)
  })
  
  #//////////////////////////////////////////////////////////////////////////////////
  
  #### Available Libraries
  
  output$libseedfilesint <- renderUI({
    if (is.null(seedlibdata())) return()
    else
      paste(input$seedlib$name)
  })
  
  output$libextfilesint <- renderUI({
    if (is.null(extlibdata())) return()
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
  
  
  observeEvent(input$cleanupdate, {
    output$seedlibsummary <- renderText({
      seeddata <- seedcleanlibdata()
      if(!is.null(seeddata)) {
        summarydata <- libSummary(seeddata)
        
        paste0("Seed assay library contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  observeEvent(input$apply, {
    output$seedlibsummary <- renderText({
      seeddata <- seedlibdata()
      if(!is.null(seeddata)) {
        summarydata <- libSummary(seeddata)
        
        paste0("Seed assay library contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  
  #//////////////////////////////////////////////////////////////////////////////////
  
  
  #### Report File Summary
  
  observeEvent(input$report, {
    output$reportSummary <- renderText({
      
      data = reportread()
      
      if(!is.null(data)) {
        repsummarydata <- reportSummary(data)
        paste0("Report File contains \n", " Proteins = ", formatC(repsummarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(repsummarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
               " Modified Peptides = ", formatC(repsummarydata[["modpeptides"]], format = "d", big.mark = ","))
      }
    })
  }) 
  
  
  observeEvent(input$processFile, {
    output$reportSummary <- renderText({
      
      data = reportprocess()
      
      if(!is.null(data)) {
        repsummarydata <- reportSummary(data)
        paste0("Report File contains \n", " Proteins = ", formatC(repsummarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(repsummarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
               " Modified Peptides = ", formatC(repsummarydata[["modpeptides"]], format = "d", big.mark = ","))
      }
    })
  }) 
  
  #////////////////////////////////////////////////////////////////////////////////// 
  
  #### External Library Summary
  
  observeEvent(input$cleanupdate, {
    output$extlibsummary <- renderText({
      extdata <- extcleanlibdata()
      if(!is.null(extdata)){
        summarydata <- libSummary(extdata)
        paste0("External assay library contains \n", " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","), 
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","), "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  observeEvent(input$apply, {
    output$extlibsummary <- renderText({
      extdata <- extlibdata()
      if(!is.null(extdata)){
        summarydata <- libSummary(extdata)
        paste0("External assay library contains \n", " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","), 
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","), "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  
  #////////////////////////////////////////////////////////////////////////////////// 
  
  #### Combined Library Summary
  
  output$comblibsummary <- renderText({
    if(!is.null(seedlibdata()) && !is.null(extlibdata())){
      data = combdatalib()
    }
    summarydata <- libSummary(data[[1]])
    paste0("Combined library contains \n", " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
           "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
           " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","), "\n",
           " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
    
  })
  
  #////////////////////////////////////////////////////////////////////////////////// 
  
  #### Compute RT Time reactive function
  
  ### by time tables/ graphs
  
  rtcordt <- eventReactive(input$calRT, {
    # dat1 <- lib1clean()
    # dat2 <- lib2clean()
    
    if(value() == 1)
    {dat1 = seedcleanlibdata()
    dat2 = extcleanlibdata()}
    else
    {dat1 = seedlibdata()
    dat2 = extlibdata()}
    
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
    # dat1 <- lib1clean()
    # dat2 <- lib2clean()
    
    if(value() == 1)
    {dat1 = seedcleanlibdata()
    dat2 = extcleanlibdata()}
    else
    {dat1 = seedlibdata()
    dat2 = extlibdata()}
    
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
    
    # dat1 <- lib1clean()
    # dat2 <- lib2clean()
    
    if(value() == 1)
    {dat1 = seedcleanlibdata()
    dat2 = extcleanlibdata()}
    else
    {dat1 = seedlibdata()
    dat2 = extlibdata()}
    
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
    # dat1 <- lib1clean()
    # dat2 <- lib2clean()
    
    if(value() == 1)
    {dat1 = lib1clean()
    dat2 = lib2clean()}
    else
    {dat1 = lib1read()
    dat2 = lib2read()}
    
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
    
    # dat1 <- lib1clean()
    # dat2 <- lib2clean()
    
    if(value() == 1)
    {dat1 = lib1clean()
    dat2 = lib2clean()}
    else
    {dat1 = lib1read()
    dat2 = lib2read()}
    
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
    
    # dat1 <- lib1clean()
    # dat2 <- lib2clean() 
    
    if(value() == 1)
    {dat1 = seedcleanlibdata()
    dat2 = extcleanlibdata()}
    else
    {dat1 = seedlibdata()
    dat2 = extlibdata()}
    
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
      # dat1 <- lib1clean()
      # dat2 <- lib2clean()
      
      if(value() == 1)
      {dat1 = seedcleanlibdata()
      dat2 = extcleanlibdata()}
      else
      {dat1 = seedlibdata()
      dat2 = extlibdata()}
      
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
    
    if(input$uselibs == "clean") {
      dat1 <- lib1clean()
      dat2 <- lib2clean()
    } else {
      
      dat1 <- lib1read()
      dat2 <- lib2read() }
    
    # if(value() == 1)
    # {dat1 = seedcleanlibdata()
    # dat2 = extcleanlibdata()}
    # else
    # {dat1 = seedlibdata()
    # dat2 = extlibdata()}
    
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
    if(!is.null(seedlibdata()) && !is.null(extlibdata())){
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
    
    output$dvpepvenn <- renderPlot({
      data = autocombLib()
      data2 <- data[[2]]
      if(is.null(data2))
        return()
      else print(data2$pvennpep)
    })
  })
  
  observeEvent(input$combrun, {
    
    output$dvpepvenn <- renderPlot({
      data = combdatalib()
      data2 <- data[[2]]
      if(is.null(data2))
        return()
      else print(data2$pvennpep)
    })
  })
  
  observeEvent(input$autowizrun, {
    
    output$dvprotvenn <- renderPlot({
      data = autocombLib()
      data2 <- data[[2]]
      if(is.null(data2))
        return()
      else print(data2$pvennprot)
    })
  })
  
  observeEvent(input$combrun, {
    
    output$dvprotvenn <- renderPlot({
      data = combdatalib()
      data2 <- data[[2]]
      if(is.null(data2))
        return()
      else print(data2$pvennprot)
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
    } , deleteFile = F)
  })
  
  observeEvent(input$combrun, {
    
    output$dvlibinfo <- renderImage({
      # data = combdatalib()
      # dat1 = lib1clean()
      # reliabilityCheckLibrary(dat1, data)
      filepath <- getwd()
      filepath2 <- paste0(filepath,"/graphs/plots_of_library_info.png")
      list(src = filepath2, width = "100%", height = "100%")
    } , deleteFile = F)
  })
  
  
  # ///////////////////////////////////////////////////////////////////////////////////
  
  #### Report DotP Graphs
  
  output$dotpPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          cols <- dat[, str_which(colnames(dat), pattern = "Library.Dot.Product")]
          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Peptide", colnames(cols))
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          dat <- gather(dat, Replicates, Dot.Product, -Peptide)
          dat$Dot.Product <- as.numeric(dat$Dot.Product) 
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dat[dat == 0] <- NA
          
          dp <- ggplot(data = dat, aes(y = Dot.Product, x = "", fill = Replicates)) +
            geom_boxplot(alpha = 1,
                         notch = TRUE, notchwidth = 0.8,
                         outlier.colour = "red", outlier.fill = "red", outlier.size = 1) +
            # scale_fill_manual(values=c("#d4d0d0", "#bebbbb", "#a9a6a6", "#949191",  "#7f7c7c")) +
            # scale_fill_brewer(palette="Dark2") +
            scale_fill_manual(values = getPalette(length(unique(dat$Replicates)))) +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(face = "bold"),
                  legend.justification = "center") +
            xlab(paste(input$inputreport$name)) +
            ylab("Dot Product")
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  #### Report Intensity Distribution Graph
  
  output$intensityPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          cols <- dat[, str_which(colnames(dat), pattern = "Normalized.Area")]
          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Peptide", colnames(cols))
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          dat <- gather(dat, Replicates, Area, -Peptide)
          dat$Area <- log(as.numeric(format(dat$Area, digits = 7)), base = 2)
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dp <- ggplot(data = dat, aes(y = Area, x = "", fill = Replicates)) +
            geom_boxplot(alpha = 1,
                         notch = TRUE, notchwidth = 0.8,
                         outlier.colour = "red", outlier.fill = "red", outlier.size = 1) +
            # scale_fill_manual(values=c("#d4d0d0", "#bebbbb", "#a9a6a6", "#949191",  "#7f7c7c")) +
            # scale_fill_brewer(palette="Dark2") +
            scale_fill_manual(values = getPalette(length(unique(dat$Replicates)))) +
            facet_wrap(~Replicates, drop = FALSE, scales = "free_y") +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(face = "bold"),
                  legend.justification = "center",
                  legend.position = "None") +
            xlab(paste(input$inputreport$name)) +
            ylab("Log Intensity")
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  ###### Mass Error Plots
  
  output$massEPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          cols <- dat[, str_which(colnames(dat), pattern = "Mass.Error.PPM")]
          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Peptide", colnames(cols))
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          dat <- gather(dat, Replicates, Mass.Error, -Peptide)
          dat$Mass.Error <- as.numeric(dat$Mass.Error) 
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dp <- ggplot(data = dat, aes(x = Mass.Error, fill = Replicates)) +
            geom_histogram(binwidth = 0.4, color = "black",
                           # fill = "#ff4945", 
                           alpha = 0.9) +
            labs(x = "Mass Error (PPM)", y = "Count") +
            # scale_fill_manual(values=c("#d4d0d0", "#bebbbb", "#a9a6a6", "#949191",  "#7f7c7c")) +
            # scale_fill_brewer(palette="Dark2") +
            scale_fill_manual(values = getPalette(length(unique(dat$Replicates)))) +
            facet_wrap(~Replicates, drop = FALSE, scales = "free_y") +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(face = "bold"),
                  legend.justification = "center", 
                  legend.position = "None") +
            xlab("Mass Error (PPM)") +
            ylab("Count")
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  ##### CV Plots
  
  output$cvPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          dat <- mutate(dat, "CV.Normalized" = as.character(dat$Cv.Total.Area.Normalized))
          dat$CV.Normalized <- str_replace(dat$CV.Normalized, pattern = "%", replacement = "")
          dat$CV.Normalized <- round(as.numeric(dat$CV.Normalized))
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dp <- ggplot(data = dat, aes(y = CV.Normalized, x = "")) +
            geom_violin(alpha = 0.9, color = "black",
                        trim = F)+
            # scale_fill_manual(values = c("#ff4945")) +
            scale_fill_brewer(palette="Dark2") +
            geom_boxplot(width = 0.1, fill = "white") +
            theme(legend.position = "none") +
            labs(x = paste(input$inputreport$name), y = "Normalized CVs (%)") +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  aspect.ratio = 1
            )
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  ###### Q-Value Plots
  
  output$qvalPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          cols <- dat[, str_which(colnames(dat), pattern = "Detection.Q.Value")]
          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Peptide", colnames(cols))
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          # dat$Detection.Q.Value <- as.numeric(format(dat$Area, digits = 4))
          dat <- gather(dat, Replicates, Q_Value, -Peptide)
          dat$Q_Value <- as.numeric(dat$Q_Value)
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dp <- ggplot(data = dat) +
            geom_density(aes(x = Q_Value, y = ..scaled.., fill = Replicates), alpha=.8, stat = "density", position = "identity", linetype = "solid") +
            # scale_fill_manual(values=c("#d4d0d0", "#bebbbb", "#a9a6a6", "#949191",  "#7f7c7c")) +
            # scale_fill_brewer(palette="Dark2") +
            scale_fill_manual(values = getPalette(length(unique(dat$Replicates)))) +
            facet_wrap(~Replicates, drop = FALSE, scales = "free_y") +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"),
                  axis.title.y = element_text(color="black", size=12, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 11),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(face = "bold"),
                  legend.justification = "center",
                  legend.position = "None") +
            # xlab("Cattle") +
            ylab("Density")
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  ###### DDA - DIA Retention time Plots
  
  
  report_lib_file <- reactive({
    
    inFile <- input$report_lib
    
    if(is.null(inFile)) {
      
      return(NULL)
    } else
      
      readLibFile(inFile$datapath, input$report_lib_format, "spectrum", clean = FALSE)
  })
  
  observeEvent(input$drtplot, {
  output$drtPlot <- renderPlot({
    if(!is.null(reportread()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]
          
          dat_rt <- dat[, c("Protein.Name", "Peptide", "Modified.Sequence",  "Average.Measured.Retention.Time")]
          
          lib_rt <- report_lib_file()
          lib_rt <- lib_rt[!duplicated(lib_rt$modification_sequence), c("modification_sequence", "RT_detected")]
          
          comm_pep <- intersect(dat$Modified.Sequence, lib_rt$modification_sequence)
          
          dat_rt <- dat_rt[dat_rt$Modified.Sequence %in% comm_pep,]
          lib_rt <- lib_rt[lib_rt$modification_sequence %in% comm_pep,]
          
          rt_cor_data <- merge(dat_rt, lib_rt, by.x = 3, by.y = 1, all = TRUE)
          
          rt_cor_data$Average.Measured.Retention.Time <- as.numeric(rt_cor_data$Average.Measured.Retention.Time)
          
          
          drt <- ggplot(rt_cor_data, aes(x = Average.Measured.Retention.Time, y = RT_detected)) +
            geom_point(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
            annotate("text", x = 20, y = 50,
                     label = deparse(bquote(italic(R)^2 ==. (format(
                       round(cor(rt_cor_data$Average.Measured.Retention.Time, rt_cor_data$RT_detected), digits = 2), digits = 2)))),
                     color = "black", parse = T) +
            stat_smooth(method = "auto", col = "red", se = T, size = 0.8)+
            scale_y_continuous("DDA RT (min)", limits = c(min(rt_cor_data$RT_detected), max(rt_cor_data$RT_detected)),
                               expand = c(0,0)) +
            scale_x_continuous("DIA RT (min)", limits = c(min(rt_cor_data$Average.Measured.Retention.Time), max(rt_cor_data$Average.Measured.Retention.Time)),
                               expand = c(0,0)) +
          theme(legend.position = "none") +
            # labs(x = paste(input$inputreport$name), y = "Normalized CVs (%)") +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  aspect.ratio = 1
            )
          
          print(drt)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  })
  
  
  # ###### PROTEIN INTENSITIES TABLE AND PLOT
  # 
  output$intensity_table <- DT::renderDataTable({

    if(!is.null(reportread())) {
      withProgress(message = "Generating Table", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]

          cols <- dat[, str_which(colnames(dat), pattern = "Normalized.Area")]
          cols <- apply(cols, 2, as.numeric)
          cols <- as.data.frame(cols)

          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Protein.Name","Peptide", "Modified.Sequence")
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          dat <- cbind(dat, cols)

          dat_int <- dplyr::group_by(dat, Protein.Name)
          dat_int2 <- dat_int %>% summarise_at(vars(colnames(cols)), funs(sum))

          dat_int2[is.na(dat_int2)] <- 0


          

          incProgress(0.1, detail = "Generating Table")
          Sys.sleep(0.25)
        }
      })
    
    DT::datatable(data = dat_int2,
                  options = list(pageLength = 10,
                                 scrollX = T,
                                 scrollY = "500px",
                                 scrollCollapse = T,
                                 autoWidth = T,
                                 scroller.loadingIndicator = T)
                  # caption = input$seedlib$name
    )
    }

  })
  
  
  # 
  output$intensity_plot <- renderPlot({

    if(!is.null(reportread())) {
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat <- reportread()
          dat <- dat[!duplicated(dat$Modified.Sequence),]

          cols <- dat[, str_which(colnames(dat), pattern = "Normalized.Area")]
          cols <- apply(cols, 2, as.numeric)
          cols <- as.data.frame(cols)

          # cols <- cbind(cols, dat[, str_which(colnames(dat), pattern = "Peptide")])
          cols_names <- c("Protein.Name","Peptide", "Modified.Sequence")
          # colnames(dat) <- c("Peptides", "Rep01", "Rep02", "Rep03", "Rep04", "Rep05")
          dat <- dat[, cols_names]
          dat <- cbind(dat, cols)

          dat_int <- dplyr::group_by(dat, Protein.Name)
          dat_int2 <- dat_int %>% summarise_at(vars(colnames(cols)), funs(sum))

          dat_int2[is.na(dat_int2)] <- 0


          dat_int3 <- dat_int2[, "Protein.Name"]
          dat_int3 <- cbind(dat_int3, apply(dat_int2[, str_which(colnames(dat_int2), pattern = "Normalized.Area")], 2, log2))

          dat_int3 <- dat_int3 %>%
            mutate(Proteins = paste0("Prot", 1:length(dat_int3$Protein.Name)))

          dat_int3 <- dat_int3[, -c(1)]
          colnames(dat_int3) <- c(paste0("Rep", 1: (ncol(dat_int3)-1)), "Proteins")

          dat_int4 <- gather(data = dat_int3, key = Replicates, value = Intensity, -c(ncol(dat_int3)))
          
          textcol <- "grey40"

          rep_int <- ggplot(data = dat_int4, mapping = aes(x = Replicates, y = Proteins, fill = Intensity)) +
            geom_tile(colour="gray",size=0.0005) +
            labs(x="",y="", title = "Protein intensities comparison among replicates")+
            scale_y_discrete(expand=c(0,0))+
            scale_fill_distiller(palette = "Spectral") +
            theme_grey(base_size=10)+
            theme(legend.position="right",legend.direction="vertical",
                  legend.title=element_text(colour=textcol, face = "bold"),
                  legend.margin=margin(grid::unit(0,"cm")),
                  legend.text=element_text(colour=textcol,size=7,face="bold"),
                  legend.key.height=grid::unit(0.8,"cm"),
                  legend.key.width=grid::unit(0.2,"cm"),
                  axis.text.x=element_text(size=10,colour=textcol),
                  axis.text.y=element_text(vjust=0.2,colour=textcol),
                  axis.ticks=element_line(size=0.4),
                  plot.background=element_blank(),
                  panel.border=element_blank(),
                  # plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
                  plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
                  # aspect.ratio = 0.5
                  )

          
          print(rep_int)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
    }

  })
  
  # ///////////////////////////////////////////////////////////////////////////////////
  
  #### Report Statistics Graphs
  
  output$report_stats <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1_prot <- dat1[!duplicated(dat1$Protein.Name),"Protein.Name"]
          dat1_pep <- dat1[!duplicated(dat1$Modified.Sequence),"Modified.Sequence"]
          
          dat2_prot <- dat2[!duplicated(dat2$Protein.Name),"Protein.Name"]
          dat2_pep <- dat2[!duplicated(dat2$Modified.Sequence),"Modified.Sequence"]
          
          dat3_prot <- dat3[!duplicated(dat3$Protein.Name),"Protein.Name"]
          dat3_pep <- dat3[!duplicated(dat3$Modified.Sequence),"Modified.Sequence"]
          
          extractions_data <- data.frame("Proteins" = c(length(dat1_prot), length(dat2_prot), length(dat3_prot)),
                                         "Peptides" = c(length(dat1_pep), length(dat2_pep), length(dat3_pep)),
                                         "Data" = c(input$dataset1, input$dataset2, input$dataset3))
          
          rep_prot <- ggplot(data=extractions_data, aes(x = Data, y = Proteins, fill = Data)) +
            geom_bar(stat = "identity", alpha = 0.8, position=position_dodge(), width = 0.9) +
            geom_text(aes(label=Proteins), vjust=0, color="black", fontface = "bold",
            position = position_dodge(0.9),
            size=5) +
            # scale_y_continuous(breaks=seq(0, 300, 50))+
            #scale_fill_brewer(palette="Paired") +
            scale_fill_manual(labels = c(input$dataset1, input$dataset2, input$dataset3), values=c("#ff4945", "#2ac940", "#75a3e7"))+
            # scale_x_discrete(labels = c(input$inputreport1$name, input$inputreport2$name, input$inputreport3$name)) +
            theme(axis.title.x = element_text(color="black", size=13, face="bold"),
                  axis.title.y = element_text(color="black", size=13, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
                  # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                  #                                 colour = "white"), 
                  # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                  #                                 colour = "gray"),
                  axis.text.y = element_text(color="black", size = 12, face = "bold"),
                  axis.text.x = element_text(color="black", size = 12, face = "bold"),
                  legend.text = element_text(size = 12, face = "bold"),
                  legend.title = element_text(size = 14, face = "bold"),
                  legend.justification = "center",
                  legend.position = "None")
          
          rep_pep <- ggplot(data=extractions_data, aes(x = Data, y = Peptides, fill = Data)) +
            geom_bar(stat = "identity", alpha = 0.8, position=position_dodge(), width = 0.9) +
            geom_text(aes(label=Peptides), vjust=0, color="black", fontface = "bold",
            position = position_dodge(0.9),
            size=5) +
            # scale_y_continuous(breaks=seq(0, 300, 50))+
            #scale_fill_brewer(palette="Paired") +
            scale_fill_manual(labels = c(input$dataset1, input$dataset2, input$dataset3), values=c("#ff4945", "#2ac940", "#75a3e7"))+
            # scale_x_discrete(labels = c(input$inputreport1$name, input$inputreport2$name, input$inputreport3$name)) +
            theme(axis.title.x = element_text(color="black", size=13, face="bold"),
                  axis.title.y = element_text(color="black", size=13, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
                  # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                  #                                 colour = "white"), 
                  # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                  #                                 colour = "gray"),
                  axis.text.y = element_text(color="black", size = 12, face = "bold"),
                  axis.text.x = element_text(color="black", size = 12, face = "bold"),
                  legend.text = element_text(size = 12, face = "bold"),
                  legend.title = element_text(size = 14, face = "bold"),
                  legend.justification = "center",
                  legend.position = "None")
          
          rep_stats <- ggarrange(rep_prot, rep_pep)
          
          print(rep_stats)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # ///////////////////////////////////////////////////////////////////////////////////    
  
  #### Multireport venn diagram Rep 1 - Rep 2
  # c("#ff4945", "#2ac940", "#75a3e7")
  output$rep1_rep2 <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1_prot <- dat1[!duplicated(dat1$Protein.Name),"Protein.Name"]
          dat1_pep <- dat1[!duplicated(dat1$Modified.Sequence),"Modified.Sequence"]
          
          dat2_prot <- dat2[!duplicated(dat2$Protein.Name),"Protein.Name"]
          dat2_pep <- dat2[!duplicated(dat2$Modified.Sequence),"Modified.Sequence"]
          
          dat3_prot <- dat3[!duplicated(dat3$Protein.Name),"Protein.Name"]
          dat3_pep <- dat3[!duplicated(dat3$Modified.Sequence),"Modified.Sequence"]
          
          d1_d2_pep <- venn.diagram(list(dat1_pep, dat2_pep),
                                  category.names = c(input$dataset1, input$dataset2),
                                  resolution = 500,
                                  # height = 1000,
                                  # width = 1000,
                                  cat.default.pos = "outer",
                                  lwd = 2,
                                  fill = c("#ff4945", "#2ac940"),
                                  cat.pos = c(-170, 170),
                                  alpha = c(0.7, 0.7), 
                                  cex = 1.5,
                                  cat.fontface = 2,
                                  lty =2, 
                                  cat.fontfamily = "Arial",
                                  fontfamily = "Arial",
                                  # fontfamily ="sans",
                                  # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Peptides.tiff"
                                  filename = NULL)
          d1_d2_pro <- venn.diagram(list(dat1_prot, dat2_prot),
                                  category.names = c(input$dataset1, input$dataset2),
                                  resolution = 500,
                                  # height = 1000,
                                  # width = 1000,
                                  cat.default.pos = "outer",
                                  lwd = 2,
                                  fill = c("#ff4945", "#2ac940"),
                                  cat.pos = c(-170, 170),
                                  alpha = c(0.7, 0.7), 
                                  cex = 1.5,
                                  cat.fontface = 2,
                                  lty =2, 
                                  cat.fontfamily = "Arial",
                                  fontfamily = "Arial",
                                  # fontfamily ="sans",
                                  # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Proteins.tiff"
                                  filename = NULL)
          
          vennpro12 <- gTree(children = d1_d2_pro)
          vennpro12 <- as_ggplot(vennpro12)
          
          vennpep12 <- gTree(children = d1_d2_pep)
          vennpep12 <- as_ggplot(vennpep12)
          
          venn12 <- ggarrange(vennpro12, vennpep12, ncol = 2, nrow = 1, labels = c("A. Proteins", "B. Peptides"), font.label = list(size = 15), hjust = -0.1)
          print(venn12)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # /////////////////////////////////////////////////////////////////////////////////// 
  
  #### Multireport venn diagram Rep 2 - Rep 3
  # c("#ff4945", "#2ac940", "#75a3e7")
  output$rep2_rep3 <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1_prot <- dat1[!duplicated(dat1$Protein.Name),"Protein.Name"]
          dat1_pep <- dat1[!duplicated(dat1$Modified.Sequence),"Modified.Sequence"]
          
          dat2_prot <- dat2[!duplicated(dat2$Protein.Name),"Protein.Name"]
          dat2_pep <- dat2[!duplicated(dat2$Modified.Sequence),"Modified.Sequence"]
          
          dat3_prot <- dat3[!duplicated(dat3$Protein.Name),"Protein.Name"]
          dat3_pep <- dat3[!duplicated(dat3$Modified.Sequence),"Modified.Sequence"]
          
          d2_d3_pep <- venn.diagram(list(dat2_pep, dat3_pep),
                                    category.names = c(input$dataset2, input$dataset3),
                                    resolution = 500,
                                    # height = 1000,
                                    # width = 1000,
                                    cat.default.pos = "outer",
                                    lwd = 2,
                                    fill = c("#2ac940", "#75a3e7"),
                                    cat.pos = c(-170, 170),
                                    alpha = c(0.7, 0.7), 
                                    cex = 1.5,
                                    cat.fontface = 2,
                                    lty =2, 
                                    cat.fontfamily = "Arial",
                                    fontfamily = "Arial",
                                    # fontfamily ="sans",
                                    # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Peptides.tiff"
                                    filename = NULL)
          d2_d3_pro <- venn.diagram(list(dat2_prot, dat3_prot),
                                    category.names = c(input$dataset2, input$dataset3),
                                    resolution = 500,
                                    # height = 1000,
                                    # width = 1000,
                                    cat.default.pos = "outer",
                                    lwd = 2,
                                    fill = c("#2ac940", "#75a3e7"),
                                    cat.pos = c(-170, 170),
                                    alpha = c(0.7, 0.7), 
                                    cex = 1.5,
                                    cat.fontface = 2,
                                    lty =2, 
                                    cat.fontfamily = "Arial",
                                    fontfamily = "Arial",
                                    # fontfamily ="sans",
                                    # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Proteins.tiff"
                                    filename = NULL)
          
          vennpro23 <- gTree(children = d2_d3_pro)
          vennpro23 <- as_ggplot(vennpro23)
          
          vennpep23 <- gTree(children = d2_d3_pep)
          vennpep23 <- as_ggplot(vennpep23)
          
          venn23 <- ggarrange(vennpro23, vennpep23, ncol = 2, nrow = 1, labels = c("A. Proteins", "B. Peptides"), font.label = list(size = 15), hjust = -0.1)
          print(venn23)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # /////////////////////////////////////////////////////////////////////////////////// 
  
  #### Multireport venn diagram Rep 1 - Rep 3
  # c("#ff4945", "#2ac940", "#75a3e7")
  output$rep1_rep3 <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1_prot <- dat1[!duplicated(dat1$Protein.Name),"Protein.Name"]
          dat1_pep <- dat1[!duplicated(dat1$Modified.Sequence),"Modified.Sequence"]
          
          dat2_prot <- dat2[!duplicated(dat2$Protein.Name),"Protein.Name"]
          dat2_pep <- dat2[!duplicated(dat2$Modified.Sequence),"Modified.Sequence"]
          
          dat3_prot <- dat3[!duplicated(dat3$Protein.Name),"Protein.Name"]
          dat3_pep <- dat3[!duplicated(dat3$Modified.Sequence),"Modified.Sequence"]
          
          d1_d3_pep <- venn.diagram(list(dat1_pep, dat3_pep),
                                    category.names = c(input$dataset1, input$dataset3),
                                    resolution = 500,
                                    # height = 1000,
                                    # width = 1000,
                                    cat.default.pos = "outer",
                                    lwd = 2,
                                    fill = c("#ff4945", "#75a3e7"),
                                    cat.pos = c(-170, 170),
                                    alpha = c(0.7, 0.7), 
                                    cex = 1.5,
                                    cat.fontface = 2,
                                    lty =2, 
                                    cat.fontfamily = "Arial",
                                    fontfamily = "Arial",
                                    # fontfamily ="sans",
                                    # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Peptides.tiff"
                                    filename = NULL)
          d1_d3_pro <- venn.diagram(list(dat1_prot, dat3_prot),
                                    category.names = c(input$dataset1, input$dataset3),
                                    resolution = 500,
                                    # height = 1000,
                                    # width = 1000,
                                    cat.default.pos = "outer",
                                    lwd = 2,
                                    fill = c("#ff4945", "#75a3e7"),
                                    cat.pos = c(-170, 170),
                                    alpha = c(0.7, 0.7), 
                                    cex = 1.5,
                                    cat.fontface = 2,
                                    lty =2, 
                                    cat.fontfamily = "Arial",
                                    fontfamily = "Arial",
                                    # fontfamily ="sans",
                                    # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/Cattle_SELheep_Proteins.tiff"
                                    filename = NULL)
          
          vennpro13 <- gTree(children = d1_d3_pro)
          vennpro13 <- as_ggplot(vennpro13)
          
          vennpep13 <- gTree(children = d1_d3_pep)
          vennpep13 <- as_ggplot(vennpep13)
          
          venn13 <- ggarrange(vennpro13, vennpep13, ncol = 2, nrow = 1, labels = c("A. Proteins", "B. Peptides"), font.label = list(size = 15), hjust = -0.1)
          print(venn13)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # ///////////////////////////////////////////////////////////////////////////////////
  
  #### Multireport venn diagram Rep 1 - Rep 2 - Rep 3
  # c("#ff4945", "#2ac940", "#75a3e7")
  output$rep1_rep2_rep3 <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1_prot <- dat1[!duplicated(dat1$Protein.Name),"Protein.Name"]
          dat1_pep <- dat1[!duplicated(dat1$Modified.Sequence),"Modified.Sequence"]
          
          dat2_prot <- dat2[!duplicated(dat2$Protein.Name),"Protein.Name"]
          dat2_pep <- dat2[!duplicated(dat2$Modified.Sequence),"Modified.Sequence"]
          
          dat3_prot <- dat3[!duplicated(dat3$Protein.Name),"Protein.Name"]
          dat3_pep <- dat3[!duplicated(dat3$Modified.Sequence),"Modified.Sequence"]
          
          d123_pep <- venn.diagram(
            x = list(dat1_pep, dat2_pep, dat3_pep),
            category.names = c(input$dataset1, input$dataset2, input$dataset3),
            # height = 2000,
            # width = 2000,
            resolution = 500,
            lwd = 2,
            lty = 2,
            fill = c("#ff4945", "#2ac940", "#75a3e7"),
            cat.cex = 1,
            cex = 1.5,
            cat.fontfamily = "Arial",
            fontfamily = "Arial",
            cat.fontface = 2,
            cat.default.pos = "outer",
            cat.pos = c(-36, 30, 135),
            cat.dist = c(0.055, 0.055, 0.085),
            rotation = 1,
            alpha = rep(0.7, 3),
            # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/PBL_Peptides.tiff"
            filename = NULL
          )
          
          d123_pro <- venn.diagram(
            x = list(dat1_prot, dat2_prot, dat3_prot),
            category.names = c(input$dataset1, input$dataset2, input$dataset3),
            resolution = 500,
            # height = 1000,
            # width = 1000,
            lwd = 2,
            lty = 2,
            fill = c("#ff4945", "#2ac940", "#75a3e7"),
            cat.cex = 1,
            cex = 1.5,
            cat.fontfamily = "Arial",
            fontfamily = "Arial",
            cat.fontface = 2,
            cat.default.pos = "outer",
            cat.pos = c(-36, 30, 135),
            cat.dist = c(0.055, 0.055, 0.085),
            rotation = 1,
            alpha = rep(0.7, 3),
            # filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/allC_Proteins.tiff"
            filename = NULL
          )
          
          vennpro123 <- gTree(children = d123_pro)
          vennpro123 <- as_ggplot(vennpro123)
          
          vennpep123 <- gTree(children = d123_pep)
          vennpep123 <- as_ggplot(vennpep123)
          
          venn123 <- ggarrange(vennpro123, vennpep123, ncol = 2, nrow = 1, labels = c("A. Proteins", "B. Peptides"), font.label = list(size = 15), hjust = -0.1)
          print(venn123)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # ///////////////////////////////////////////////////////////////////////////////////
  
  #### Retention Time Correlation

  output$multirtCorPlot <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          
          dat1 <- dat1[!duplicated(dat1$Modified.Sequence),]
          dat1$Average.Measured.Retention.Time <- as.numeric(dat1$Average.Measured.Retention.Time)
          
          dat2 <- dat2[!duplicated(dat2$Modified.Sequence),]
          dat2$Average.Measured.Retention.Time <- as.numeric(dat2$Average.Measured.Retention.Time)
          
          dat3 <- dat3[!duplicated(dat3$Modified.Sequence),]
          dat3$Average.Measured.Retention.Time <- as.numeric(dat3$Average.Measured.Retention.Time)
          
          dat1_dat2_peptides <- intersect(dat1$Peptide, dat2$Peptide)
          
          dat1_dat3_peptides <- intersect(dat1$Peptide, dat3$Peptide)
          
          dat2_dat3_peptides <- intersect(dat2$Peptide, dat3$Peptide)
          
          
          d1_d2_rt <- dat1[dat1$Peptide %in% dat1_dat2_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          d2_d1_rt <- dat2[dat2$Peptide %in% dat1_dat2_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          
          d1_d3_rt <- dat1[dat1$Peptide %in% dat1_dat3_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          d3_d1_rt <- dat3[dat3$Peptide %in% dat1_dat3_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          
          d2_d3_rt <- dat2[dat2$Peptide %in% dat2_dat3_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          d3_d2_rt <- dat3[dat3$Peptide %in% dat2_dat3_peptides, c("Peptide", "Average.Measured.Retention.Time")]
          
          
          dat1_dat2 <- merge(d1_d2_rt, d2_d1_rt, by = 1)
          dat1_dat3 <- merge(d1_d3_rt, d3_d1_rt, by = 1)
          dat2_dat3 <- merge(d2_d3_rt, d3_d2_rt, by = 1)
          
          p1 <- ggplot(dat1_dat2, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
            geom_point(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
            annotate("text", x = 25, y = 45, 
                     label = deparse(bquote(italic(R)^2 ==. (format(
                       cor(as.numeric(dat1_dat2$Average.Measured.Retention.Time.x), as.numeric(dat1_dat2$Average.Measured.Retention.Time.y)), digits = 3)))),
                     color = "black", parse = T) +
            stat_smooth(method = "auto", col = "red", se = T, size = 0.8)+
            scale_y_continuous(input$dataset1, limits = c(min(dat1_dat2$Average.Measured.Retention.Time.y), max(dat1_dat2$Average.Measured.Retention.Time.y)),
                               expand = c(0,0)) +
            scale_x_continuous(input$dataset2, limits = c(min(dat1_dat2$Average.Measured.Retention.Time.x), max(dat1_dat2$Average.Measured.Retention.Time.x)),
                               expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
                  axis.text = element_text(color="black", size = 12, face = "bold"),
                  axis.title.x = element_text(color="black", size=12, face="bold"),
                  axis.title.y = element_text(color="black", size=12, face="bold"),
                  aspect.ratio = 1
            )
          
          p2 <- ggplot(dat1_dat3, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
            geom_point(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
            annotate("text", x = 25, y = 45, 
                     label = deparse(bquote(italic(R)^2 ==. (format(
                       cor(as.numeric(dat1_dat3$Average.Measured.Retention.Time.x), as.numeric(dat1_dat3$Average.Measured.Retention.Time.y)), digits = 3)))),
                     color = "black", parse = T) +
            stat_smooth(method = "auto", col = "red", se = T, size = 0.8)+
            scale_y_continuous(input$dataset1, limits = c(min(dat1_dat3$Average.Measured.Retention.Time.y), max(dat1_dat3$Average.Measured.Retention.Time.y)),
                               expand = c(0,0)) +
            scale_x_continuous(input$dataset3, limits = c(min(dat1_dat3$Average.Measured.Retention.Time.x), max(dat1_dat3$Average.Measured.Retention.Time.x)),
                               expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
                  axis.text = element_text(color="black", size = 12, face = "bold"),
                  axis.title.x = element_text(color="black", size=12, face="bold"),
                  axis.title.y = element_text(color="black", size=12, face="bold"),
                  aspect.ratio = 1
            )
          
          p3 <- ggplot(dat2_dat3, aes(x = Average.Measured.Retention.Time.x, y = Average.Measured.Retention.Time.y)) +
            geom_point(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
            annotate("text", x = 25, y = 45, 
                     label = deparse(bquote(italic(R)^2 ==. (format(
                       cor(as.numeric(dat2_dat3$Average.Measured.Retention.Time.x), as.numeric(dat2_dat3$Average.Measured.Retention.Time.y)), digits = 3)))),
                     color = "black", parse = T) +
            stat_smooth(method = "auto", col = "red", se = T, size = 0.8)+
            scale_y_continuous(input$dataset2, limits = c(min(dat2_dat3$Average.Measured.Retention.Time.y), max(dat2_dat3$Average.Measured.Retention.Time.y)),
                               expand = c(0,0)) +
            scale_x_continuous(input$dataset3, limits = c(min(dat2_dat3$Average.Measured.Retention.Time.x), max(dat2_dat3$Average.Measured.Retention.Time.x)),
                               expand = c(0,0)) +
            theme(panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
                  axis.text = element_text(color="black", size = 12, face = "bold"),
                  axis.title.x = element_text(color="black", size=12, face="bold"),
                  axis.title.y = element_text(color="black", size=12, face="bold"),
                  aspect.ratio = 1
            )
          
          pnum <- ggarrange(p1, p2, p3, ncol = 3)
          pnum
          
          print(pnum)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  # ///////////////////////////////////////////////////////////////////////////////////
  
  #### Transition Ratio Correlation
  
  # output$multitransCorPlot <- renderPlot({
  #   if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
  #     withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
  #       for (i in 1:1) {
  #         dat1 <- reportread1()
  #         dat2 <- reportread2()
  #         dat3 <- reportread3()
  #         
  #         dat1 <- dat1 %>%
  #           select("Protein.Name", "Peptide", "Modified.Sequence", "Fragment.Ion", 
  #                  str_which(colnames(dat1), pattern = "Peak.Rank"), str_which(colnames(dat1), fixed(pattern ="Area"))) 
  #         dat2 <- dat2 %>%
  #           select("Protein.Name", "Peptide", "Modified.Sequence", "Fragment.Ion", 
  #                  str_which(colnames(dat2), pattern = "Peak.Rank"), str_which(colnames(dat2), fixed(pattern ="Area"))) 
  #         dat3 <- dat3 %>%
  #           select("Protein.Name", "Peptide", "Modified.Sequence", "Fragment.Ion", 
  #                  str_which(colnames(dat3), pattern = "Peak.Rank"), str_which(colnames(dat3), fixed(pattern ="Area"))) 
  # 
  # 
  #         dat1[, str_which(colnames(dat1), pattern = "Area")] <- apply(m <- (dat1[, str_which(colnames(dat1), pattern = "Area")]), 2, as.numeric)
  #         dat1[, str_which(colnames(dat1), pattern = "Peak.Rank")] <- apply(m <- (dat1[, str_which(colnames(dat1), pattern = "Peak.Rank")]), 2, as.numeric)
  #         dat1 <- na.omit(dat1)
  # 
  #         dat2[, str_which(colnames(dat2), pattern = "Area")] <- apply(m <- (dat2[, str_which(colnames(dat2), pattern = "Area")]), 2, as.numeric)
  #         dat2[, str_which(colnames(dat2), pattern = "Peak.Rank")] <- apply(m <- (dat2[, str_which(colnames(dat2), pattern = "Peak.Rank")]), 2, as.numeric)
  #         dat2 <- na.omit(dat2)
  # 
  #         dat3[, str_which(colnames(dat3), pattern = "Area")] <- apply(m <- (dat3[, str_which(colnames(dat3), pattern = "Area")]), 2, as.numeric)
  #         dat3[, str_which(colnames(dat3), pattern = "Peak.Rank")] <- apply(m <- (dat3[, str_which(colnames(dat3), pattern = "Peak.Rank")]), 2, as.numeric)
  #         dat3 <- na.omit(dat3)
  # 
  #         ### Some rearrangment
  # 
  #         dat1_tr <- dat1[!duplicated(dat1$Modified.Sequence), c("Protein.Name", "Peptide", "Modified.Sequence")]
  #         dat2_tr <- dat2[!duplicated(dat2$Modified.Sequence), c("Protein.Name", "Peptide", "Modified.Sequence")]
  #         dat3_tr <- dat3[!duplicated(dat3$Modified.Sequence), c("Protein.Name", "Peptide", "Modified.Sequence")]
  # 
  # 
  #         area_col <- dat1[, str_which(colnames(dat1), pattern = "Area")]
  #         area_col_names <- colnames(area_col)
  #         
  # dat1_tr <- dat1_tr %>%
  # mutate(RP1.Transition.Ratio = dat1[dat1$RP1.Peak.Rank == 1, area_col_names[1]] / 
  #          dat1[dat1$RP1.Peak.Rank == 2, area_col_names[1]])
  # dat1_tr <- dat1_tr %>%
  #   mutate(RP2.Transition.Ratio = dat1[dat1$RP2.Peak.Rank == 1, area_col_names[2]] / 
  #            dat1[dat1$RP2.Peak.Rank == 2, area_col_names[2]])
  # dat1_tr <- dat1_tr %>%
  #   mutate(RP3.Transition.Ratio = dat1[dat1$RP3.Peak.Rank == 1, area_col_names[3]] / 
  #            dat1[dat1$RP3.Peak.Rank == 2, area_col_names[3]])
  # dat1_tr <- dat1_tr %>%
  #   mutate(RP4.Transition.Ratio = dat1[dat1$RP4.Peak.Rank == 1, area_col_names[4]] / 
  #            dat1[dat1$RP4.Peak.Rank == 2, area_col_names[4]])
  # dat1_tr <- dat1_tr %>%
  #   mutate(RP5.Transition.Ratio = dat1[dat1$RP5.Peak.Rank == 1, area_col_names[5]] / 
  #            dat1[dat1$RP5.Peak.Rank == 2, area_col_names[5]])
  # 
  # ################################################################
  # 
  # area_col <- dat2[, str_which(colnames(dat2), pattern = "Area")]
  # area_col_names <- colnames(area_col)
  # 
  # dat2_tr <- dat2_tr %>%
  #   mutate(RP1.Transition.Ratio = dat2[dat2$RP1.Peak.Rank == 1, area_col_names[1]] / 
  #            dat2[dat2$RP1.Peak.Rank == 2, area_col_names[1]])
  # dat2_tr <- dat2_tr %>%
  #   mutate(RP2.Transition.Ratio = dat2[dat2$RP2.Peak.Rank == 1, area_col_names[2]] / 
  #            dat2[dat2$RP2.Peak.Rank == 2, area_col_names[2]])
  # dat2_tr <- dat2_tr %>%
  #   mutate(RP3.Transition.Ratio = dat2[dat2$RP3.Peak.Rank == 1, area_col_names[3]] / 
  #            dat2[dat2$RP3.Peak.Rank == 2, area_col_names[3]])
  # dat2_tr <- dat2_tr %>%
  #   mutate(RP4.Transition.Ratio = dat2[dat2$RP4.Peak.Rank == 1, area_col_names[4]] / 
  #            dat2[dat2$RP4.Peak.Rank == 2, area_col_names[4]])
  # dat2_tr <- dat2_tr %>%
  #   mutate(RP5.Transition.Ratio = dat2[dat2$RP5.Peak.Rank == 1, area_col_names[5]] / 
  #            dat2[dat2$RP5.Peak.Rank == 2, area_col_names[5]])
  # 
  # ################################################################
  # 
  # area_col <- dat3[, str_which(colnames(dat3), pattern = "Area")]
  # area_col_names <- colnames(area_col)
  # 
  # dat3_tr <- dat3_tr %>%
  #   mutate(RP1.Transition.Ratio = dat3[dat3$RP1.Peak.Rank == 1, area_col_names[1]] / 
  #            dat3[dat3$RP1.Peak.Rank == 2, area_col_names[1]])
  # dat3_tr <- dat3_tr %>%
  #   mutate(RP2.Transition.Ratio = dat3[dat3$RP2.Peak.Rank == 1, area_col_names[2]] / 
  #            dat3[dat3$RP2.Peak.Rank == 2, area_col_names[2]])
  # dat3_tr <- dat3_tr %>%
  #   mutate(RP3.Transition.Ratio = dat3[dat3$RP3.Peak.Rank == 1, area_col_names[3]] / 
  #            dat3[dat3$RP3.Peak.Rank == 2, area_col_names[3]])
  # dat3_tr <- dat3_tr %>%
  #   mutate(RP4.Transition.Ratio = dat3[dat3$RP4.Peak.Rank == 1, area_col_names[4]] / 
  #            dat3[dat3$RP4.Peak.Rank == 2, area_col_names[4]])
  # dat3_tr <- dat3_tr %>%
  #   mutate(RP5.Transition.Ratio = dat3[dat3$RP5.Peak.Rank == 1, area_col_names[5]] / 
  #            dat3[dat3$RP5.Peak.Rank == 2, area_col_names[5]])
  # 
  # ################################################################
  # 
  # dat1_tr <- dat1_tr %>%
  #   mutate(Average.Transition.Ratio = apply(m <- (dat1_tr[, c("RP1.Transition.Ratio", "RP2.Transition.Ratio", 
  #                                                             "RP3.Transition.Ratio", "RP4.Transition.Ratio", "RP5.Transition.Ratio")]), 1, mean, na.rm = T))
  # dat2_tr <- dat2_tr %>%
  #   mutate(Average.Transition.Ratio = apply(m <- (dat2_tr[, c("RP1.Transition.Ratio", "RP2.Transition.Ratio", 
  #                                                             "RP3.Transition.Ratio", "RP4.Transition.Ratio", "RP5.Transition.Ratio")]), 1, mean, na.rm = T))
  # dat3_tr <- dat3_tr %>%
  #   mutate(Average.Transition.Ratio = apply(m <- (dat3_tr[, c("RP1.Transition.Ratio", "RP2.Transition.Ratio", 
  #                                                             "RP3.Transition.Ratio", "RP4.Transition.Ratio", "RP5.Transition.Ratio")]), 1, mean, na.rm = T))
  # 
  # 
  # dat1_dat2_pep <- intersect(dat1_tr$Modified.Sequence, dat2_tr$Modified.Sequence)
  # 
  # dat1_dat3_pep <- intersect(dat1_tr$Modified.Sequence, dat3_tr$Modified.Sequence)
  # 
  # dat2_dat3_pep <- intersect(dat2_tr$Modified.Sequence, dat3_tr$Modified.Sequence)
  # 
  # 
  # dat1_dat2_dat3 <- intersect(intersect(dat1_tr$Modified.Sequence, dat2_tr$Modified.Sequence), dat3_tr$Modified.Sequence)
  # 
  # dat1_dat2_tr <- dat1_tr[dat1_tr$Modified.Sequence %in% dat1_dat2_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # dat2_dat1_tr <- dat2_tr[dat2_tr$Modified.Sequence %in% dat1_dat2_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # 
  # dat1_dat3_tr <- dat1_tr[dat1_tr$Modified.Sequence %in% dat1_dat3_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # dat3_dat1_tr <- dat3_tr[dat3_tr$Modified.Sequence %in% dat1_dat3_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # 
  # dat2_dat3_tr <- dat2_tr[dat2_tr$Modified.Sequence %in% dat2_dat3_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # dat3_dat2_tr <- dat3_tr[dat3_tr$Modified.Sequence %in% dat2_dat3_pep, c("Protein.Name", "Modified.Sequence", "Average.Transition.Ratio")]
  # 
  # 
  # 
  # dat11_dat22_tr <- merge(dat1_dat2_tr, dat2_dat1_tr, by = "Modified.Sequence")
  # dat22_dat33_tr <- merge(dat2_dat3_tr, dat3_dat2_tr, by = "Modified.Sequence")
  # dat11_dat33_tr <- merge(dat1_dat3_tr, dat3_dat1_tr, by = "Modified.Sequence")
  # 
  # 
  # dat11_dat22_tr <- dat11_dat22_tr %>%
  #   filter(Average.Transition.Ratio.x <= 20) %>%
  #   filter(Average.Transition.Ratio.y <= 20)
  # dat22_dat33_tr <- dat22_dat33_tr %>%
  #   filter(Average.Transition.Ratio.x <= 20) %>%
  #   filter(Average.Transition.Ratio.y <= 20)
  # dat11_dat33_tr <- dat11_dat33_tr %>%
  #   filter(Average.Transition.Ratio.x <= 20) %>%
  #   filter(Average.Transition.Ratio.y <= 20)
  # 
  # 
  #         ##### Some plot
  # tr1 <- ggplot(dat11_dat22_tr, aes(x = Average.Transition.Ratio.x, y = Average.Transition.Ratio.y)) +
  #   geom_jitter(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
  #   annotate("text", x = 6, y = 18, 
  #            label = deparse(bquote(italic(R)^2 ==. (format(
  #              cor(dat11_dat22_tr$Average.Transition.Ratio.x, dat11_dat22_tr$Average.Transition.Ratio.y), digits = 2)))),
  #            color = "black", parse = T) +
  #   stat_smooth(method = "auto", col = "red", se = T, size = 1)+
  #   scale_y_continuous("Sheep", breaks=seq(0, 20, 3)) +
  #   scale_x_continuous("Cattle", breaks=seq(0, 20, 3))+
  #   theme(panel.background = element_rect(fill = "white",
  #                                         colour = "white",
  #                                         size = 0.5, linetype = "solid"),
  #         panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
  #         axis.text = element_text(color="black",size = 13, face = "bold"),
  #         axis.title.x = element_text(color="black", size=12, face="bold"),
  #         axis.title.y = element_text(color="black", size=12, face="bold") 
  #   )
  # 
  # tr2 <- ggplot(dat11_dat33_tr, aes(x = Average.Transition.Ratio.x, y = Average.Transition.Ratio.y)) +
  #   geom_jitter(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
  #   annotate("text", x = 6, y = 13.5, 
  #            label = deparse(bquote(italic(R)^2 ==. (format(
  #              cor(dat11_dat33_tr$Average.Transition.Ratio.x, dat11_dat33_tr$Average.Transition.Ratio.y), digits = 2)))),
  #            color = "black", parse = T) +
  #   stat_smooth(method = "auto", col = "red", se = T, size = 1)+
  #   scale_y_continuous("Giraffe", breaks=seq(0, 15, 3)) +
  #   scale_x_continuous("Cattle", breaks=seq(0, 21, 3))+
  #   theme(panel.background = element_rect(fill = "white",
  #                                         colour = "white",
  #                                         size = 0.5, linetype = "solid"),
  #         panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
  #         axis.text = element_text(color="black",size = 13, face = "bold"),
  #         axis.title.x = element_text(color="black", size=12, face="bold"),
  #         axis.title.y = element_text(color="black", size=12, face="bold")
  #   )
  # 
  # tr3 <- ggplot(dat11_dat33_tr, aes(x = Average.Transition.Ratio.x, y = Average.Transition.Ratio.y)) +
  #   geom_jitter(alpha = 0.7, col = "black", fill = "black", size  = 1.5, shape = 15) +
  #   annotate("text", x = 6, y = 13.5, 
  #            label = deparse(bquote(italic(R)^2 ==. (format(
  #              cor(dat11_dat33_tr$Average.Transition.Ratio.x, dat11_dat33_tr$Average.Transition.Ratio.y), digits = 2)))),
  #            color = "black", parse = T) +
  #   stat_smooth(method = "auto", col = "red", se = T, size = 1)+
  #   scale_y_continuous("Sheep", breaks=seq(0, 20, 3)) +
  #   scale_x_continuous("Giraffe", breaks=seq(0, 20, 3))+
  #   theme(panel.background = element_rect(fill = "white",
  #                                         colour = "white",
  #                                         size = 0.5, linetype = "solid"),
  #         panel.border = element_blank(), axis.line = element_line(linetype = "solid"),
  #         axis.text = element_text(color="black",size = 13, face = "bold"),
  #         axis.title.x = element_text(color="black", size=12, face="bold"),
  #         axis.title.y = element_text(color="black", size=12, face="bold") 
  #   )
  # 
  # tr <- ggarrange(tr1, tr2, tr3, ncol = 3, nrow = 1, labels = "AUTO", hjust = -1.2, vjust = 1.7)
  #         
  # tr
  #         print(tr)
  # 
  #         incProgress(0.1, detail = "plotting")
  #         Sys.sleep(0.25)
  #       }
  #     })
  # })
  # 
  # 
  # ///////////////////////////////////////////////////////////////////////////////////
  
  
  
  ##### Multi CV Plots
  
  output$multiCVPlot <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()
          
          dat1 <- dat1[!duplicated(dat1$Modified.Sequence),]
          dat1 <- mutate(dat1, "CV.Normalized" = as.character(dat1$Cv.Total.Area))
          dat1$CV.Normalized <- str_replace(dat1$CV.Normalized, pattern = "%", replacement = "")
          dat1$CV.Normalized <- round(as.numeric(dat1$CV.Normalized))
          
          dat2 <- dat2[!duplicated(dat2$Modified.Sequence),]
          dat2 <- mutate(dat2, "CV.Normalized" = as.character(dat2$Cv.Total.Area))
          dat2$CV.Normalized <- str_replace(dat2$CV.Normalized, pattern = "%", replacement = "")
          dat2$CV.Normalized <- round(as.numeric(dat2$CV.Normalized))
          
          dat3 <- dat3[!duplicated(dat3$Modified.Sequence),]
          dat3 <- mutate(dat3, "CV.Normalized" = as.character(dat3$Cv.Total.Area))
          dat3$CV.Normalized <- str_replace(dat3$CV.Normalized, pattern = "%", replacement = "")
          dat3$CV.Normalized <- round(as.numeric(dat3$CV.Normalized))
          
          
          dat1 <- mutate(dat1, "Data" = input$dataset1)
          dat2 <- mutate(dat2, "Data" = input$dataset2)
          dat3 <- mutate(dat3, "Data" = input$dataset3)
          
          
          
          dat1 <- dat1[!duplicated(dat1$Modified.Sequence),]
          dat2 <- dat2[!duplicated(dat2$Modified.Sequence),]
          dat3 <- dat3[!duplicated(dat3$Modified.Sequence),]
          
         dat1 <- dat1[, c(1:10, 56, 57)]
         dat2 <- dat2[, c(1:10, 56, 57)]
         dat3 <- dat3[, c(1:10, 56, 57)]
          
          
          pv_dat <- rbind(dat1, dat2)
          pv_dat <- rbind(pv_dat, dat3)
          
          pv_dat <- pv_dat %>%
            filter(CV.Normalized <= 100)
          
          ## some plot
          getPalette = colorRampPalette(brewer.pal(9, "Dark2"))
          
          dp <- ggplot(data = pv_dat, aes(y = CV.Normalized, x = Data, fill = Data)) +
            geom_violin(alpha = 0.9, color = "black",
                        trim = F)+
            # scale_fill_manual(values = c("#ff4945")) +
            scale_fill_brewer(palette="Dark2") +
            geom_boxplot(width = 0.1, fill = "white") +
            theme(legend.position = "none") +
            labs(x = paste(input$inputreport$name), y = "Normalized CVs (%)") +
            theme(axis.title.x = element_text(color="black", size=16, face="bold"),
                  axis.title.y = element_text(color="black", size=16, face="bold"), 
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.5, linetype = "solid"),
                  panel.border = element_blank(), axis.line = element_line(linetype = "solid"), 
                  axis.text = element_text(color="black",size = 14),
                  aspect.ratio = 1
            )
          
          print(dp)
          
          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
      })
  })
  
  
  
  #### Protein Intensities Correlation
  
  output$multiprotIntPlot <- renderPlot({
    if(!is.null(reportread1()) & !is.null(reportread2()) & !is.null(reportread3()))
      withProgress(message = "Generating Plot", style = "notification", value = 0.1, {
        for (i in 1:1) {
          dat1 <- reportread1()
          dat2 <- reportread2()
          dat3 <- reportread3()

          dat1 <- dat1[!duplicated(dat1$Modified.Sequence),]
          dat2 <- dat2[!duplicated(dat2$Modified.Sequence),]
          dat3 <- dat3[!duplicated(dat3$Modified.Sequence),]

          dat1_dat2_comm_pep <- intersect(dat1$Modified.Sequence, dat2$Modified.Sequence)

          dat1_dat3_comm_pep <- intersect(dat1$Modified.Sequence, dat3$Modified.Sequence)

          dat2_dat3_comm_pep <- intersect(dat2$Modified.Sequence, dat3$Modified.Sequence)

          dat1_dat2_dat3_comm_pep <- intersect(intersect(dat1_dat2_comm_pep, dat1_dat3_comm_pep), dat2_dat3_comm_pep)

          # 
          dat1_int <- dat1[dat1$Peptide %in% dat1_dat2_dat3_comm_pep, ]
          dat2_int <- dat2[dat2$Peptide %in% dat1_dat2_dat3_comm_pep, ]
          dat3_int <- dat3[dat3$Peptide %in% dat1_dat2_dat3_comm_pep,]
          # 
          dat1_int2 <- dat1_int[order(dat1_int$Peptide),]
          dat2_int2 <- dat2_int[order(dat2_int$Peptide),]
          dat3_int2 <- dat3_int[order(dat3_int$Peptide),]
          # 
          # 
          dat1_int2[, str_which(colnames(dat1_int2), pattern = "Normalized.Area")] <- apply(m <- (dat1_int2[, str_which(colnames(dat1_int2), pattern = "Normalized.Area")]), 2, as.numeric)

          agg_dat1 <- dplyr::group_by(dat1_int2, Protein.Name)
          agg_dat1 <- agg_dat1 %>% summarise_at(vars(ends_with("Normalized.Area")),sum)

          dat1 <- merge(dat1_int2, agg_dat1, by = "Protein.Name", all = TRUE)
          dat1 <- dat1[order(dat1$Peptide),]

          dat1 <- dat1[!duplicated(dat1$Protein.Name),]

          dat1[is.na(dat1)] <- 0
          # ##
          # 
          dat2_int2[, str_which(colnames(dat2_int2), pattern = "Normalized.Area")] <- apply(m <- (dat2_int2[, str_which(colnames(dat2_int2), pattern = "Normalized.Area")]), 2, as.numeric)

          agg_dat2 <- dplyr::group_by(dat2_int2, Protein.Name)
          agg_dat2 <- agg_dat2 %>% summarise_at(vars(ends_with("Normalized.Area")),sum)

          dat2 <- merge(dat2_int2, agg_dat2, by = "Protein.Name", all = TRUE)
          dat2 <- dat2[order(dat2$Peptide),]

          dat2 <- dat2[!duplicated(dat2$Protein.Name),]

          dat2[is.na(dat2)] <- 0
          # ##
          # 
          dat3_int2[, str_which(colnames(dat3_int2), pattern = "Normalized.Area")] <- apply(m <- (dat3_int2[, str_which(colnames(dat3_int2), pattern = "Normalized.Area")]), 2, as.numeric)

          agg_dat3 <- dplyr::group_by(dat3_int2, Protein.Name)
          agg_dat3 <- agg_dat3 %>% summarise_at(vars(ends_with("Normalized.Area")),sum)

          dat3 <- merge(dat3_int2, agg_dat3, by = "Protein.Name", all = TRUE)
          dat3 <- dat3[order(dat3$Peptide),]

          dat3 <- dat3[!duplicated(dat3$Protein.Name),]

          dat3[is.na(dat3)] <- 0
          # ####
          # 
          # 
          dat1 <- dat1 %>%
            mutate(Average.Intensity.dat1 = apply(m <- (dat1[, str_which(colnames(dat1), pattern = "Normalized.Area.y")]), 1, mean, na.rm = TRUE))
          dat2 <- dat2 %>%
            mutate(Average.Intensity.dat2 = apply(m <- (dat2[, str_which(colnames(dat2), pattern = "Normalized.Area.y")]), 1, mean, na.rm = TRUE))
          dat3 <- dat3 %>%
            mutate(Average.Intensity.dat3 = apply(m <- (dat3[, str_which(colnames(dat3), pattern = "Normalized.Area.y")]), 1, mean, na.rm = TRUE))

          # 
          dat1_dat2 <- merge(dat1, dat2, by = 1)
          dat1_dat2_dat3 <- merge(dat1_dat2, dat3, by = 1)

          dat1_dat2_dat3 <- dat1_dat2_dat3 %>%
            select(c("Protein.Name", "Average.Intensity.dat1", "Average.Intensity.dat2", "Average.Intensity.dat3"))
          # 
          dat1_dat2_dat3$Average.Intensity.dat1 <- log(x = dat1_dat2_dat3$Average.Intensity.dat1, base = 2)
          dat1_dat2_dat3$Average.Intensity.dat2 <- log(x = dat1_dat2_dat3$Average.Intensity.dat2, base = 2)
          dat1_dat2_dat3$Average.Intensity.dat3 <- log(x = dat1_dat2_dat3$Average.Intensity.dat3, base = 2)
          # 
          # 
          row.names(dat1_dat2_dat3) <- dat1_dat2_dat3$Protein.Name
          dat1_dat2_dat3_t <- dat1_dat2_dat3[, -1]
          # 
          dat1_dat2_dat3_max <- as.matrix(dat1_dat2_dat3_t)

          annot_names <- data.frame("Dataset" = c(input$dataset1, input$dataset2, input$dataset3))
          # 
          # # annot_names <- read.delim("E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/my_df_to_use.csv", na.strings = T, sep = ",")
          rownames(annot_names) <- c("Average.Intensity.dat1", "Average.Intensity.dat2", "Average.Intensity.dat3")
          # 
          # 
          col = RColorBrewer::brewer.pal(6, "YlOrRd")
          # 
          # 
          # ann_color = list(
          #   Dataset = c(input$inputreport1$name = "#ff4945", input$inputreport2$name = "#2ac940", input$inputreport3$name = "#75a3e7")
          # )
          # 
          # 
          int_plot <- pheatmap(dat1_dat2_dat3_max, cluster_cols = T, cluster_rows = T,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                   annotation_col = annot_names,
                   cellheight = 12,
                   cellwidth =40,
                   fontsize_row = 10, fontsize_col = 10, border_color = NA,
                   fontsize = 17,
                   treeheight_row =20,
                   show_colnames = F, color = col,
                   # annotation_colors = ann_color[1],
                   annotation_legend = T,
                   legend = T,
                   treeheight_col = 10,
                   annotation_names_col = T
                   ,clustering_method = "average"
                   # scale = "row"
                   # ,legend_labels = c("-3", "-2", "-1", "0", "1", "-2", "-3", "Intensity")
                   # ,filename = "E:/QUT-Pawel-Data/DIA_Extractions_27Nov18/Extractions_Figures2/protein_intensity2.tiff"
                   # , width = 10, height = 10
                   # , main = "Heatmap of protein intensities (log)"
                   , res = 300
          )

          ppp <- grid.arrange(grobs = list(int_plot[[4]]))
          print(ppp)

          incProgress(0.1, detail = "plotting")
          Sys.sleep(0.25)
        }
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
  #     write.table(x = peakviewFormat(data), file = file)
  #     
  #   } else if(input$inputextlibformat == "OpenSwath")
  #   {
  #     write.table(x = OswathFormat(data), file = file)
  #     
  #   } else if(input$inputextlibformat == "Skyline")
  #   {
  #     write.table(x = skylineFormat(data), file = file)
  #     
  #   } else # spectronaut
  #   {
  #     write.table(x = spectronautFormat(data), file = file)
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
  
  #### Reading and Combining Multiple libraries
  
  #### Selecting only one option at a time
  
  observeEvent(input$multireadthree, {
    updateCheckboxInput(session, "multireadfour", value = FALSE)
  })
  
  observeEvent(input$multireadfour, {
    updateCheckboxInput(session, "multireadthree", value = FALSE)
  })
  
  
  ### Combining three libraries
  
  ### Seed library read
  
  tmseedlibdata <- eventReactive (input$tmultiapply, {
    inFile <- input$tmultilib1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib1 <- readLibFile(inFile$datapath, input$tmultilib1format, "spectrum", clean = FALSE)
    return(lib1)
  })
  
  ### External library 1 read
  
  tmextlibdata1 <- eventReactive (input$tmultiapply, {
    inFile <- input$tmultilib2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib2 <- readLibFile(inFile$datapath, input$tmultilib2format, "spectrum", clean = FALSE)
    return(lib2)
  })
  
  ### External library 2 read
  
  tmextlibdata2 <- eventReactive (input$tmultiapply, {
    inFile <- input$tmultilib3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib3 <- readLibFile(inFile$datapath, input$tmultilib3format, "spectrum", clean = FALSE)
    return(lib3)
  })
  
  ### Seed library reading reactive
  
  tmseedlibread <- reactive({
    inFile <- input$tmultilib1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$tmultilib1format, "spectrum", clean = FALSE)
    
  })
  
  ### External Library 1 reading reactive
  
  tmextlibread1 <- reactive({
    inFile <- input$tmultilib2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$tmultilib2format, "spectrum", clean = FALSE)
  })
  
  ### External Library 2 reading reactive
  
  tmextlibread2 <- reactive({
    inFile <- input$tmultilib3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$tmultilib3format, "spectrum", clean = FALSE)
  })
  
  
####################################################################################
  ###############################################################################
  
  ### Combining four libraries
  
  ### Seed library read
  
  fmseedlibdata <- eventReactive (input$fmultiapply, {
    inFile <- input$fmultilib1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib1 <- readLibFile(inFile$datapath, input$fmultilib1format, "spectrum", clean = FALSE)
    return(lib1)
  })
  
  ### External library 1 read
  
  fmextlibdata1 <- eventReactive (input$fmultiapply, {
    inFile <- input$fmultilib2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib2 <- readLibFile(inFile$datapath, input$fmultilib2format, "spectrum", clean = FALSE)
    return(lib2)
  })
  
  ### External library 2 read
  
  fmextlibdata2 <- eventReactive (input$fmultiapply, {
    inFile <- input$fmultilib3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib3 <- readLibFile(inFile$datapath, input$fmultilib3format, "spectrum", clean = FALSE)
    return(lib3)
  })
  
  ### External library 3 read
  
  fmextlibdata3 <- eventReactive (input$fmultiapply, {
    inFile <- input$fmultilib4
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      lib4 <- readLibFile(inFile$datapath, input$fmultilib4format, "spectrum", clean = FALSE)
    return(lib4)
  })
  
  ### Seed library reading reactive
  
  fmseedlibread <- reactive({
    inFile <- input$fmultilib1
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$fmultilib1format, "spectrum", clean = FALSE)
    
  })
  
  ### External Library 1 reading reactive
  
  fmextlibread1 <- reactive({
    inFile <- input$fmultilib2
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$fmultilib2format, "spectrum", clean = FALSE)
  })
  
  ### External Library 2 reading reactive
  
  tmextlibread2 <- reactive({
    inFile <- input$fmultilib3
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$fmultilib3format, "spectrum", clean = FALSE)
  })
  
  ### External Library 3 reading reactive
  
  tmextlibread3 <- reactive({
    inFile <- input$fmultilib4
    
    if(is.null(inFile)) {
      return(NULL)
    } else
      readLibFile(inFile$datapath, input$fmultilib4format, "spectrum", clean = FALSE)
  })
  
  #############################################################################
  ###############################################################################
  
  #### Enable tapply button    
  output$treadapply <- renderUI({
    
    inFile1 <- input$tmultilib1
    inFile2 <- input$tmultilib2
    inFile3 <- input$tmultilib3
    if(is.null(inFile1) && is.null(inFile2) && is.null(inFile3)) return()
    
    # if (is.null(lib1read()) && is.null(lib2read())) return()
    if(!is.null(inFile1) && !is.null(inFile2) && !is.null(inFile3))
    shinyjs::enable("tmultiapply")
  })  
  
  #### Enable fapply button    
  output$freadapply <- renderUI({
    
    inFile1 <- input$fmultilib1
    inFile2 <- input$fmultilib2
    inFile3 <- input$fmultilib3
    inFile4 <- input$fmultilib4
    if(is.null(inFile1) && is.null(inFile2) && is.null(inFile3) && is.null(inFile4)) return()
    
    # if (is.null(lib1read()) && is.null(lib2read())) return()
    if(!is.null(inFile1) && !is.null(inFile2) && !is.null(inFile3) && !is.null(inFile4))
    shinyjs::enable("fmultiapply")
  })  
  
  ######################################################################
  #######################################################################
  
  ### Libraries Summaries (three)
  observeEvent(input$tmultiapply, {
    output$multiseedlibsummary <- renderText({
      seeddata <- tmseedlibdata()
      if(!is.null(seeddata)) {
        summarydata <- libSummary(seeddata)
        
        paste0("Seed assay library contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
    output$multiextlibsummary1 <- renderText({
      extdata1 <- tmextlibdata1()
      if(!is.null(extdata1)) {
        summarydata <- libSummary(extdata1)
        
        paste0("External library 1 contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
    output$multiextlibsummary2 <- renderText({
      extdata2 <- tmextlibdata2()
      if(!is.null(extdata2)) {
        summarydata <- libSummary(extdata2)
        
        paste0("External library 2 contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  ###################################################################
  
  ### Libraries Summaries (four)
  observeEvent(input$fmultiapply, {
    output$multiseedlibsummary <- renderText({
      seeddata <- fmseedlibdata()
      if(!is.null(seeddata)) {
        summarydata <- libSummary(seeddata)
        
        paste0("Seed assay library contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
    output$multiextlibsummary1 <- renderText({
      extdata1 <- fmextlibdata1()
      if(!is.null(extdata1)) {
        summarydata <- libSummary(extdata1)
        
        paste0("External library 1 contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
    output$multiextlibsummary2 <- renderText({
      extdata2 <- fmextlibdata2()
      if(!is.null(extdata2)) {
        summarydata <- libSummary(extdata2)
        
        paste0("External library 2 contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
    output$multiextlibsummary3 <- renderText({
      extdata3 <- fmextlibdata3()
      if(!is.null(extdata3)) {
        summarydata <- libSummary(extdata3)
        
        paste0("External library 3 contains \n",  " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
               "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","),  "\n",
               " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","),  "\n",
               " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))
      }
    })
  })
  
  #################### plot for multi library corr three and corrstats table
  
  multirtcordt3 <- reactive({
    
    if(!is.null(tmseedlibdata()) && !is.null(tmextlibdata1()) && !is.null(tmextlibdata2())) {
      withProgress(message = "Generating correlation stats table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib1 <- tmseedlibdata()
          lib2 <- tmextlibdata1()
          lib3 <- tmextlibdata2()
          
          data <- multiCorrLibThree(lib1, lib2, lib3, label1 = input$tmultilib1$name, label2 = input$tmultilib2$name, label3 = input$tmultilib3$name)
          incProgress(0.1, detail = "computing..")
          Sys.sleep(0.25)
        }
      })
      return(data)
    }
    
   else return(NULL)
    
    
  })
  
  multirtcordt4 <- reactive({
    
    if(!is.null(fmseedlibdata()) && !is.null(fmextlibdata1()) && !is.null(fmextlibdata2()) && !is.null(fmextlibdata3())){
      withProgress(message = "Generating correlation stats table", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib1 <- fmseedlibdata()
          lib2 <- fmextlibdata1()
          lib3 <- fmextlibdata2()
          lib4 <- fmextlibdata3()
          
          data <- multiCorrLibFour(lib1, lib2, lib3, lib4, label1 = input$fmultilib1$name, label2 = input$fmultilib2$name, label3 = input$fmultilib3$name, label4 = input$fmultilib4$name)
          incProgress(0.1, detail = "computing..")
          Sys.sleep(0.25)
        }
      })
      return(data)
    }
    
    else return(NULL)
    
    
  })

  multirtcordt <- eventReactive(input$corcompute, {
    
    if(input$multireadthree == TRUE)
      data <- multirtcordt3()
    else if(input$multireadfour == TRUE)
      data <- multirtcordt4()

    return(data)
  })
  
  observeEvent(input$corcompute, {
    output$librariescorrstats <- DT::renderDataTable({
      
      data <- multirtcordt()
      
      if(is.null(data)) 
        return ("data has no value")
      else{
      
      # data <- multirtcordt4()
      
      DT::datatable(data = data[[1]],
                    options = list(
                      # pageLength = 10,
                                   # scrollX = T,
                                   # scrollY = "500px",
                                   scrollCollapse = T,
                                   autoWidth = T,
                                   scroller.loadingIndicator = T),
                    caption = paste("Correlation with Seed library")
      )
    }
    })
  })
  
  observeEvent(input$corcompute, {
    output$multilibrtcorplot <- renderPlot({
     
      data <- multirtcordt()
      
      if(is.null(data)) 
        return ()
      else {
        # data <- multirtcordt4()
        data <- data[[2]]
        
        print(data)
      }
    })
  })
  
  
  ########################################
  ########################
  
  ## Combining libraries
  
  multicombdt3 <- reactive({
    
    if(!is.null(tmseedlibdata()) && !is.null(tmextlibdata1()) && !is.null(tmextlibdata2())) {
      withProgress(message = "Combining Three Libraries", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib1 <- tmseedlibdata()
          lib2 <- tmextlibdata1()
          lib3 <- tmextlibdata2()
          
          data <- buildSpectralLibThree(lib1, lib2, lib3, method = "time", clean = FALSE, merge = TRUE, nomod = FALSE,
                                        consolidateAccession = input$multiconsacc, plot = input$multigplots, recalibrate = input$multirecalrt,
                                        cutoff.r2 = input$multicutoffr2, cutoff.size = input$multicutofftsize, label1 = input$tmultilib1$name, 
                                        label2 = input$tmultilib2$name, label3 = input$tmultilib3$name)
          incProgress(0.1, detail = "computing..")
          Sys.sleep(0.25)
        }
      })
      return(data)
    }
    
    else return(NULL)
  })
  
  multicombdt4 <- reactive({

    if(!is.null(fmseedlibdata()) && !is.null(fmextlibdata1()) && !is.null(fmextlibdata2()) && !is.null(fmextlibdata3())){
      withProgress(message = "Combining Four Libraries", style = "notification", value = 0.1, {
        for(i in 1:1) {
          lib1 <- fmseedlibdata()
          lib2 <- fmextlibdata1()
          lib3 <- fmextlibdata2()
          lib4 <- fmextlibdata3()

          data <- buildSpectralLibFour(lib1, lib2, lib3, lib4, method = "time", clean = FALSE, merge = TRUE, nomod = FALSE,
                                        consolidateAccession = input$multiconsacc, plot = input$multigplots, recalibrate = input$multirecalrt,
                                        cutoff.r2 = input$multicutoffr2, cutoff.size = input$multicutofftsize, label1 = input$fmultilib1$name, 
                                       label2 = input$fmultilib2$name, label3 = input$fmultilib3$name, label4 = input$fmultilib4$name)
          incProgress(0.1, detail = "computing..")
          Sys.sleep(0.25)
        }
      })
      return(data)
    }

    else return(NULL)
  })

  multicombdt <- eventReactive(input$multicombrun, {

    if(input$multireadthree == TRUE)
      data <- multicombdt3()
    else if(input$multireadfour == TRUE)
      data <- multicombdt4()

    return(data)
  })
  
  
  observeEvent(input$multicombrun, {
    output$multicombineLib <- DT::renderDataTable({
      data <- multicombdt()
      
      # if(is.null(data)) 
      #   return ("data has no value")
      # else{ 
        
        withProgress(message = "Generating combined library table", style = "notification", value = 0.1, {
          for(i in 1:1) {
            
            data <- data[[1]]
            
            incProgress(0.1, detail = "computing..")
            Sys.sleep(0.25)
          }
        })
        
        DT::datatable(data,
                      options = list(
                        pageLength = 10,
                        scrollX = T,
                        scrollY = "500px",
                        scrollCollapse = T,
                        autoWidth = T,
                        scroller.loadingIndicator = T),
                      caption = paste("Combined Library")
        )
      # }
    })
  })
  
  observeEvent(input$multicombrun, {
    output$multicombplot <- renderPlot({

      data <- multicombdt()
      
      if(!is.null(data)) {
        withProgress(message = "Generating combined library plot", style = "notification", value = 0.1, {
          for(i in 1:1) {
            # data = combdatalib()
            data <- data[[1]]

            data2 <- data[!duplicated(data$modification_sequence),]
            x <- aggregate(x = data2$modification_sequence, by = list(Protein = data2$uniprot_id, charge = data2$prec_z),
                           FUN = function(x) {length(x)})
            x$charge <- as.character(x$charge)
            breaks1 <- pretty(range(x$x), n = nclass.FD(x$x), min.n = 1)
            bwidth1 <- breaks1[3] - breaks1[1]

            x2 <- aggregate(x = data$modification_sequence, by = list(Peptide = data$modification_sequence,
                                                                      frg_type = data$frg_type), FUN = function(x2) {length(x2)})
            breaks2 <- pretty(range(x2$x), n = nclass.FD(x2$x), min.n = 1)
            bwidth2 <- breaks2[3] - breaks2[1]

            incProgress(0.1, detail = "plotting")
            Sys.sleep(0.25)
          }
        })

        p1 <- ggplot(data = x, aes(x = x, fill = charge)) +
          geom_histogram(binwidth = bwidth1) +
          labs(x ="Number of Peptides", y = "Number of Proteins") +
          # scale_y_continuous(name = "Frequency") +
          scale_fill_brewer("Precursor Charge", palette = "YlOrRd")+
          theme_classic() +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 11, face = "bold"),
                axis.text.x = element_text(face = "bold", size = 12),
                axis.text.y = element_text(face= "bold", size = 12))

        p2 <- ggplot(data = x2, aes(x = x, fill = frg_type)) +
          geom_histogram(binwidth = bwidth2) +
          labs(x ="Number of Ions", y = "Number of Peptides") +
          scale_fill_brewer("Fragment Type", palette = "PuRd") +
          # scale_fill_manual(values = c("#eeba30", "#ae0001")) +
          # scale_y_continuous(name = "Frequency") +
          # scale_fill_gradient("Frequency", low = "blue", high = "red")+
          theme_classic() +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 11, face = "bold"),
                axis.text.x = element_text(face = "bold", size = 12),
                axis.text.y = element_text(face= "bold", size = 12))

        p <- annotate_figure(ggarrange(p1, p2, ncol = 2), top = "Frequency distribution of Peptides per Protein and Ions per Peptide")
        print(p)

        filepath <- getwd()
        filepath2 <- paste0(filepath,"/graphs")
        ggsave("Combined_Multi_library_plot.png", width = 8, height = 5, path = filepath2)

      } else return()
    })
  })

  observeEvent(input$multicombrun, {
    output$multicomblibsummary <- renderText({
      data <- multicombdt()
    data <- data[[1]]
    summarydata <- libSummary(data)
    paste0("Combined library contains \n", " Proteins = ", formatC(summarydata[["proteins"]], format = "d", big.mark = ","),
           "\n", " Unmodified Peptides = ", formatC(summarydata[["unmodpeptides"]], format = "d", big.mark = ","), "\n",
           " Modified Peptides = ", formatC(summarydata[["modpeptides"]], format = "d", big.mark = ","), "\n",
           " Transitions = ", formatC(summarydata[["transitions"]], format = "d", big.mark = ","))

  })
  })
  
  #///////////////////////////////////////////////////////
  #########################################
  
  ## Downloading Combined Library
  
  ####  Output / Combined Libraries Download 
  output$multidownloadoutputLib <- downloadHandler(filename = function () 
  {paste("Combined", Sys.Date(), ".txt", sep = "_")
  },
  content = function (file) {
    
    data <- multicombdt()
    data <- data[[1]]
    
    if(input$multioutputlibformat == "PeakView"){
      write.table(x = peakviewFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$multioutputlibformat == "OpenSwath")
    {
      write.table(x = OswathFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else if(input$multioutputlibformat == "Skyline")
    {
      write.table(x = skylineFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    } else # spectronaut
    {
      write.table(x = spectronautFormat(data), file = file, row.names = F, na = " ", sep = "\t", quote = FALSE)
    }
    
  }
  )
  
  ##########################################################
  ######################################################
  
  ## Multi Libraries Venn Diagrams
  
  observeEvent(input$multicombrun, {
    output$multivenn1 <- renderPlot({
      data <- multicombdt()
      data <- data[[2]]
      
      if(!is.null(data))
      {
        print(data[[1]])
      }
    })
  })
  observeEvent(input$multicombrun, {
    output$multivenn2 <- renderPlot({
      data <- multicombdt()
      data <- data[[2]]
      
      if(!is.null(data))
      {
        print(data[[2]])
      }
    })
  })
  observeEvent(input$multicombrun, {
    output$multivenn3 <- renderPlot({
      data <- multicombdt()
      data <- data[[2]]
      
      if(!is.null(data))
      {
        print(data[[3]])
      }
    })
  })
  
  ########################################################################################
  ##########################################################
})   ##### Server Function Ends Here
