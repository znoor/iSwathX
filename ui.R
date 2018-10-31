#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 

#    http://shiny.rstudio.com/
#


packages <- c("shiny", "ggplot2", "dplyr", "stringr", "readr", "DT", "tools",
              "shinydashboard", "utils", "e1071", "shinyjs", "shinythemes", "shinyBS", "pryr",
              "graphics", "ggthemes", "grid", "ggpubr", "plyr", "VennDiagram")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}


library(shiny)
library(ggplot2)
library(ggthemes)
library(stringr)
library(readr)
library(DT)
library(tools)
library(shinydashboard)
library(utils)
library(e1071)
library(shinyjs)
library(shinythemes)
library(shinyBS)

# Define UI for application 
shinyUI(fluidPage(
  
  tags$script("$(document).on('shiny:connected', function(event) {
var myWidth = $(window).width();
              Shiny.onInputChange('shiny_width',myWidth)
              
              });"),

  tags$script("$(document).on('shiny:connected', function(event) {
              var myHeight = $(window).height();
              Shiny.onInputChange('shiny_height',myHeight)
              
              });"),
  # theme = shinytheme("flatly"),
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")),
  useShinyjs(),
  
  # Application title
  #titlePanel(title = "iSwathX Web Application", windowTitle = "iSwathX"),
  
  #Application Dashboad
  fluidRow(
    column(12, 
           dashboardPage(
    dashboardHeader(title = "iSwathX",
                    # title = img(src='title2.png', align = "centre", width = 160, height = 50),
                    disable = FALSE,
                    titleWidth = 300),
    dashboardSidebar(width = 300,
                     # sidebarSearchForm(textId = "search", label = "Search...", icon = shiny::icon("search"), buttonId = "searchbutton"),
                     sidebarMenu(id = "sidemenu",
                                 menuItem(text = "Home",
                                          icon = icon("home", lib = "glyphicon"), 
                                          #badgeLabel = 1, 
                                          #badgeColor = "aqua",
                                          tabName = "introduction",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE,
                                          startExpanded = FALSE),
                                 menuItem(text = "iSwathX Wizard",
                                          icon = icon("road", lib = "glyphicon"), 
                                          #badgeLabel = 1, 
                                          #badgeColor = "aqua",
                                          tabName = "autowizard",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE,
                                          startExpanded = FALSE),
                                 menuItem(text = "Manual Wizard",
                                          icon = icon("industry", "fa-1x"),
                                          selected = FALSE,
                                          startExpanded = TRUE,
                                 menuItem(text = "Libraries Input and Format Conversion",
                                          icon = icon("list", lib = "glyphicon", "fa-1x"), 
                                          #badgeLabel = 2, 
                                          #badgeColor = "aqua",
                                          tabName = "datainput",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE
                                          # startExpanded = FALSE
                                          ),
                                 menuItem(text = "Libraries Comparison",
                                          icon = icon("duplicate", lib = "glyphicon", "fa-1x"), 
                                          #badgeLabel = 3, 
                                          #badgeColor = "aqua",
                                          tabName = "comparelib",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE
                                          # startExpanded = FALSE
                                          ),
                                 menuItem(text = "Building Combined Libraries",
                                          icon = icon("asterisk", lib = "glyphicon", "fa-1x"), 
                                          #badgeLabel = "new", 
                                          #badgeColor = "green",
                                          tabName = "librarycombination",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE
                                          # startExpanded = FALSE
                                          )
                                 ),
                                 menuItem(text = "Statistical Analysis",
                                          icon = icon("stats", lib = "glyphicon", "fa-1x"),
                                          #badgeLabel = 4,
                                          #badgeColor = "green",
                                          tabName = "statanalysis",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE,
                                          startExpanded = FALSE),
                                 menuItem(text = "Data Visualization",
                                          icon = icon("picture", lib = "glyphicon", "fa-1x"), 
                                          #badgeLabel = 4, 
                                          #badgeColor = "green",
                                          tabName = "visualization",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE,
                                          startExpanded = FALSE),
                                 menuItem(text = "Help",
                                          icon = icon("edit", lib = "glyphicon", "fa-1x"), 
                                          #badgeLabel = 4, 
                                          #badgeColor = "green",
                                          tabName = "help",
                                          #href = "links",
                                          newtab = TRUE,
                                          selected = FALSE,
                                          startExpanded = FALSE)),
                     disable = FALSE,
                     collapsed = FALSE),
    dashboardBody(
      tags$style(HTML("

                      
                      .box.box-solid.box-warning>.box-header {
                      color:#fff;
                      background:#eaa932
                      }
                      
                      .box.box-solid.box-warning{
                      border-bottom-color:#eaa932;
                      border-left-color:#eaa932;
                      border-right-color:#eaa932;
                      border-top-color:##eaa932;
                      }
                      
                      ")),
      tabItems(
        tabItem(tabName = "introduction",
                wellPanel(fluidRow(
                  column(6, fluidRow(img(src='title1.png', align = "left", width = "65%", height = "65%")),
                         fluidRow(strong(h4(HTML("A web tool to generate extended peptide reference MS/MS libraries for use in DIA mass spectrometry.")))),
                         fluidRow(strong(HTML("Main functions performed by SwathXtend application include: <br />
                                                  <br />
                                                  1. Library format conversions <br />
                                                  2. Libraries comparisons <br />
                                                  3. Building the combined or extended libraries <br />
                                                  4. Data plotting or Visualization <br />"), style = "background-color: #ffffff;")),
                         fluidRow(column(4, 
                                         strong(h3(tags$a(icon("play"), "iSwathX Wizard", onclick = "openTab('autowizard')"))),
                                  tags$script(HTML("
                                                              var openTab = function(tabName){
                                                              $('a', $('.sidebar')).each(function(){
                                                              if(this.getAttribute('data-value') == tabName) {
                                                              this.click()
                                                              };
                                                              });
                                                              }
                                                              "))),
                                  column(4, 
                                         strong(h3(tags$a(icon("play"), "Manual Mode", onclick = "openTab('datainput')"))),
                                         tags$script(HTML("
                                                          var openTab = function(tabName){
                                                          $('a', $('.sidebar')).each(function(){
                                                          if(this.getAttribute('data-value') == tabName) {
                                                          this.click()
                                                          };
                                                          });
                                                          }
                                                          "))),
                                  column(4, 
                                         strong(h3(tags$a(icon("edit", lib = "glyphicon"), "Help", onclick = "openTab('help')"))),
                                         tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                         )
                                  )
                         ),
                   column(6, img(src='titleimage.png', align = "right", width = "90%", height = "90%"))
                  ),
                  style = "background-color: #ffffff;"),
                splitLayout(cellWidths = c("45%", "55%"),
                            wellPanel(h3("Library Generation Workflow"),
                          wellPanel(fluidRow(column(6,
                                          strong(h4("Step 1")),
                                          img(src='fconvert.png', align = "center", width = "85%", height = "85%"),
                                          h4(tags$a("MSConvert", 
                                                           href = "http://proteowizard.sourceforge.net/tools.shtml", 
                                                    target = "blank",
                                                           align = "center")),
                                          h4(tags$a("qtofpeakpicker", 
                                                           href = "http://proteowizard.sourceforge.net/tools/qtofpeakpicker.html", 
                                                    target = "blank",
                                                           align = "center", span = TRUE))
                                          ),
                                   column(6,
                                          strong(h4("Step 2")),
                                          img(src='dsearch.png', align = "center", width = "80%", height = "80%"),
                                          h4(tags$a("X!Tandem", 
                                                    href = "http://www.thegpm.org/tandem/instructions.html", 
                                                    target = "blank",
                                                    align = "center")),
                                          h4(tags$a("SpectraST", 
                                                    href = "http://tools.proteomecenter.org/wiki/index.php?title=SpectraST", 
                                                    target = "blank",
                                                    align = "center")),
                                          h4(tags$a("Comet", 
                                                    href = "http://comet-ms.sourceforge.net/", 
                                                    target = "blank",
                                                    align = "center", span = TRUE))
                                          ),
                                   hr(width = 1)
                ),
                fluidRow(column(6,
                                strong(h4("Step 3")),
                                img(src='psm.png', align = "center", width = "85%", height = "85%"),
                                h4(tags$a("PeptideProphet", 
                                          href = "http://tools.proteomecenter.org/wiki/index.php?title=Software:PeptideProphet",
                                          target = "blank",
                                          align = "center")),
                                h4(tags$a("iProphet", 
                                          href = "http://tools.proteomecenter.org/wiki/index.php?title=TPP_Tutorial#8._Further_peptide-level_validation_iProphet", 
                                          target = "blank",
                                          align = "center", span = TRUE)),
                                h4(tags$a("ProteinProphet", 
                                          href = "http://tools.proteomecenter.org/wiki/index.php?title=Software:ProteinProphet", 
                                          target = "blank",
                                          align = "center", span = TRUE))
                                ),
                         column(6,
                                strong(h4("Step 4")),
                                img(src='festimation.png', align = "center", width = "80%", height = "80%"),
                                h4(tags$a("Mayu", 
                                          href = "http://proteomics.ethz.ch/muellelu/web/LukasReiter/Mayu/", 
                                          target = "blank",
                                          align = "center"))
                                )
                ),
                fluidRow(column(6,
                                strong(h4("Step 5")),
                                img(src='libgen.png', align = "center", width = "80%", height = "80%"),
                                h4(tags$a("SpectraST", 
                                          href = "http://tools.proteomecenter.org/wiki/index.php?title=SpectraST", 
                                          target = "blank",
                                          align = "center"))
                                ),
                         column(6,
                                strong(h4("Step 6")),
                                img(src='assaygen.png', align = "center", width = "85%", height = "85%"),
                                h4(tags$a("SpectraST", 
                                          href = "http://tools.proteomecenter.org/wiki/index.php?title=SpectraST",
                                          target = "blank",
                                          align = "center")),
                                h4(tags$a("Anaconda", 
                                          href = "https://www.anaconda.com/what-is-anaconda/", 
                                          target = "blank",
                                          align = "center", span = TRUE)),
                                h4(tags$a("OpenMS", 
                                          href = "ftp://ftp.mi.fu-berlin.de/pub/OpenMS/release2.0-documentation/html/index.html", 
                                          target = "blank",
                                          align = "center", span = TRUE))
                                )
                ),style = "background-color: #ffffff;"
                            ),
                style = "background-color: #ffffff;"
                ),
                wellPanel(h3("Libraries Combination Workflow"),
                          style = "background-color: #ffffff;",
                          wellPanel(style = "background-color: #ffffff;",
                                    fluidRow(column(5, h4(strong("iSwathX Wizard"))),
                                             column(3, offset = 4,
                                                    strong(h3(tags$a("Start..", onclick = "openTab('autowizard')"))),
                                             tags$script(HTML("
                                                              var openTab = function(tabName){
                                                              $('a', $('.sidebar')).each(function(){
                                                              if(this.getAttribute('data-value') == tabName) {
                                                              this.click()
                                                              };
                                                              });
                                                              }
                                                              "))
                                             )
                                             
                                      ),
                                    fluidRow(column(12,
                                                    img(src='autoworkflow.png', align = "center", width = "60%", height = "60%"))
                                    )
                                    ),
                          wellPanel(style = "background-color: #ffffff;",
                                    fluidRow(column(5,h4(strong("Manual Wizard"))),
                                      column(3, offset = 4,
                                             strong(h3(tags$a("Start..", onclick = "openTab('datainput')"))),
                                             tags$script(HTML("
                                                              var openTab = function(tabName){
                                                              $('a', $('.sidebar')).each(function(){
                                                              if(this.getAttribute('data-value') == tabName) {
                                                              this.click()
                                                              };
                                                              });
                                                              }
                                                              "))
                                             )
                                    ),
                                    fluidRow(column(12,
                                                    img(src='manualworkflow.png', align = "center", width = "70%", height = "70%"))
                                    )
                          )
                                    
                          )
                )
                
                ),
        tabItem(tabName = "autowizard",
                h3(tags$i("Load standard assay reference libraries")),
                wellPanel(style = "background-color: #ffffff;",
                          fluidRow(column(3, box(selectInput(inputId = "autoformlib1",
                                                      label = "Library Format",
                                                      choices = c("PeakView" = "peakview",
                                                                  "OpenSwath" = "openswath",
                                                                  "Skyline" = "skyline",
                                                                  "Spectronaut" = "spectronaut"),
                                                      selected = NULL,
                                                      multiple = F,
                                                      selectize = T,
                                                      width = '100%'),
                                          fileInput(inputId = "autoseedlib",
                                                    label = "Choose .txt / .csv / .tsv file",
                                                    multiple = TRUE,
                                                    buttonLabel = "Choose file",
                                                    placeholder = "Libraries...",
                                                    accept = c("text/csv/tsv",
                                                               "text/comma-separated-values, text/plain",
                                                               ".csv"),
                                                    width = '100%'),
                                            title = "Seed Libary", solidHeader = TRUE, status = "warning", width = 12
                                          )
                                          ),
                                   column(3, box(selectInput(inputId = "autoformlib2",
                                                      label = "Library Format",
                                                      choices = c("PeakView" = "peakview",
                                                                  "OpenSwath" = "openswath",
                                                                  "Skyline" = "skyline",
                                                                  "Spectronaut" = "spectronaut"),
                                                      selected = NULL,
                                                      multiple = F,
                                                      selectize = T,
                                                      width = '100%'),
                                          fileInput(inputId = "autoextlib",
                                                    label = "Choose .txt / .csv / .tsv file",
                                                    multiple = TRUE,
                                                    buttonLabel = "Choose file",
                                                    placeholder = "Libraries...",
                                                    accept = c("text/csv/tsv",
                                                               "text/comma-separated-values, text/plain",
                                                               ".csv"),
                                                    width = '100%'),
                                          title = "External Libary", solidHeader = TRUE, status = "warning", width = 12
                                   )
                                          ),
                                   column(3,
                                          #uiOutput(outputId = "autoreadapply"),
                                          
                                          fluidRow(tags$head(
                                            tags$style(HTML('#autowizrun{color: #6d7983;border-color:#eaa932}'))
                                          ),
                                          actionButton(inputId = "autowizrun",
                                                       label = "Run",
                                                       icon = icon("refresh", lib = "glyphicon"),
                                                       width = '100%')),
                                          br(),
                                          fluidRow(box(selectInput(inputId = "autooutputlibformat",
                                                                   label = "Library Format",
                                                                   choices = c("PeakView", "OpenSwath", "Skyline", "Spectronaut"),
                                                                   selected = NULL,
                                                                   multiple = F,
                                                                   selectize = T,
                                                                   width = '100%'),
                                                       tags$head(
                                                         tags$style(HTML('#wizdownloadoutputLib{color: #6d7983;border-color:#008000}'))
                                                       ),
                                                       downloadButton(outputId = "wizdownloadoutputLib", label = "Download", width = '100%'),
                                                       title = "Download Combined Library", solidHeader = TRUE, status = "success", width = 16
                                          )
                                          )
                                          ),
                                   column(1),
                                   column(2, 
                                          fluidRow(
                                            strong(h3(tags$a(icon("picture", lib = "glyphicon"), "Data Visualization", onclick = "openTab('visualization')"))),
                                          tags$script(HTML("
                                                           var openTab = function(tabName){
                                                           $('a', $('.sidebar')).each(function(){
                                                           if(this.getAttribute('data-value') == tabName) {
                                                           this.click()
                                                           };
                                                           });
                                                           }
                                                           "))
                                          ),
                                          br(),
                                          fluidRow(
                                            strong(h3(tags$a(icon("edit", lib = "glyphicon"), "Help", onclick = "openTab('help')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                          ),
                                          br(),
                                          fluidRow(
                                            strong(h3(tags$a(icon("home", lib = "glyphicon"), "Home", onclick = "openTab('introduction')"))),
                                            tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                            
                                          )
                                          )
                                   
                            
                          )),
                wellPanel(style = "background-color: #ffffff;",
                          fluidRow(column(7,
                                          wellPanel(style = "background-color: #ffffff;",
                                            h4(tags$i("Combined Library")),
                                          DT::dataTableOutput(outputId = "autocombineLib")
                                          )
                                          ),
                                   column(5,
                                          wellPanel(style = "background-color: #ffffff;",
                                                    h4(tags$i("Parameters Used")),
                                                    verbatimTextOutput("autocombparameters"),
                                                    hr(),
                                                    h4(tags$i("Log")),
                                                    
                                                    verbatimTextOutput("autocomblog")
                                                    )
                                          )
                            
                          )
                          )
                # wellPanel(style = "background-color: #ffffff;",
                #           fluidRow(column(4,
                #           # h4(tags$i("Download Combined Library")),
                #           box(selectInput(inputId = "autooutputlibformat",
                #                       label = "Library Format",
                #                       choices = c("PeakView", "OpenSwath", "Skyline", "Spectronaut"),
                #                       selected = NULL,
                #                       multiple = F,
                #                       selectize = T),
                #               tags$head(
                #                 tags$style(HTML('#wizdownloadoutputLib{color: #6d7983;border-color:#008000}'))
                #               ),
                #           downloadButton(outputId = "wizdownloadoutputLib", label = "Download", width = '300px'),
                #           title = "Download Combined Library", solidHeader = TRUE, status = "success", width = 16
                #           )
                #           )))
          
        ),
                           
                           tabItem(tabName = "datainput",
                                   fluidRow(column(4, h3(tags$i("Reading and Writing Libraries"))),
                                            column(2, offset = 4,
                                                   strong(h3(tags$a(icon("home", lib = "glyphicon"), "Home", onclick = "openTab('introduction')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   ),
                                                   column(2,
                                                          strong(h3(tags$a("Libraries Comparison", onclick = "openTab('comparelib')", icon("angle-double-right")))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   )
                                                   ),
                                   sidebarLayout(
                                     sidebarPanel(
                                       wellPanel(style = "background-color: #ffffff;",
                                         # h4(tags$i("Reading Libraries")),
                                         box(selectInput(inputId = "libformats11",
                                                     label = "Seed Library Format",
                                                     choices = c("PeakView" = "peakview",
                                                                 "OpenSwath" = "openswath",
                                                                 "Skyline" = "skyline",
                                                                 "Spectronaut" = "spectronaut"),
                                                     selected = NULL,
                                                     multiple = F,
                                                     selectize = T,
                                                     width = '100%'),
                                       fileInput(inputId = "seedlib",
                                                 label = "Choose .txt / .csv / .tsv file",
                                                 multiple = TRUE,
                                                 buttonLabel = "Choose file",
                                                 placeholder = "Libraries...",
                                                 accept = c("text/csv/tsv",
                                                            "text/comma-separated-values, text/plain",
                                                            ".csv"),
                                                 width = '100%'),
                                       selectInput(inputId = "libformats12",
                                                   label = "External Library Format",
                                                   choices = c("PeakView" = "peakview",
                                                               "OpenSwath" = "openswath",
                                                               "Skyline" = "skyline",
                                                               "Spectronaut" = "spectronaut"),
                                                   selected = NULL,
                                                   multiple = F,
                                                   selectize = T,
                                                   width = '100%'),
                                       fileInput(inputId = "extlib",
                                                 label = "Choose .txt / .csv / .tsv file",
                                                 multiple = TRUE,
                                                 buttonLabel = "Choose file",
                                                 placeholder = "Libraries...",
                                                 accept = c("text/csv/tsv",
                                                            "text/comma-separated-values, text/plain",
                                                            ".csv"),
                                                 width = '100%'),
                                       uiOutput(outputId = "readapply",
                                                width = '100%'),
                                       tags$head(
                                         tags$style(HTML('#apply{color: #6d7983;border-color:#eaa932}'))
                                       ),
                                       actionButton(inputId = "apply",
                                                    label = "Apply",
                                                    icon = icon("toggle-off")),
                                       title = "Reading Libraries", solidHeader = TRUE, status = "warning", width = 14
                                       )
                                       ),
                                       wellPanel(style = "background-color: #ffffff;",
                                         # h4(tags$i("Cleaning Libraries")),
                                                 box(checkboxInput(inputId = "cleanlib",
                                                               label = "Clean Library",
                                                               value = FALSE,
                                                               width = '100%'),
                                                 conditionalPanel(condition = "input.cleanlib == true",
                                                                  sliderInput(inputId = "intensity",
                                                                              label =  "Intensity Cutoff",
                                                                              min = 0,
                                                                              max = 20,
                                                                              value = 5,
                                                                              step = 2,
                                                                              ticks = TRUE,
                                                                              animate = TRUE,
                                                                              width = '100%'),
                                                                  sliderInput(inputId = "confidence",
                                                                              label = "Confidence Cutoff",
                                                                              min = 0,
                                                                              max = 1,
                                                                              value = 0.99,
                                                                              step = 0.01,
                                                                              ticks = TRUE,
                                                                              animate = TRUE,
                                                                              width = '100%'),
                                                                  checkboxInput(inputId = "modified",
                                                                                label = "Remove modified peptide",
                                                                                value = FALSE,
                                                                                width = '100%'),
                                                                  checkboxInput(inputId = "misscleavage",
                                                                                label = "Remove Missed Cleavages",
                                                                                value = FALSE,
                                                                                width = '100%'),
                                                                  conditionalPanel(condition = "input.misscleavage == true",
                                                                                   selectInput(inputId = "enzyme",
                                                                                               label = "Enzyme",
                                                                                               choices = c("Trypsin"="trypsin",
                                                                                                           "GluC" = "gluc",
                                                                                                           "Chymotrypsin" = "chym"),
                                                                                               selected = "Trypsin",
                                                                                               selectize = TRUE,
                                                                                               width = '100%')
                                                                                   ),
                                                                  tags$head(
                                                                    tags$style(HTML('#cleanupdate{color: #6d7983;border-color:#eaa932}'))
                                                                  ),
                                                                  actionButton(inputId = "cleanupdate",
                                                                               label = "Clean")
                                                                  ),
                                         bsTooltip(id = "cleanlib", title = "Filtering libraries by removing low-intensity and low-confident peptides",
                                                    placement = "right", trigger = "hover"),
                                         title = "Cleaning Libraries", solidHeader = TRUE, status = "warning", width = 14)
                                                ),
                                       wellPanel(style = "background-color: #ffffff;",
                                                 # h4(tags$i("Downloading Libraries")),
                                                 #strong("Select file format"),
                                                 box(selectInput(inputId = "inputseedlibformat",
                                                             label = "Seed Library Format",
                                                             choices = c("PeakView", "OpenSwath", "Skyline", "Spectronaut"),
                                                             selected = NULL,
                                                             multiple = F,
                                                             selectize = T,
                                                             width = '100%'),
                                                     tags$head(
                                                       tags$style(HTML('#downloadseedinputLibs{color: #6d7983;border-color:#008000}'))
                                                     ),
                                                     uiOutput(outputId = "downloadseed"),
                                                 downloadButton(outputId = "downloadseedinputLibs", label = "Download"),
                                                 
                                                 selectInput(inputId = "inputextlibformat",
                                                             label = "External Library Format",
                                                             choices = c("PeakView", "OpenSwath", "Skyline", "Spectronaut"),
                                                             selected = NULL,
                                                             multiple = F,
                                                             selectize = T,
                                                             width = '100%'),
                                                 tags$head(
                                                   tags$style(HTML('#downloadextinputLibs{color: #6d7983;border-color:#008000}'))
                                                 ),
                                                 uiOutput(outputId = "downloadext"),
                                                 downloadButton(outputId = "downloadextinputLibs", label = "Download"),
                                                 title = "Downloading Libraries", solidHeader = TRUE, status = "success", width = 14)
                                       )
                                       # ,
                                       # strong(h3(tags$a("Next", onclick = "openTab('comparelib')"))),
                                       # tags$script(HTML("
                                       #                  var openTab = function(tabName){
                                       #                  $('a', $('.sidebar')).each(function(){
                                       #                  if(this.getAttribute('data-value') == tabName) {
                                       #                  this.click()
                                       #                  };
                                       #                  });
                                       #                  }
                                       #                  "))
                                       ),
                                     mainPanel(
                                       shinyjs::hidden(wellPanel(id = "libreadpanel",
                                                    tabBox(
                                         id = "libraries",
                                         # type = "tabs",
                                         selected = FALSE,
                                         width = '100%',
                                         tabPanel("Seed Library", 
                                                  wellPanel(h5("Library Summary"),
                                                            
                                                            # verbatimTextOutput(outputId = "seedlibsummary", inline = TRUE, container = div),
                                                            verbatimTextOutput(outputId = "seedlibsummary", placeholder = FALSE),
                                                  style = "background-color: #ffffff;"),
                                                  br(),
                                                  DT::dataTableOutput("seedlibcontents"),
                                                  br(),br(),br(),
                                                  plotOutput("lib1plot", width = '100%'
                                                             # click = "seed_click",
                                                             # dblclick = "seed_dblclick",
                                                             # hover = "seed_hover",
                                                             # brush = "seed_brush"
                                                             )
                                                  # verbatimTextOutput("seed_info")
                                                  ),
                                         tabPanel("External Library", 
                                                  wellPanel(h5("Library Summary"),
                                                  
                                                  # htmlOutput(outputId = "extlibsummary", inline = TRUE, container = div),
                                                  verbatimTextOutput(outputId = "extlibsummary", placeholder = FALSE),
                                                  style = "background-color: #ffffff;"),
                                                  br(),
                                                  DT::dataTableOutput("extlibcontents"), 
                                                  br(),br(),br(),
                                                  plotOutput("lib2plot", width = '100%'
                                                             # click = "ext_click",
                                                             # dblclick = "ext_dblclick",
                                                             # hover = "ext_hover",
                                                             # brush = "ext_brush"
                                                             )
                                                  # verbatimTextOutput("ext_info")
                                                  )
                                       )
                                       )
                                       )
                                      ),
                                     fluid = TRUE)),
                           tabItem(tabName = "comparelib",
                                   fluidRow(column(4, h3(tags$i("Analyzing the integration quality of two libraries"))),
                                            column(2, offset = 4,
                                                   strong(h3(tags$a(icon("angle-double-left"),"Reading Libraries", onclick = "openTab('datainput')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   ),
                                            column(2,
                                                   strong(h3(tags$a("Libraries Integration", onclick = "openTab('librarycombination')", icon("angle-double-right")))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   )
                                            ),
                                   navbarPage(
                                              title = "Quality Checks",
                                              tabPanel(title = "Retention Time Correlation", value = "rtcor",
                                                       fluidRow(
                                                         column(6,
                                                                # wellPanel(style = "background-color: #ffffff;",
                                                                          # h4(tags$i("Retention time correlation computations")),
                                                                          box(h5("Select the method for calculating correlation between two libraries"),
                                                                          radioButtons(inputId = "rttime", 
                                                                                       label = "Correlation Method",
                                                                                       choices = c("Time" = "time",
                                                                                                      "Time & Hydrophobicity" = "hydro"),
                                                                                       selected = NULL,
                                                                                       width = '100%'),
                                                                          conditionalPanel(condition = "input.rttime == 'hydro'",
                                                                                           fileInput(inputId = "hydrofile",
                                                                                                     label = "Upload hydrophobicity index file",
                                                                                                     multiple = TRUE,
                                                                                                     buttonLabel = "Choose file",
                                                                                                     placeholder = "hydroindex",
                                                                                                     accept = c("text/csv",
                                                                                                                "text/comma-separated-values, text/plain",
                                                                                                                ".csv"),
                                                                                                     width = '100%')
                                                                                           ),
                                                                          title = "Retention time correlation computations", solidHeader = TRUE, status = "warning", width = 12)
                                                                          
                                                                         
                                                                            # )
                                                                ),
                                                         column(3, fluidRow(
                                                                # wellPanel(style = "background-color: #ffffff;",
                                                                #           h4("Available libraries:"),
                                                                          box(uiOutput(outputId = "libseedfilesrt", inline = TRUE, container = div),
                                                                          
                                                                          title = "Seed Library", solidHeader = TRUE, status = "warning", width = 10, height = 100)
                                                                          
                                                                          # )
                                                         ),
                                                         fluidRow(
                                                           # wellPanel(style = "background-color: #ffffff;",
                                                           #           h4("Available libraries:"),
                                                           box(uiOutput(outputId = "libextfilesrt", inline = TRUE, container = div),
                                                               
                                                               title = "Addon Library", solidHeader = TRUE, status = "warning", width = 10, height = 100)
                                                           
                                                           # )
                                                         )
                                                         ),
                                                         column(3, fluidRow(
                                                           tags$head(
                                                             tags$style(HTML('#calRT{color: #6d7983;border-color:#eaa932}'))
                                                           ),
                                                           actionButton("calRT", "Compute", width = '300px', icon = icon("refresh", lib = "glyphicon"))
                                                         ),
                                                         br(), br(), br(), br(), br(),
                                                                fluidRow(
                                                           # strong(h3(tags$a("Next", onclick = "openTab('librarycombination')"))),
                                                           # tags$script(HTML("
                                                           #                  var openTab = function(tabName){
                                                           #                  $('a', $('.sidebar')).each(function(){
                                                           #                  if(this.getAttribute('data-value') == tabName) {
                                                           #                  this.click()
                                                           #                  };
                                                           #                  });
                                                           #                  }
                                                           #                  "))
                                                         )
                                                                 )
                                                         
                                                       ),
                                                       fluidRow(
                                                         column(6,
                                                                conditionalPanel(condition = "input.rttime == 'time'",
                                                                                 wellPanel(
                                                                          h4(strong("Correlated Data")),
                                                                          DT::dataTableOutput("rtcordata"))
                                                         ),
                                                         
                                                         conditionalPanel(condition = "input.rttime == 'hydro'",
                                                                          wellPanel(
                                                                            h4(strong("Correlated Data")),
                                                                            DT::dataTableOutput("hydrocordata"))
                                                         )
                                                                
                                                                ),
                                                         column(6,
                                                                conditionalPanel(condition = "input.rttime == 'time'",
                                                                                 wellPanel(
                                                                          h4(strong("Correlation Plots")),
                                                                                   plotOutput("rtcorplot", click = "plot_click", width = '100%'),
                                                                          
                                                                                   plotOutput("rtresdplot", click = "plot_click", width = '100%')
                                                                          
                                                                          
                                                                          )
                                                         ),
                                                         
                                                         conditionalPanel(condition = "input.rttime == 'hydro'",
                                                                          wellPanel(
                                                                            h4(strong("Correlation Plots")),
                                                                            plotOutput("hydrocorplot", click = "plot_click", width = '100%'),
                                                                            
                                                                            plotOutput("hydroresdplot", click = "plot_click", width = '100%')
                                                                            
                                                                            
                                                                          )
                                                         )
                                                                )
                                                       )
                                                       ),
                                              tabPanel(title = "Relative Intensity Correlation", value = "ricor",
                                                         fluidRow(
                                                         column(6,
                                                                # wellPanel(style = "background-color: #ffffff;",
                                                                #           h4(tags$i("Relative Intensity correlation computations")),
                                                                          box(h5("Select the method for calculating correlation between two libraries"),
                                                                          radioButtons(inputId = "rintensity",
                                                                                       label = "Correlation Method",
                                                                                       choices = c("Spearman" = "spearman",
                                                                                                   "Kendall" = "kendall"),
                                                                                       selected = NULL,
                                                                                       width = '100%'),
                                                                          title = "Relative Intensity correlation computations", solidHeader = TRUE, status = "warning", width = 12   )
                                                                          # #actionButton("calcrt", "Calculate")

                                                                # )
                                                                ),
                                                         column(3, fluidRow(
                                                                # wellPanel(style = "background-color: #ffffff;",
                                                                #           h4("Available libraries:"),
                                                                          box(uiOutput(outputId = "libseedfilesint", inline = TRUE, container = div),
                                                                          # actionButton("calRI", "Compute"),
                                                                          title = "Seed Library", solidHeader = TRUE, status = "warning", width = 10, height = 100)
                                                   
                                                                # )
                                                                ),
                                                                fluidRow(
                                                                  # wellPanel(style = "background-color: #ffffff;",
                                                                  #           h4("Available libraries:"),
                                                                  box(uiOutput(outputId = "libextfilesint", inline = TRUE, container = div),
                                                                      # actionButton("calRI", "Compute"),
                                                                      title = "Addon Library", solidHeader = TRUE, status = "warning", width = 10, height = 100)
                                                                  
                                                                  # )
                                                                )
                                                                ),
                                                    column(3, fluidRow(
                                                      tags$head(
                                                        tags$style(HTML('#calRI{color: #6d7983;border-color:#eaa932}'))
                                                      ),
                                                      actionButton("calRI", "Compute", width = '300px', icon = icon("refresh", lib = "glyphicon"))
                                                    ),
                                                    br(), br(), br(), br(), br()
                                                    # fluidRow(
                                                    #               strong(h3(tags$a("Next", onclick = "openTab('librarycombination')"))),
                                                    #               tags$script(HTML("
                                                    #   var openTab = function(tabName){
                                                    #   $('a', $('.sidebar')).each(function(){
                                                    #   if(this.getAttribute('data-value') == tabName) {
                                                    #   this.click()
                                                    #   };
                                                    #   });
                                                    #   }
                                                    #   "))
                                                    # )
                                                    )
                                                                
                                                       # #
                                                         ),
                                                        fluidRow(
                                                          column(6,
                                                                 wellPanel(h4(strong("Correlated Data")),
                                                                           DT::dataTableOutput("intensitycordata"))),
                                                          column(6,
                                                                 wellPanel(h4(strong("Correlation Plot")),
                                                                           plotOutput("intensitycorplot", width = '100%')))
                                                        )
                                                       ),
                                              id = "quality")
                                   
                                   ),
                           tabItem(tabName = "librarycombination",
                                   fluidRow(column(4, h3(tags$i("Building an extended reference library by integrating a pair of spectrum libraries"))),
                                            column(2, offset = 4,
                                                   strong(h3(tags$a(icon("angle-double-left"),"Libraries Comparison", onclick = "openTab('comparelib')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   ),
                                            column(2,
                                                   strong(h3(tags$a( "Visualization", onclick = "openTab('visualization')", icon("angle-double-right")))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   )
                                            ),
                                   sidebarLayout(
                                     sidebarPanel(wellPanel(style = "background-color: #ffffff;",
                                                            # h4(tags$i("Parameters")),
                                             box(selectInput(inputId = "combmethod",
                                                         label = "Method",
                                                         choices = c("Time" = "time",
                                                                     "Hydrophobicity" = "hydro",
                                                                     "Hydrophobicity n Sequence" = "hydrosequence"),
                                                         selected = NULL,
                                                         multiple = F,
                                                         selectize = T,
                                                         width = '100%'),
                                             conditionalPanel(condition = "input.combmethod != 'time'",
                                                              fileInput(inputId = "combhydrofile",
                                                                        label = "Upload hydrophobicity index file",
                                                                        multiple = F,
                                                                        buttonLabel = "Choose file",
                                                                        placeholder = "hydroindex",
                                                                        accept = c("text/csv",
                                                                                   "text/comma-separated-values, text/plain",
                                                                                   ".csv")),
                                                              width = '100%'),
                                             h5(strong("Other parameters")),
                                             uiOutput(outputId = "includelength"),
                                             checkboxInput(inputId = "inclength", 
                                                          label = "Include Length",
                                                          value = FALSE,
                                                          width = '100%'),
                                             # checkboxInput(inputId = "parseacc", 
                                             #               label = "Parse Accessions",
                                             #               value = FALSE),
                                             checkboxInput(inputId = "consacc", 
                                                           label = "Consolidate Accessions",
                                                           value = FALSE,
                                                           width = '100%'),
                                             checkboxInput(inputId = "gplots", 
                                                           label = "Generate plots",
                                                           value = TRUE,
                                                           width = '100%'),
                                             checkboxInput(inputId = "mergelib", 
                                                           label = "Merge Libraries",
                                                           value = TRUE,
                                                           width = '100%'),
                                             checkboxInput(inputId = "recalrt", 
                                                           label = "Recalibrate Retention Times",
                                                           value = FALSE,
                                                           width = '100%'),
                                             numericInput(inputId = "cutoffr2",
                                                          label = "Retention Time Correlation Cutoff",
                                                          value = 0.8,
                                                          min = 0,
                                                          max = 1,
                                                          step = 0.1,
                                                          width = '100%'),
                                             numericInput(inputId = "cutofftsize",
                                                          label = "Training Set Size Cutoff",
                                                          value = 50,
                                                          min = 20,
                                                          step = 5,
                                                          width = '100%'),
                                             actionButton(inputId = "combrun",
                                                          label = "Run",
                                                          icon = icon("refresh", lib = "glyphicon"),
                                                          width = '200px'),
                                             tags$head(
                                               tags$style(HTML('#combrun{color: #6d7983;border-color:#eaa932}'))
                                             ),
                                             title = "Parameters", solidHeader = TRUE, status = "warning", width = 14)
                                             
                                             ),
                                   wellPanel(style = "background-color: #ffffff;",
                                     # h4(tags$i("Download Combined Library")),
                                     #strong("Select file format"),
                                    box(selectInput(inputId = "outputlibformat",
                                                 label = "Library Format",
                                                 choices = c("PeakView", "OpenSwath", "Skyline", "Spectronaut"),
                                                 selected = NULL,
                                                 multiple = F,
                                                 selectize = T,
                                                 width = '100%'),
                                      downloadButton(outputId = "downloadoutputLib", label = "Download"),
                                      tags$head(
                                        tags$style(HTML('#downloadoutputLib{color: #6d7983;border-color:#008000}'))
                                      ),
                                     title = "Download Library", solidHeader = TRUE, status = "success", width = 14)
                                     )
                                   # strong(h3(tags$a("Next", onclick = "openTab('visualization')"))),
                                   # tags$script(HTML("
                                   #                  var openTab = function(tabName){
                                   #                  $('a', $('.sidebar')).each(function(){
                                   #                  if(this.getAttribute('data-value') == tabName) {
                                   #                  this.click()
                                   #                  };
                                   #                  });
                                   #                  }
                                   #                  "))
                                   ),
                                   mainPanel(
                                     wellPanel(style = "background-color: #ffffff;",
                                               h4(tags$i("Combined Library")),
                                             DT::dataTableOutput(outputId = "combineLib")))
                                   )
                                   ),
                           # tabItem(tabName = "statanalysis",
                           #         h3(tags$i("Statistical evaluation of extended Libraries")),
                           #         navbarPage(title = "Reliability Checks",
                           #                    tabPanel(title = "Library Checking", value = "libcheck",
                           #                             wellPanel(style = "background-color: #ffffff;",
                           #                                       fluidRow(
                           #                                         column(6,
                           #                                                tags$i(h4("Input Files")),
                           #                                                selectInput(inputId = "inputfiles",
                           #                                                            label = "Select files to analyze",
                           #                                                            choices = c("Use Current Files"="current",
                           #                                                                        "Upload New Files" = "new"),
                           #                                                            selected = "Use Current Files",
                           #                                                            selectize = TRUE,
                           #                                                            width = 300),
                           #                                                conditionalPanel(condition = "input.inputfiles == 'new'",
                           #                                                                 fileInput(inputId = "seedfile",
                           #                                                                           label = "Select seed file (.txt / .csv / .tsv)",
                           #                                                                           multiple = F,
                           #                                                                           buttonLabel = "Choose file",
                           #                                                                           placeholder = "seed",
                           #                                                                           accept = c("text/csv",
                           #                                                                                      "text/comma-separated-values, text/plain",
                           #                                                                                      ".csv")),
                           #                                                                 fileInput(inputId = "extendedfile",
                           #                                                                           label = "Select combined/extended file (.txt / .csv / .tsv)",
                           #                                                                           multiple = F,
                           #                                                                           buttonLabel = "Choose file",
                           #                                                                           placeholder = "extended",
                           #                                                                           accept = c("text/csv",
                           #                                                                                      "text/comma-separated-values, text/plain",
                           #                                                                                      ".csv"))
                           #                                                )),
                           #                                         column(6,
                           #                                                tags$i(h4("Available libraries:")),
                           #                                                uiOutput(outputId = "statlibnames", inline = TRUE, container = div),
                           #                                                actionButton(inputId = "librelcheck",
                           #                                                             label = "Compute",
                           #                                                             icon = icon("toggle-off")
                           #                                                ))
                           #                                       )
                           #                             ),
                           #                             wellPanel(style = "background-color: #ffffff;",width = 600, height = 800,
                           #                                       # uiOutput(outputId = "rellibcheck", inline = TRUE, container = div),
                           #                                       plotOutput("libcheckplot"))),
                           #                    tabPanel(title = "SWATH Results Checking", value = "swathcheck",
                           #                             sidebarPanel(
                           #                               wellPanel(style = "background-color: #ffffff;", 
                           #                                        tags$i(h4("Swath Result Files")),
                           #                                                fileInput(inputId = "seedswathfile",
                           #                                                          label = "Seed library results (.xlsx)",
                           #                                                          multiple = F,
                           #                                                          buttonLabel = "Choose file",
                           #                                                          placeholder = "seed",
                           #                                                          accept = c("text/csv",
                           #                                                                     "text/comma-separated-values, text/plain",
                           #                                                                     ".csv")),
                           #                                                fileInput(inputId = "extendedswathfile",
                           #                                                          label = "Extended library results (.xlsx)",
                           #                                                          multiple = F,
                           #                                                          buttonLabel = "Choose file",
                           #                                                          placeholder = "extended",
                           #                                                          accept = c("text/csv",
                           #                                                                     "text/comma-separated-values, text/plain",
                           #                                                                     ".csv")),
                           #                                                numericInput(inputId = "fdrpass",
                           #                                                             label = "Max. no. of samples passing FDR threshold (<0.01)",
                           #                                                             value = 8,
                           #                                                             min = 1,
                           #                                                             step = 1),
                           #                                                numericInput(inputId = "maxfdrpep",
                           #                                                             label = "Upper threshold value for no. of peptides in a protein ",
                           #                                                             value = 2,
                           #                                                             min = 1,
                           #                                                             step = 1),
                           #                                                actionButton(inputId = "swathr",
                           #                                                             label = "Compute",
                           #                                                             icon = icon("toggle-off"))
                           #                                       )
                           #                             ),
                           #                             mainPanel(
                           #                               shinyjs::hidden(wellPanel(id = "swathtables",
                           #                                                          style = "background-color: #ffffff;",
                           #                                                          tabsetPanel(id = "swathrestables",
                           #                                                                      type = "tabs",
                           #                                                                      selected = FALSE,
                           #                                                                      tabPanel("FDR Bins",
                           #                                                                               DT::dataTableOutput("fdrbinsseed"),
                           #                                                                               DT::dataTableOutput("fdrbinsext")
                           #                                                                               ),
                           #                                                                      tabPanel("Peptide/Proteins Number",
                           #                                                                               DT::dataTableOutput("datcomb")
                           #                                                                      )
                           #                                                            
                           #                                                          )
                           #                                       
                           #                                       )
                           #                               )
                           #                    )
                           #                             )
                           #                    ),
                           #         strong(h3(tags$a("Next", onclick = "openTab('visualization')"))),
                           #         tags$script(HTML("
                           #                          var openTab = function(tabName){
                           #                          $('a', $('.sidebar')).each(function(){
                           #                          if(this.getAttribute('data-value') == tabName) {
                           #                          this.click()
                           #                          };
                           #                          });
                           #                          }
                           #                          "))
                           #         ),
                           tabItem(tabName = "visualization",
                                   fluidRow(column(4, h3(tags$i("Data Visualization"))),
                                            column(2, offset = 2,
                                                   strong(h3(tags$a(icon("home", lib = "glyphicon"),"Home", onclick = "openTab('introduction')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   ),
                                            column(2, strong(h3(tags$a(icon("edit", lib = "glyphicon"), "Help", onclick = "openTab('help')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    "))
                                                   )
                                            ),
                                  fluidRow(tabBox(width = 12,
                                    tabPanel(title = "Correlation Plots",
                                             fluidRow(
                                             column(4, conditionalPanel(condition = "input.rttime == 'time'",
                                                                        wellPanel(style = "background-color: #ffffff;",
                                                                 h4(strong("Retention time correlation"),
                                                                    br(), br(),
                                                                     plotOutput("dvrtcor", width = '100%')))),
                                                    conditionalPanel(condition = "input.rttime == 'hydro'",
                                                                     wellPanel(style = "background-color: #ffffff;",
                                                                               h4(strong("Retention time correlation"),
                                                                                  br(), br(),
                                                                                  plotOutput("dvhydrocor", width = '100%'))))
                                                    ),
                                             column(4, conditionalPanel(condition = "input.rttime == 'time'",
                                                                        wellPanel(style = "background-color: #ffffff;",
                                                                 h4(strong("Retention time correlation residuals"),
                                                                    br(), br(),
                                                                    plotOutput("dvrtresd", width = '100%')))),
                                                    conditionalPanel(condition = "input.rttime == 'hydro'",
                                                                     wellPanel(style = "background-color: #ffffff;",
                                                                               h4(strong("Retention time correlation residuals"),
                                                                                  br(), br(),
                                                                                  plotOutput("dvhydroresd", width = '100%'))))
                                                    ),
                                             column(4, wellPanel(style = "background-color: #ffffff;",
                                                                 h4(strong("Intensity correlation"),
                                                                    br(), br(),
                                                                    plotOutput("dvricor", width = '100%')))
                                                    )
                                    )
                                    # ,
                                    # fluidRow(column(4,
                                    #                 wellPanel(style = "background-color: #ffffff;",
                                    #                    
                                    #                           box(selectInput(inputId = "libcorplots",
                                    #                                    label = "Correlation Plots",
                                    #                                    choices = c("Retention Time Correlation", "Retention Time Correlation Residuals",
                                    #                                                "Relative Intensity Correlation"),
                                    #                                    selected = NULL,
                                    #                                    multiple = F,
                                    #                                    selectize = T,
                                    #                                    width = '100%'),
                                    #                        downloadButton(outputId = "downloadcorplots", label = "Download"),
                                    #                        tags$head(
                                    #                          tags$style(HTML('#downloadcorplots{color: #6d7983;border-color:#008000}'))
                                    #                        ),
                                    #                        title = "Download Correlation Plots", solidHeader = TRUE, status = "success", width = 14)
                                    # )
                                    # )
                                    #   
                                    # )
                      
                                             ),
                                    tabPanel(title = "Libraries Stats Plots",
                                             fluidRow(
                                               column(4, wellPanel(style = "background-color: #ffffff;",
                                                                 h4(strong("Histogram of peptides and ions"),
                                                                    br(), br(),
                                                                    plotOutput("dvhisto", width = '100%')))
                                                    ),
                                             column(4, wellPanel(style = "background-color: #ffffff;",
                                                                 h4(strong("Density plots"),
                                                                    br(), br(),
                                                                    plotOutput("dvdensity", width = '100%')))
                                                    ),
                                             # plotOutput("dvpepvenn"),
                                             # plotOutput("dvprotvenn"),
                                             column(4,wellPanel(style = "background-color: #ffffff;",
                                                                h4(strong("Peptide protein numbers"),
                                                                   br(), br(),
                                                                   plotOutput("dvpepprotnum", width = '100%')))
                                                    )
                                             ),
                                             fluidRow(
                                               column(4, wellPanel(style = "background-color: #ffffff;",
                                                                   h4(strong("Peptide Venn Diagram"),
                                                                      br(), br(),
                                                                      plotOutput("dvpepvenn", width = '100%')))
                                               ),
                                               column(4, wellPanel(style = "background-color: #ffffff;",
                                                                   h4(strong("Protein Venn Diagram"),
                                                                      br(), br(),
                                                                      plotOutput("dvprotvenn", width = '100%')))
                                               ),
                                               column(4,wellPanel(style = "background-color: #ffffff;",
                                                                  h4(strong("Libraries Coverage"),
                                                                     br(), br(),
                                                                     plotOutput("dvlibinfo", width = '100%')))
                                               )
                                             )
                                             # ,
                                             # fluidRow(column(4,
                                             #                 wellPanel(style = "background-color: #ffffff;",
                                             #                           
                                             #                           box(selectInput(inputId = "libanalysisplots",
                                             #                                           label = "Libraries Stats Plots",
                                             #                                           choices = c("Histograms", "Density plots", "Peptide protein numbers",
                                             #                                                       "Peptide venn diagram", "Protein venn diagram", "Libraries coverage"),
                                             #                                           selected = NULL,
                                             #                                           multiple = F,
                                             #                                           selectize = T,
                                             #                                           width = '100%'),
                                             #                               downloadButton(outputId = "downloadlibplots", label = "Download"),
                                             #                               tags$head(
                                             #                                 tags$style(HTML('#downloadlibplots{color: #6d7983;border-color:#008000}'))
                                             #                               ),
                                             #                               title = "Download Library Stats Plots", solidHeader = TRUE, status = "success", width = 14)
                                             #                 )
                                             # )
                                             # 
                                             # )
                                             )
                                    # tabPanel(title = "Reliability Checks Plots",
                                    #          plotOutput("dvlibinfo"))
                                     
                                   )
                                   )
                                   # fluidRow(
                                   # strong(h3(tags$a("Next", onclick = "openTab('help')"))),
                                   # tags$script(HTML("
                                   #                  var openTab = function(tabName){
                                   #                  $('a', $('.sidebar')).each(function(){
                                   #                  if(this.getAttribute('data-value') == tabName) {
                                   #                  this.click()
                                   #                  };
                                   #                  });
                                   #                  }
                                   #                  "))
                                   # )
                                  ),
                           tabItem(tabName = "help",
                                   fluidRow(column(5, h2(strong("iSwathX help and tutorial"))),
                                            
                                            column(3, offset = 2, 
                                                   strong(h3(tags$a(icon("home", lib = "glyphicon"),"Home", onclick = "openTab('introduction')"))),
                                                   tags$script(HTML("
                                                                    var openTab = function(tabName){
                                                                    $('a', $('.sidebar')).each(function(){
                                                                    if(this.getAttribute('data-value') == tabName) {
                                                                    this.click()
                                                                    };
                                                                    });
                                                                    }
                                                                    ")))
                                            
                                            ),
                                   wellPanel(style = "background-color: #ffffff;",
                                                br(),
                                                htmlOutput("helppage")
                                                
                                   )
                                   )
                           ),
                  tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
      # tags$style(
      #   "body{
      #   height: 200px;
      #   max-width: auto;
      #   margin: auto;
      #   }"
      # )
      ),
      
      style = "background-color: #ffffff;"),
   # title = "iSwathX",
    skin = "blue"
  ))),
  fluidRow(
    column(12,
           wellPanel(HTML("This application is built with"), strong(tags$a("RShiny", href = "http://shiny.rstudio.com/")), HTML(" and implements"), strong(tags$i("iSwathX")), HTML("package "),
                      HTML("which is based on previously available package "), 
                     strong(tags$a("SwathXtend", href = "http://bioconductor.org/packages/SwathXtend/")), HTML("("),
                     strong(tags$a("Wu et al., 2016", href = "http://www.mcponline.org/content/15/7/2501.full")), HTML(")."),
                     style = "background-color: #ffffff;")
           )
  )
))

