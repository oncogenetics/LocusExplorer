# Created 13/01/2015 by Tokhir Dadaev
# User interface file for shiny

# Workspace ---------------------------------------------------------------
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "tidyr", "ggplot2",
#                    "knitr", "markdown", "stringr","DT"),
#                  dependencies = TRUE)
library(shiny)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(markdown)
library(stringr)
library(DT)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "org.Hs.eg.db"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols

# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    title = "Locus Explorer v0.2",
    windowTitle = "Locus Explorer",
    fluid = FALSE,
    position = "fixed-top",
    inverse = TRUE,
    # Data --------------------------------------------------------------------
    tabPanel(
      "Input Data",
      sidebarPanel(
        #push it down 70px to avoid going under navbar
        tags$style(type="text/css", "body {padding-top: 70px;}"),
        #Choose data type
        radioButtons("dataType", h4("Input Data:"),
                     c("Prostate"="Prostate",
                       "Custom"="Custom",
                       "Example"="Example"),
                     selected = "Example"),
        conditionalPanel("input.dataType == 'Prostate'",
                         #Select CHR, 1:23 excluding chromosomes with no hit regions.
                         selectInput("Chr",label=h5("Chr"),
                                     choices=paste0("chr",c(1:14,16:22,"X")),
                                     selected="chr17"),
                         uiOutput("RegionID")
        ),#conditionalPanel - prostate
        
        conditionalPanel("input.dataType == 'Custom'",
                         fileInput("FileStats", "Association File (required)"),
                         fileInput("FileLD", "LD File (required)"),#recommended
                         fileInput("FileLNCAP", "LNCAP File"),
                         fileInput("FileEQTL", "eQTL File")
        ),#conditionalPanel- Custom
        
        conditionalPanel("input.dataType == 'Example'"
        )#conditionalPanel- CustomExample
        
      ),#sidebarPanel
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
                   h4("Summary"),
                   
                   #if Prostate data is selected then link to Ali finemapping paper
                   conditionalPanel("input.dataType == 'Prostate'",
                   hr(),
                   includeMarkdown("Markdown/FinemappingPaperAbstract.md"),
                   hr()
                   ),
                   
                   helpText("From association file the region - UCSC link:"),
                   htmlOutput("SummaryRegion"),
                   hr(),
                   helpText("From LD file hit SNPs are as below - NCBI link:"),
                   dataTableOutput("SummaryHits"),
                   hr(),
                   helpText("Input file number of rows and columns"),
                   dataTableOutput("SummaryFileNrowNcol")),
          tabPanel("Association",
                   h4("Association"),
                   dataTableOutput("SummaryStats")),
          tabPanel("LD",
                   h4("Linkage Disequilibrium"),
                   dataTableOutput("SummaryLD")),
          tabPanel("LNCaP",
                   h4("Prostate Cancer Cells"),
                   dataTableOutput("SummaryLNCAP")),
          tabPanel("eQTL",
                   h4("Expression Quantitative Trait Loci"),
                   dataTableOutput("SummaryEQTL")),
          tabPanel("Input File Format",
                   h4("Input File Format"),
                   includeMarkdown("Markdown/InputFileFormat.md"))
        )#tabsetPanel
      )#mainPanel
    ),#tabPanel - Data
    # Plot --------------------------------------------------------------------  
    tabPanel("Plot Settings",
             sidebarPanel(
               h4("SNP Filters:"),
               h6("Use sliders to set required threshold for P-value and LD. Filtered SNPs will not be plotted."),
               sliderInput("FilterMinPlog",h5("-Log10(PValue)"),
                           min = 0, max = 5, value = 0, step = 0.5),
               sliderInput("FilterMinLD", h5("LD"),
                           min = 0, max = 0.9, value = 0,step = 0.05),
               hr(),
               h4("Zoom region:"),
               h6("Use sliders to zoom in to required region, or enter region as text, e.g.: chr1:36020000-36140000"),
               uiOutput("BPrange"),
               #Zoom using text: chr1:36020000-36140000
               textInput("RegionZoom", label = h5("Region zoom"),
                         value = "chr:start-end"),
               selectInput("Flank", label = h5("Region Flank"), 
                           choices = list("10KB"=10000,
                                          "50KB"=50000,
                                          "1MB"=100000,
                                          "2MB"=200000), 
                           selected = 10000),
               hr(),
               uiOutput("HitSNPs"),
               h6("Maximum of 5 hit SNPs can be plotted."),
               hr(),
               checkboxGroupInput("ShowHideTracks", h4("Tracks:"),
                                  c("Chromosome"="Chromosome",
                                    "Manhattan"="Manhattan",
                                    "LDSmooth"="LDSmooth",
                                    "LD"="LD",
                                    "SNPType"="SNPType",
                                    "LNCAP"="LNCAP",
                                    "eQTL"="eQTL",
                                    "Gene"="Gene"),
                                  selected=c("Manhattan","LD","LDSmooth")),
               h6("Recommneded to hide tracks until final zoom region is decided."),
               #actionButton("resetInput", "Reset inputs",icon = icon("undo"))
               actionButton("resetInput", "Reset inputs",icon = icon("ambulance"),
                            style = "background-color:#C9DD03")
             ), #sidebarPanel
             mainPanel(
               conditionalPanel("input.ShowHideTracks.indexOf('Chromosome')>-1",
                                plotOutput("PlotChromosome",width=800,height=70)),
               conditionalPanel("input.ShowHideTracks.indexOf('Manhattan')>-1",
                                plotOutput("PlotManhattan",width=800,height=500)),
               conditionalPanel("input.ShowHideTracks.indexOf('LD')>-1",
                                plotOutput("PlotSNPLD",width=800,height=110)),
               conditionalPanel("input.ShowHideTracks.indexOf('SNPType')>-1",
                                plotOutput("PlotSNPType",width=800,height=90)),
               conditionalPanel("input.ShowHideTracks.indexOf('LNCAP')>-1",
                                plotOutput("PlotLNCAP",width=800,height=70)),
               conditionalPanel("input.ShowHideTracks.indexOf('eQTL')>-1",
                                plotOutput("PlotEQTL",width=800,height=70)),
               conditionalPanel("input.ShowHideTracks.indexOf('Gene')>-1",
                                plotOutput("PlotGene",width=800,height=350))
             ) #mainPanel
             #  Legend Floating --------------------------------------------------------
             #              absolutePanel(id = "Legend", right = 20, bottom = 20, 
             #                            class = "modal", fixed = FALSE, draggable = TRUE,
             #                            wellPanel(
             #                            h5("Legend")),
             #                            style = "opacity: 0.92")
             
             
    ), #tabPanel - "Plot"
    # Final Plot --------------------------------------------------------------  
    tabPanel("Final Plot",
             sidebarPanel(
               # Choose Title for merged plot
               uiOutput("downloadPlotTitle"),
               # Choose download filename.
               uiOutput("downloadPlotFileName"),
               h6("Use PDF or SVG format for further image editing."),
               selectInput(
                 inputId = "downloadPlotType",
                 label   = h4("File type"),
                 choices = list(
                   "PDF"  = "pdf",
                   "SVG" = "svg",
                   "JPEG" = "jpeg",
                   "TIFF"  = "tiff"),
                 selected = "pdf"),
               # Allow the user to set the height and width of the plot download.
               numericInput(
                 inputId = "downloadPlotHeight",label = h4("Height (inches)"),
                 value = 14,min = 1,max = 100),
               
               numericInput(
                 inputId = "downloadPlotWidth",label = h4("Width (inches)"),
                 value = 10,min = 1,max = 100),
               p(),
               #if jpeg or tiff then set resolution
               conditionalPanel("(input.downloadPlotType == 'jpeg' || 
                                input.downloadPlotType == 'tiff')",
                                sliderInput("downloadPlotResolution",
                                            h4("Resolution"),
                                            min = 100, max = 600, value = 100,
                                            step = 50)),
               # File downloads when this button is clicked.
               downloadButton(outputId = "downloadPlot", label = "Download Plot")
    ), #sidebarPanel plot Download section
    mainPanel(
      #plot merge height is dynamic, based on seleceted tracks
      uiOutput("plotMergeUI")
    )
    
  ), #tabPanel - "Final Plot"
  
  tabPanel("Make LD file",
           sidebarPanel(
             fileInput("FileLDlink", "Unprocessed LDlink output file"),
             
             # File downloads when this button is clicked.
             downloadButton(outputId = "downloadLDFile", label = "Download LD file")
             
           ),
           mainPanel(
             tabsetPanel(
               tabPanel("LDlink file",
                        dataTableOutput("SummaryLDlink")),
               tabPanel("Processed LDlink file",
                        dataTableOutput("SummaryLDlinkProcess")),
               tabPanel("LD Tutorial",
                        includeMarkdown("Markdown/LDTutorial.md"))
             )#tabsetPanel
           )#mainPanel
  ),#tabPanel - Make LD file
  
  
  
  tabPanel("Help",
           mainPanel(
             tabsetPanel(
               tabPanel("About",
                        #includeMarkdown("Markdown/HelpAbout.md")
                        includeMarkdown("README.md")),
               tabPanel("Example Plot",
                        h4("Sample JPEG output"),
                        h5("res = 100, hight = 1200px, width = 1000px"),
                        imageOutput("ExamplePlotJPEG")
                        ),
               tabPanel("R Session Info",
                        includeMarkdown("Markdown/RSessionInfo.md")),
               tabPanel("LD Tutorial",
                        includeMarkdown("Markdown/LDTutorial.md"))
             )#tabsetPanel
           )#mainPanel
  )#tabPanel - Help
  )#navbarPage
)#shinyUI
