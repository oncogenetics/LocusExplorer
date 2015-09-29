# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt


# About -------------------------------------------------------------------
# User interface file for shiny

# Workspace ---------------------------------------------------------------
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "tidyr", "ggplot2",
#                    "knitr", "markdown", "stringr","DT","seqminer",
#                    "lattice","cluster"),
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
library(seqminer)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "org.Hs.eg.db","rtracklayer"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols
library(rtracklayer) # bigwig

#Map fonts to Windows
windowsFonts(Courier=windowsFont("TT Courier New"))
setInternet2(TRUE)

# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    title = "Locus Explorer v0.3",
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
                         fileInput("FileLD", "LD File (recommended)"),
                         fileInput("FileBED", "BED File"),
                         uiOutput("FileBEDName")
                         ),#conditionalPanel- Custom
        
        conditionalPanel("input.dataType == 'Example'"
        )#conditionalPanel- CustomExample
        
      ),#sidebarPanel
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
                   h4("Summary"),
                   hr(),
                   #if Prostate data is selected then link to Ali finemapping paper
                   conditionalPanel("input.dataType == 'Prostate'",
                   hr(),
                   includeMarkdown("Data/ProstateData/README.md"),
                   hr()
                   ),
                   helpText("UCSC link to selected region:"),
                   htmlOutput("SummaryRegion"),
                   hr(),
                   helpText("NCBI link to hit SNPs:"),
                   dataTableOutput("SummaryHits"),
                   hr()
                   #dataTableOutput("tempSummaryplotDatLD")
                   ),
          tabPanel("Association",
                   h4("Association"),
                   hr(),
                   dataTableOutput("SummaryStats")),
          tabPanel("LD",
                   h4("Linkage Disequilibrium"),
                   hr(),
                   dataTableOutput("SummaryLD")),
          tabPanel("BED File",
                   h4("BED File"),
                   hr(),
                   dataTableOutput("SummaryBED")),
          tabPanel("ENCODE",
                   h4("ENCODE"),
                   hr(),
                   includeMarkdown("Data/wgEncodeBroadHistone/README.md"),
                   helpText("Scores filtered at 5+, and rounded and set maximum value to 100."),
                   dataTableOutput("SummarywgEncodeBroadHistone"),
                   hr(),
                   includeMarkdown("Data/wgEncodeRegDnaseClustered/README.md"),
                   dataTableOutput("SummarywgEncodeRegDnaseClustered"),
                   hr(),
                   includeMarkdown("Data/ProstateLNCAP/README.md"),
                   dataTableOutput("SummaryLNCAP")
                   ),
          tabPanel("Input File Format",
                   h4("Input File Format"),
                   hr(),
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
                           min = 0, max = 0.9, value = 0.2, step = 0.05),
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
                                    "Recombination"="Recombination",
                                    "SNPType"="SNPType",
                                    "LDSmooth"="LDSmooth",
                                    "LD"="LD",
                                    "BED"="BED",
                                    "wgEncodeBroadHistone"="wgEncodeBroadHistone",
                                    "wgEncodeRegDnaseClustered"="wgEncodeRegDnaseClustered",
                                    "LNCaP Prostate"="LNCAP",
                                    "Gene"="Gene"),
                                  selected=c("Manhattan","Recombination")
                                  ),
               
               h6("Recommneded to hide tracks until final zoom region is decided."),
               #actionButton("resetInput", "Reset inputs",icon = icon("undo"))
               actionButton("resetInput", "Reset inputs",icon = icon("ambulance"),
                            style = "background-color:#C9DD03")
             ), #sidebarPanel
             mainPanel(
               # Plot title zoomed region, link to UCSC
               htmlOutput("plotTitle", align = "center"),
               # Conditional plots
               conditionalPanel("input.ShowHideTracks.indexOf('Chromosome')>-1",
                                plotOutput("PlotChromosome",width=800,height=70)),
               conditionalPanel("input.ShowHideTracks.indexOf('Manhattan')>-1",
                                plotOutput("PlotManhattan",width=800,height=500)),
               conditionalPanel("input.ShowHideTracks.indexOf('SNPType')>-1",
                                plotOutput("PlotSNPType",width=800,height=70)),
               conditionalPanel("input.ShowHideTracks.indexOf('LD')>-1",
                                plotOutput("PlotSNPLD",width=800,height=110)),
               conditionalPanel("input.ShowHideTracks.indexOf('BED')>-1",
                                plotOutput("PlotBED",width=800,height=60)),
               conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeBroadHistone')>-1",
                                plotOutput("PlotwgEncodeBroadHistone",width=800,height=90)),
               conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeRegDnaseClustered')>-1",
                                plotOutput("PlotwgEncodeRegDnaseClustered",width=800,height=60)),
               conditionalPanel("input.ShowHideTracks.indexOf('LNCAP')>-1",
                                plotOutput("PlotLNCAP",width=800,height=60)),
               conditionalPanel("input.ShowHideTracks.indexOf('Gene')>-1",
                                plotOutput("PlotGene",width=800,height=350)) #,
               
               #Legend Floating --------------------------------------------------------
#                absolutePanel(id = "Legend", class = "panel panel-default", fixed = TRUE,
#                              draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
#                              width = 200, height = "auto",
#                              h4("Legend"),
#                              helpText("Coming soon..."),
#                              style = "opacity: 0.75")
               ) #mainPanel
    ), #tabPanel - "Plot"
    # Final Plot --------------------------------------------------------------  
    tabPanel("Final Plot",
             sidebarPanel(
               # Choose Title for merged plot
               uiOutput("downloadPlotTitle"),
               radioButtons("PlotTheme","Plot theme",
                            choices = list("Gray" = 1,
                                           "Yellow" = 2,
                                           "Green" = 3,
                                           "Classic" = 4,
                                           "None" = 5),
                            selected = 1),
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
                                            min = 100, max = 600, value = 120,
                                            step = 20)),
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
                        h4("LDlink file"),
                        hr(),
                        dataTableOutput("SummaryLDlink")),
               tabPanel("Processed LDlink file",
                        h4("Processed LDlink file"),
                        hr(),
                        dataTableOutput("SummaryLDlinkProcess")),
               tabPanel("LD Tutorial",
                        h4("LD Tutorial"),
                        hr(),
                        includeMarkdown("Markdown/LDTutorial.md"))
             )#tabsetPanel
           )#mainPanel
  ),#tabPanel - Make LD file
  
  
  
  tabPanel("Help",
           mainPanel(
             tabsetPanel(
               tabPanel("About",
                        h4("About"),
                        hr(),
                        includeMarkdown("README.md")),
               tabPanel("Publication",
                        h4("Publication"),
                        hr(),
                        h4("Submitted to Oxford University Press: Bioinformatics"),
                        h4("Status: Awaiting Reviewer Assignment"),
                        img(src="BioinformaticsPaper.PNG")),
               tabPanel("Input File Format",
                        h4("Input File Format"),
                        hr(),
                        includeMarkdown("Markdown/InputFileFormat.md")),
               tabPanel("Example Plot",
                        h4("Example Plot"),
                        hr(),
                        h4("Sample JPEG output"),
                        h5("res = 100, hight = 1200px, width = 1000px"),
                        imageOutput("ExamplePlotJPEG")
                        ),
               tabPanel("R Session Info",
                        h4("R Session Info"),
                        hr(),
                        includeMarkdown("Markdown/RSessionInfo.md")),
               tabPanel("LD Tutorial",
                        h4("LD Tutorial"),
                        hr(),
                        includeMarkdown("Markdown/LDTutorial.md"))
               
             )#tabsetPanel
           )#mainPanel
  )#tabPanel - Help
  )#navbarPage
)#shinyUI
