# Created 13/01/2015 by Tokhir Dadaev
# User interface file for shiny

# Workspace ---------------------------------------------------------------
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "ggplot2", "ggbio", "knitr", "markdown", "stringr"),dependencies = TRUE)
library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(markdown)
library(stringr)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db"))
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
  # Data --------------------------------------------------------------------
  tabPanel(
    "Input Data",
           sidebarPanel(
             tags$style(type="text/css", "body {padding-top: 70px;}"),
             #tags$head(tags$link(rel="shortcut icon", href="www/favicon.ico")),
             radioButtons("dataType", h4("Input Data"),
                          c(#"Prostate"="Prostate",
                            "Custom"="Custom",
                            "Example"="Example"),
                          "Example"),
             #"CustomExample"),
             
             conditionalPanel("input.dataType == 'Prostate'",
                              #Select CHR, 1:23 excluding chromosomes with no hit regions.
                              selectInput("Chr",label=h5("Chr"),
                                          choices=c(1:23)[c(-15,-16,-21)],
                                          selected=1),
                              uiOutput("RegionID")
             ),#conditionalPanel - prostate
             
             conditionalPanel("input.dataType == 'Custom'",
                              fileInput("FileStats", "Association File"),
                              fileInput("FileLD", "LD File"),
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
                        tableOutput("SummaryRegion"),
                        tableOutput("SummaryFileNrowNcol")
                        #tableOutput("SummaryROIdatEQTL"),
                        #textOutput("SummaryRegionFlank")
                        ),
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
  tabPanel("Plot",
           sidebarPanel(
             h4("SNP Filter:"),
             sliderInput("FilterMinPlog",h5("-Log10(PValue)"),
                         min = 0, max = 5, value = 3, step = 0.5),
             sliderInput("FilterMinLD", h5("LD"),
                         min = 0, max = 0.9, value = 0.3,step = 0.05),
             uiOutput("BPrange"),
             uiOutput("RegionZoom"),
             selectInput("Flank", label = h4("Region Flank"), 
                         choices = list("10KB"=10000,
                                        "50KB"=50000,
                                        "1MB"=100000,
                                        "2MB"=200000), 
                         selected = 1),
             uiOutput("HitSNPs"),
             checkboxGroupInput("ShowHideTracks", "Tracks",
                                c("Chromosome"="Chromosome",
                                  "Manhattan"="Manhattan",
                                  "LD"="LD",
                                  "SNPType"="SNPType",
                                  "LNCAP"="LNCAP",
                                  "eQTL"="eQTL",
                                  "Gene"="Gene"),
                                selected=c("Manhattan","LD")),
             h5("Recommneded to hide gene track until final zoom region is decided.")
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
  ), #tabPanel - "Plot"
  # Plot Download -----------------------------------------------------------  
  tabPanel("Plot Download",
           sidebarPanel(
             # Choose Title for merged plot
             uiOutput("downloadPlotTitle"),
             # Choose download filename.
             uiOutput("downloadPlotFileName"),
             selectInput(
               inputId = "downloadPlotType",
               label   = "File type",
               choices = list(
                 "PDF"  = "pdf",
                 "BMP"  = "bmp",
                 "JPEG" = "jpeg",
                 "PNG"  = "png",
                 "SVG" = "svg"),
               selected = "pdf"),
             # Allow the user to set the height and width of the plot download.
             numericInput(
               inputId = "downloadPlotHeight",label = h5("Height (inches)"),
               value = 14,min = 1,max = 100),
             
             numericInput(
               inputId = "downloadPlotWidth",label = h5("Width (inches)"),
               value = 10,min = 1,max = 100),
             p(),
             # File downloads when this button is clicked.
             downloadButton(outputId = "downloadPlot", label = "Download Plot")
             ), #sidebarPanel plot Download section
           mainPanel(
             #plot merge height is dynamic, based on seleceted tracks
             uiOutput("plotMergeUI")
             )
           
           ), #tabPanel - "Plot Download"
  
  tabPanel("Help",
           mainPanel(
             tabsetPanel(
               tabPanel("About",
                        includeMarkdown("Markdown/HelpAbout.md")
               ),
               tabPanel("Example Plot",
                        h4("Coming soon...")),
               tabPanel("R Session Info",
                        includeMarkdown("Markdown/RSessionInfo.md")
               )
             )#tabsetPanel
           )#mainPanel
  )#tabPanel - Help
)#navbarPage
)#shinyUI
