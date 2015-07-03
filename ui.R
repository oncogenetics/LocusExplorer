# Created 13/01/2015 by Tokhir Dadaev
# User interface file for shiny

# Workspace ---------------------------------------------------------------
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "ggplot2", "ggbio", "knitr", "markdown"),dependencies = TRUE)
library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(markdown)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols



# Define UI ---------------------------------------------------------------
shinyUI(navbarPage(
  # Application title
  title = "Locus Explorer v0.2",
  windowTitle = "Locus Explorer",
  
  # Data --------------------------------------------------------------------
  tabPanel("Input Data",
           sidebarPanel(
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
                        tableOutput("SummaryFileNrowNcol")),
               tabPanel("Association",
                        h4("Association"),
                        dataTableOutput("SummaryStats")),
               tabPanel("LD",
                        h4("LD"),
                        dataTableOutput("SummaryLD")),
               tabPanel("LNCAP",
                        h4("LNCAP"),
                        dataTableOutput("SummaryLNCAP")),
               tabPanel("eQTL",
                        h4("eQTL"),
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
             uiOutput("HitSNPs"),
             checkboxGroupInput("ShowHideTracks", "Tracks",
                                c("Chromosome"="Chromosome",
                                  "Manhattan"="Manhattan",
                                  "LD"="LD",
                                  "SNPType"="SNPType",
                                  "LNCAP"="LNCAP",
                                  "eQTL"="eQTL",
                                  "Gene"="Gene"),
                                selected=c("Chromosome","Manhattan","LD","SNPType","LNCAP","eQTL")),
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
             
             #summary Data
             #textOutput("SummaryZoom")
           ) #mainPanel
  ), #tabPanel("Plot"
  # Plot Download -----------------------------------------------------------  
  tabPanel("Plot Download",
           h4("Coming soon...")),
  
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
