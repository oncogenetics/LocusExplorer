# Workspace ---------------------------------------------------------------
require(shiny)
#require(data.table)
#require(sqldf)
# require(tcltk)
# require(ggbio)
# require(ggplot2)
# require(GenomicFeatures)
# require(RColorBrewer)
# require(reshape2)
# require(grid) #arrow function for ggplot annotation
# require(shinyIncubator) #progress bar for plots

#data for plotting chromosomes
#data(hg19IdeogramCyto, package = "biovizBase")


#load("plotGenes.RData")

# Define UI ---------------------------------------------------------------
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(title = h3("Locus Explorer - v0.1")),
  
  #Select region
  sidebarPanel(width=3,
               #Select CHR, 1:23 excluding chromosomes with no hit regions.
               selectInput("Chr",label=h5("Chr"),
                           choices=c(1:23)[c(-15,-16,-21)],
                           selected=2),
               uiOutput("RegionID")
               
               
               
  ),#sidebarPanel
  
  mainPanel(
    #  Output Tabs ------------------------------------------------------------
    tabsetPanel(id="conditionedPanels",type="pills",selected="Plot-Static",
                tabPanel("Plot-Static",
                         imageOutput("PlotStatic",width=800,height=800)
                ),
                tabPanel("Plot-Interactive",
                         helpText(h4("Coming soon..."))
                ),
                tabPanel("Data",
                         helpText(h4("Coming soon...")),
                         tags$hr(),
                         helpText(h4("SNP Data")),
                         tags$hr(),
                         helpText(h4("SNP LD Data"))
                ),
                tabPanel("Help",
                         imageOutput("Help",width=1200,height=800)
                ),
                
                tabPanel("About",
                         includeHTML("www/help_about.html")
                )
    )#tabsetPanel
  )#mainPanel
)#pageWithSidebar
)#shinyUI
