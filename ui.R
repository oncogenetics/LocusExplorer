# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny


# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    id = "navBarPageID",
    title = "LocusExplorer v0.7",
    windowTitle = "LocusExplorer",
    fluid = FALSE,
    position = "fixed-top",
    inverse = TRUE,
    # Data --------------------------------------------------------------------
    tabPanel(
      "1.Input Data",
      sidebarPanel(
        #push it down 70px to avoid going under navbar
        tags$style(type="text/css", "body {padding-top: 70px;}"),
        
        #hide red error messages
        # tags$style(type="text/css",
        #            ".shiny-output-error { visibility: hidden; }",
        #            ".shiny-output-error:before { visibility: hidden; }"),
        
        #Choose data type
        radioButtons("dataType", h4("Input data:"),
                     c("Prostate iCOGS" = "iCOGS",
                       "Prostate OncoArrayFineMapping" = "OncoArrayFineMapping",
                       "Prostate OncoArrayMeta" = "OncoArrayMeta",
                       "Custom" = "Custom",
                       "Example" = "Example"
                       ),
                     selected = "OncoArrayFineMapping"),
        conditionalPanel("input.dataType == 'iCOGS' ||
                         input.dataType == 'OncoArrayFineMapping' ||
                         input.dataType == 'OncoArrayMeta'",
                         uiOutput("Chr"),
                         uiOutput("RegionID"),
                         selectInput("Flank", label = h5("Region Flank"), 
                                     choices = list( "0 kb" = 1,
                                                     "10 kb" = 10000,
                                                     "50 kb" = 50000,
                                                     "100 kb" = 100000,
                                                     "200 kb" = 200000), 
                                     selected = 10000)
                         ),#conditionalPanel - prostate
        
        conditionalPanel("input.dataType == 'Custom'",
                         fileInput("FileStats", "Association File (required)"),
                         fileInput("FileLD", "LD File (recommended)"),
                         fileInput("FileBedGraph", "bedGraph File"),
                         uiOutput("FileBedGraphName")
        ),#conditionalPanel- Custom
        
        conditionalPanel("input.dataType == 'Example'"
        ),#conditionalPanel- CustomExample
        actionButton("goToPlotSettings", "2.Plot Settings",icon = icon("cogs"),
                     style = "background-color:#85E7FF")
      ),#sidebarPanel
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
                   h4("Summary"),
                   hr(),
                   # if Prostate data is selected then link to relevant paper
                   # Abstract pudmedID
                   conditionalPanel("input.dataType == 'iCOGS' ||
                                     input.dataType == 'OncoArrayFineMapping' ||
                                     input.dataType == 'OncoArrayMeta'",
                                    uiOutput("refProstatePaper"),
                                    hr()
                   ),
                   helpText("UCSC link to selected region:"),
                   htmlOutput("SummaryRegion"),
                   hr(),
                   helpText("NCBI link to hit SNPs:"),
                   dataTableOutput("SummaryHits")),
          #hr() #dataTableOutput("tempSummaryplotDatLD")
          tabPanel("Association",
                   h4("Association"),
                   hr(),
                   dataTableOutput("SummaryStats"),
                   h4("Hit stats"),
                   hr(),
                   dataTableOutput("SummaryHitSNPStats")
                   ),
          tabPanel("LD",
                   h4("Linkage Disequilibrium"),
                   hr(),
                   dataTableOutput("SummaryLD")),
          tabPanel("BedGraph",
                   h4("BedGraph"),
                   hr(),
                   includeMarkdown("Markdown/bedGraphFileFormat.md"),
                   hr(),
                   helpText("bedGraph input data summary"),
                   dataTableOutput("SummaryBedGraph"),
                   hr(),
                   helpText("TCGA eQTL"),
                   dataTableOutput("SummaryROIdatAnnotEQTL")
                   ),
          tabPanel("ENCODE",
                   h4("ENCODE"),
                   hr(),
                   includeMarkdown("Data/wgEncodeBroadHistone/README.md"),
                   helpText("Scores filtered at 5+, and rounded and set maximum value to 100."),
                   dataTableOutput("SummarywgEncodeBroadHistone"),
                   # hr(),
                   # includeMarkdown("Data/wgEncodeRegDnaseClustered/README.md"),
                   # dataTableOutput("SummarywgEncodeRegDnaseClustered"),
                   # hr(),
                   # includeMarkdown("Data/ProstateLNCAP/README.md"),
                   # dataTableOutput("SummaryLNCAP")
                   hr(),
                   includeMarkdown("Data/Annotation/README.md"),
                   dataTableOutput("SummaryAnnotSmooth")
                   
          ),
          tabPanel("Input File Format",
                   h4("Input File Format"),
                   hr(),
                   includeMarkdown("Markdown/InputFileFormat.md"))
        )#tabsetPanel
      )#mainPanel
    ),#tabPanel - Data
    # Plot --------------------------------------------------------------------  
    tabPanel("2.Plot Settings",
             sidebarPanel(
               h4("SNP Filters:"),
               h6("Use sliders to set required threshold for P-value and LD. Filtered SNPs will not be plotted."),
               conditionalPanel("input.plotTypePanelSelected=='Manhattan'",
                                sliderInput("FilterMinPlog",h5("-Log10(PValue)"),
                                            min = 0, max = 5, value = 0.5, step = 0.5),
                                sliderInput("FilterMinLD", h5("LD"),
                                            min = 0, max = 0.9, value = 0, step = 0.05)
                                ),
               
               conditionalPanel("input.plotTypePanelSelected=='Annotation'",
                                selectInput("FilterMinBVS_BF", h5("BVS Bayes Factor:"), 
                                            choices = c(3,5,10,100,1000,10000,30000,50000,100000))
                                ),
               
               conditionalPanel("input.plotTypePanelSelected=='Network'",
                                checkboxGroupInput("haploreg", h5("Haploreg annotation"),
                                             choices = list("SiPhy_cons",
                                                            "EnhancerHistoneMarks",
                                                            "DNAse",
                                                            "ProteinBound",
                                                            "MotifsChanged",
                                                            "GERP_cons",
                                                            "eQTL",
                                                            "Func_INT",
                                                            "Func_NSM",
                                                            "Func_SYN"), 
                                             selected = "SiPhy_cons")
                                ),
               
               
               
               conditionalPanel("input.plotTypePanelSelected=='Manhattan'",
                                #Add suggestiveline and genomewideline
                                numericInput("suggestiveLine",h5("Suggestive Line -Log10(Pvalue)"),
                                             value = 5, min = 0, max = 15, step = 0.1),
                                numericInput("genomewideLine",h5("Genomewide Line -Log10(Pvalue)"),
                                             value = 8, min = 0, max = 15, step = 0.1),
                                #if SNP labels are overlapping then adjust using repulsion force.
                                checkboxInput("adjustLabels","Adjust SNP Lables",TRUE),
                                conditionalPanel("input.adjustLabels",
                                                 sliderInput("repFact",
                                                             h5("Repulsion force factor"),
                                                             min = 0, max = 50, 
                                                             value = 20, step = 0.5)),
                                
                                h4("Zoom region:"),
                                h6("Use sliders to zoom in to required region, or enter region as text, e.g.: chr1:36020000-36140000"),
                                uiOutput("BPrange"),
                                #Zoom using text: chr1:36020000-36140000
                                textInput("RegionZoom", label = h5("Region zoom"),
                                          value = "chr:start-end"),
                                hr(),
                                radioButtons("HitSNPsType", h4("Hit SNPs Type:"),
                                                   c("JAM" = "JAM",
                                                     "Stepwise Forward" = "StepwiseForward"),
                                                   selected = c("StepwiseForward")
                                                   ),
                                hr(),
                                uiOutput("HitSNPs"),
                                h6("Maximum of 5 hit SNPs can be plotted."),
                                hr(),
                                checkboxGroupInput("ShowHideTracks", h4("Tracks:"),
                                                   c("Chromosome"="Chromosome",
                                                     "Manhattan"="Manhattan",
                                                     "Manhattan: Recombination"="Recombination",
                                                     "Manhattan: LDSmooth"="LDSmooth",
                                                     "SNPType"="SNPType",
                                                     "Custom LD"="LD",
                                                     "Custom BedGraph"="BedGraph",
                                                     "wgEncodeBroadHistone"="wgEncodeBroadHistone",
                                                     "wgEncodeRegDnaseClustered"="wgEncodeRegDnaseClustered",
                                                     #"LNCaP Prostate"="LNCAP",
                                                     "Annotation"="annotSmooth",
                                                     "Gene"="Gene"),
                                                   selected=c("Manhattan","Recombination")
                                                   ),
                                
                                h6("Recommneded to hide tracks until final zoom region is decided."),
                                actionButton("resetInput", "Reset inputs",icon = icon("ambulance"),
                                             style = "background-color:#C9DD03"),
                                hr(),
                                actionButton("goToFinalPlot", "3.Final Plot",icon = icon("area-chart"),
                                             style = "background-color:#85E7FF"))
               
               
             ), #sidebarPanel
             
             mainPanel(
               tabsetPanel(id = "plotTypePanelSelected",
                 # tabPanel("Network",
                 #          #OncoArray Graph and Arc plots for SNP relationships
                 #          conditionalPanel("input.dataType == 'OncoArray'",
                 #                           tabsetPanel(type = c("pills"),
                 #                                       tabPanel("Hits vs methods",
                 #                                                forceNetworkOutput("plotGraphStats"),
                 #                                                dataTableOutput("SummaryHaploreg")),
                 #                                       tabPanel("Hits vs other SNPs - R2",
                 #                                                hr(), h4("Hits from all methods"),
                 #                                                forceNetworkOutput("plotGraphLD_All"),
                 #                                                hr(), h4("Stepwise Forward"),
                 #                                                forceNetworkOutput("plotGraphLD_StepwiseForward"),
                 #                                                hr(), h4("ElasticNet"),
                 #                                                forceNetworkOutput("plotGraphLD_ElasticNet"),
                 #                                                hr(), h4("BVS"),
                 #                                                forceNetworkOutput("plotGraphLD_BVS")
                 #                                                
                 #                                                # trying to add borders....FAIL
                 #                                                # splitLayout ...
                 #                                                # verticalLayout(
                 #                                                #   style = "border: 1px solid silver;",
                 #                                                #   cellArgs = list(style = "padding: 6px"),
                 #                                                #   forceNetworkOutput("plotGraphLD_All"),
                 #                                                #   forceNetworkOutput("plotGraphLD_StepwiseForward")#,
                 #                                                #   # forceNetworkOutput("plotGraphLD_StepwiseBackward"),
                 #                                                #   # forceNetworkOutput("plotGraphLD_ElasticNet"),
                 #                                                #   # forceNetworkOutput("plotGraphLD_BVS")  
                 #                                                # )
                 #                                                ),
                 #                                       tabPanel("Hits vs other SNPs vs Methods - R2",
                 #                                                forceNetworkOutput("plotGraphLDStats")),
                 #                                       tabPanel("Hit SNPs - R2",
                 #                                                hr(), h4("All"),
                 #                                                plotOutput("plotArcLD_All",width=800,height=400),
                 #                                                hr(), h4("Forward"),
                 #                                                plotOutput("plotArcLD_StepwiseForward",width=800,height=400),
                 #                                                hr(), h4("ElasticNet"),
                 #                                                plotOutput("plotArcLD_ElasticNet",width=800,height=400),
                 #                                                hr(), h4("BVS"),
                 #                                                plotOutput("plotArcLD_BVS",width=800,height=400)
                 #                                                )
                 #                                       ) #tabsetPanel
                 #                           ) # conditionalPanel - input.dataType == 'OncoArray'
                 #          ), #tabPanel - Network Plots
                 
                 tabPanel("Manhattan",
                          h4("Manhattan"),
                          hr(),
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
                          conditionalPanel("input.ShowHideTracks.indexOf('BedGraph')>-1",
                                           plotOutput("PlotBedGraph",width=800,height=200)),
                          conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeBroadHistone')>-1",
                                           plotOutput("PlotwgEncodeBroadHistone",width=800,height=90)),
                          conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeRegDnaseClustered')>-1",
                                           plotOutput("PlotwgEncodeRegDnaseClustered",width=800,height=60)),
                          # conditionalPanel("input.ShowHideTracks.indexOf('LNCAP')>-1",
                          #                  plotOutput("PlotLNCAP",width=800,height=60)),
                          conditionalPanel("input.ShowHideTracks.indexOf('annotSmooth')>-1",
                                           plotOutput("PlotAnnotSmooth",width=800,height=220)),
                          conditionalPanel("input.ShowHideTracks.indexOf('Gene')>-1",
                                           plotOutput("PlotGene",width=800,height=350))
                          ), # tabPanel - Main Plot
                 # tabPanel("Annotation",
                 #          h4("Annotation"),
                 #          hr()
                 #          ), # tabPanel - Annotation
                 tabPanel("Annotation", 
                          h4("Annotation"),
                          hr(),
                          tabsetPanel(type = "pills",
                                      tabPanel("Annotation 1", 
                                               h4("Annotation 1"),
                                               dataTableOutput("SummarydatAnnotEQTL")),
                                      tabPanel("Annotation 2", 
                                               h4("Sankey (Spaghetti)"),
                                               htmlOutput("PlotSankeyAnnotation")),
                                      tabPanel("Annotation 2 - Data", 
                                               dataTableOutput("SummaryROIdatSankeyAnnotation"))
                                      # tabPanel("BVS Top Models", 
                                      #          htmlOutput("PlotSankeyBVS_TopModels")),
                                      
                                      ) # END tabsetPanel
                          
                          
                          
                 ) # END tabPanel("Annotation",
                 
                 ) # tabsetPanel
               
               
               
               
               
               
               
               #Legend Floating --------------------------------------------------------
               #                absolutePanel(id = "Legend", class = "panel panel-default", fixed = TRUE,
               #                              draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
               #                              width = 200, height = "auto",
               #                              h4("Legend"),
               #                              helpText("Coming soon..."),
               #                              style = "opacity: 0.75")
             ) #mainPanel
    ), #tabPanel - "2.Plot Settings"
    # Final Plot --------------------------------------------------------------  
    tabPanel("3.Final Plot",
             sidebarPanel(
               # Choose Title for merged plot
               uiOutput("downloadPlotTitle"),
               
               
               radioButtons("PlotTheme","Plot theme",
                            choices = list("Shades - Colour picker" = 1,
                                           "Classic - Borders only" = 2,
                                           "None - No colours no borders" = 3),
                            selected = 1),
               
               conditionalPanel("input.PlotTheme == 1",
                                colourInput("PlotThemeColour1",
                                            "Plot theme shade 1", "#C2C2C2"),
                                colourInput("PlotThemeColour2",
                                            "Plot theme shade 2", "#E5E5E5")
               ),
               
               # Choose download filename.
               uiOutput("downloadPlotFileName"),
               helpText("Use PDF or SVG format for further image editing, e.g.: Photoshop, Inkscape."),
               selectInput("downloadPlotType","File type",
                           choices = list(
                             "PDF"  = "pdf",
                             "SVG" = "svg",
                             "JPEG" = "jpeg",
                             "TIFF"  = "tiff"),
                           selected = "jpeg"),
               checkboxInput("downloadPlotAdvancedSettings","Advanced settings",
                             value = FALSE),
               
               # Advanced settings
               conditionalPanel(
                 "(input.downloadPlotAdvancedSettings)",
               wellPanel(style = "background-color: #F6FEAF;",
                         numericInput("downloadPlotWidth",
                                      "Width (pixels)",
                                      value = 1000,
                                      step = 100),
                         
                         numericInput("downloadPlotHeight", 
                                     "Height (pixels)",
                                     value = 1200,
                                     step = 100),
                         
                         sliderInput("downloadPlotPointSize",
                                     "Point size",
                                     value = 12, min = 1, max = 100),
                         
                         conditionalPanel(
                           "(input.downloadPlotType == 'pdf')",
                           selectInput("downloadPlotPaper",
                                       "Paper",
                                       list("a4","letter","legal","us",
                                            "executive","a4r","USr",
                                            "special"),
                                       selected = "special")),
                         
                         conditionalPanel("(input.downloadPlotType == 'jpeg' ||
                                          input.downloadPlotType == 'tiff')",
                                          sliderInput("downloadPlotRes",
                                                      "Resolution (ppi)",
                                                      min = 100, max = 600,
                                                      value = 100,
                                                      step = 10),
                                          #if JPEG input quality
                                          conditionalPanel(
                                            "(input.downloadPlotType == 'jpeg')",
                                            sliderInput("downloadPlotQuality",
                                                        "Quality (%)",
                                                        min = 1, max = 100,
                                                        value = 100,
                                                        step = 5)),
                                          
                                          selectInput("downloadPlotTypeJPEG",
                                                      "Type",
                                                      list("windows","cairo",
                                                           "Xlib","quartz"),
                                                      selected = "cairo"),
                                          
                                          #if TIFF input compression
                                          conditionalPanel(
                                            "(input.downloadPlotType == 'tiff')",
                                            selectInput("downloadPlotCompression",
                                                        "Compression",
                                                        list("none", "rle",
                                                             "lzw", "jpeg",
                                                             "zip", "lzw+p",
                                                             "zip+p"),
                                                        selected = "lzw"))
                                          ),#END JPEG/TIFF 
                         actionButton("resetDownloadPlotSettings",
                                      "Reset inputs",icon = icon("ambulance"),
                                      style = "background-color:#C9DD03;")
                         )#END Advanced settings wellPane
               ),#END Advanced settings 
             
             
             #              pdf(file = ifelse(onefile, "Rplots.pdf", "Rplot%03d.pdf"),
             #                  width, height, onefile, family, title, fonts, version,
             #                  paper, encoding, bg, fg, pointsize, pagecentre, colormodel,
             #                  useDingbats, useKerning, fillOddEven, compress)
             
             #              svg(filename = if(onefile) "Rplots.svg" else "Rplot%03d.svg",
             #                  width = 7, height = 7, pointsize = 12,
             #                  onefile = FALSE, family = "sans", bg = "white",
             #                  antialias = c("default", "none", "gray", "subpixel"))
             
             
             # File downloads when this button is clicked.
             downloadButton("downloadPlot", "Download Plot")
             ), #sidebarPanel plot Download section
             mainPanel(
               #plot merge height is dynamic, based on seleceted tracks
               uiOutput("plotMergeUI")
             )
             
             ), #tabPanel - "Final Plot"
    
    navbarMenu("Help",
               tabPanel("About",
                        h4("About"),
                        hr(),
                        includeMarkdown("README.md")),
               tabPanel("Publication",
                        h4("Publication"),
                        hr(),
                        #PDF
                        img(src="BioinformaticsPaper_p1.PNG"),
                        img(src="BioinformaticsPaper_p2.PNG")),
               tabPanel("Input File Format",
                        h4("Input File Format"),
                        hr(),
                        includeMarkdown("Markdown/InputFileFormat.md")),
               #                tabPanel("Example Plot",
               #                         h4("Example Plot"),
               #                         hr(),
               #                         h4("Sample JPEG output"),
               #                         h5("res = 100, hight = 1200px, width = 1000px"),
               #                         imageOutput("ExamplePlotJPEG")),
               tabPanel("R Session Info",
                        h4("R Session Info"),
                        hr(),
                        includeMarkdown("Markdown/RSessionInfo.md")),
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
               
               tabPanel("FAQ",
                        h4("Frequently Asked Quesitons"),
                        hr(),
                        includeMarkdown("Markdown/FAQ.md"))
    )#navbarMenu - Help
  )#navbarPage
  )#shinyUI
