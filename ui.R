# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny


tweaks <- 
  list(
    tags$head(tags$style(HTML("
                              .multicol {
                              -webkit-column-count: 2; /* Chrome, Safari, Opera */
                              -moz-column-count: 2;    /* Firefox */
                              column-count: 2;
                              -moz-column-fill: auto;
                              -column-fill: auto;
                              }
                              "))),
    #push it down 70px to avoid going under navbar
    tags$style(type = "text/css", "body {padding-top: 70px;}"),
    tags$head(
      tags$style(HTML("
                      #graphplotSNP_LDnetwork {
                        border: 1px solid grey;
                      }
                      "))),
    #hide red error messages
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }")
    )

# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    id = "navBarPageID",
    title = div(h4("LocusExplorer v0.7.2",
                   style = "margin-top: 0px;"),
                img(src = "icr_logo_white_on_black.PNG", height = "70px",
                    style = "position: relative; top: -60px; right: -800px;")),
    windowTitle = "LocusExplorer",
    fluid = FALSE,
    position = "fixed-top",
    inverse = TRUE,
    
    # 1.Input Data ------------------------------------------------------------
    tabPanel(
      # ~~~ sidebarPanel ------------------------------------------------------
      "1.Input Data",
      sidebarPanel(
        tweaks,
        #Choose data type
        radioButtons("dataType", h4("Input data:"),
                     c("Prostate OncoArray Fine-mapping" = "OncoArrayFineMapping",
                       "Custom" = "Custom"
                     ),
                     selected = "OncoArrayFineMapping"),
        
        conditionalPanel("input.dataType == 'OncoArrayFineMapping'", 
                         uiOutput("ui_Chr"),
                         uiOutput("ui_RegionID")
        ),#conditionalPanel - prostate
        
        conditionalPanel("input.dataType == 'Custom'",
                         fileInput("FileStats", "Association File (required)"),
                         fileInput("FileLD", "LD File (recommended)")#,
        ),#conditionalPanel- Custom
        
        conditionalPanel("input.dataType == 'Example'"
        ),#conditionalPanel- CustomExample
        
        actionButton("goToPlotSettings", "2.Plot Settings", icon = icon("cogs"),
                     style = "background-color:#85E7FF")
      ),#sidebarPanel
      # ~~~ mainPanel---------------------------------------------------------
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
                   h4("Summary"),
                   hr(),
                   uiOutput("ui_refProstatePaper"),
                   hr()
          ),
          #hr() #dataTableOutput("tempSummaryplotDatLD")
          tabPanel("Association",
                   h4("Association"),
                   hr(),
                   dataTableOutput("SummaryStats"),
                   h4("Hit stats"),
                   hr(),
                   dataTableOutput("SummaryHitSNPStats")),
          tabPanel("LD",
                   h4("Linkage Disequilibrium"),
                   hr(),
                   dataTableOutput("SummaryLD")),

          tabPanel("Annotation",
                   h4("Annotation"),
                   hr(),
                   helpText("Prostate OncoArray Fine-mapping: annotation"),
                   dataTableOutput("SummaryROIdatAnnot"),
                   hr(),
                   helpText("Prostate OncoArray Fine-mapping: TCGA eQTL"),
                   dataTableOutput("SummaryROIdatAnnotEQTL")
                   #includeMarkdown("Data/Annotation/README.md")
                   ),
          
          tabPanel("ENCODE",
                   h4("ENCODE"),
                   hr(),
                   includeMarkdown("Data/wgEncodeBroadHistone/README.md")
                   #helpText("Scores filtered at 5+, and rounded and set maximum value to 100.")
                   ),
          
          tabPanel("Input File Format",
                   h4("Input File Format"),
                   hr(),
                   includeMarkdown("Markdown/InputFileFormat.md"))
        )#tabsetPanel
        )#mainPanel
    ),#tabPanel - Data
    # 2.Plot Settings ---------------------------------------------------------  
    tabPanel("2.Plot Settings",
             #icon = icon("cogs"),
             # ~~~ sidebarPanel ----------------------------------------------
             sidebarPanel(
               h4("SNP Filters:"),
               h6("Use sliders to set required threshold for P-value and LD."),
               #conditionalPanel("input.plotTypePanelSelected=='Manhattan'",
               splitLayout(sliderInput("FilterMinPlog",h5("-Log10(PValue)"),
                                       min = 0, max = 5, value = 0.5, step = 1),
                           sliderInput("FilterMinLD", h5("LD"),
                                       min = 0, max = 0.9, value = 0, step = 0.1)),

               conditionalPanel("input.plotTypePanelSelected=='Manhattan'",
                                #Add suggestiveline and genomewideline
                                splitLayout(
                                  numericInput("suggestiveLine", 
                                               h5("Suggestive"),
                                               value = 5, min = 0, max = 15, step = 0.1),
                                  numericInput("genomewideLine",
                                               h5("Genomewide"),
                                               value = 8, min = 0, max = 15, step = 0.1)),
                                sliderInput("Flank", label = h5("Region Flank"), 
                                            min = 0,
                                            max = 100,
                                            #step = 10,
                                            round = TRUE,
                                            value = 10, post = "kb"),
                                h4("Zoom region:"),
                                uiOutput("ui_BPrange")
               ),
               hr(),
               uiOutput("ui_HitSNPs"),
               conditionalPanel("input.plotTypePanelSelected=='Manhattan'",
                                hr(),
                                uiOutput("ui_otherHits"),
                                hr(),
                                splitLayout(
                                  checkboxGroupInput("ShowHideManhattanPvalues", "Manhattan: Pvalues",
                                                     c("Manhattan" = "Manhattan",
                                                       "Recombination" = "Recombination",
                                                       "LD smooth" = "LDSmooth",
                                                       "Effect" = "Effect"),
                                                     selected = c("Manhattan", "Recombination", "LDSmooth")),
                                  conditionalPanel("input.dataType == 'OncoArrayFineMapping'",
                                                   checkboxGroupInput("ShowHideManhattanPostProbs", "Manhattan: PostProbs",
                                                                      c("Manhattan" = "Manhattan",
                                                                        "Recombination" = "Recombination",
                                                                        "LD smooth" = "LDSmooth"),
                                                                      selected = "Manhattan")
                                  )
                                                   ),
                                checkboxGroupInput("ShowHideTracks", "Tracks:",
                                                   c("Chromosome" = "Chromosome",
                                                     "SNP type" = "SNPType",
                                                     "Hit SNPs LD" = "LD",
                                                     "ENCODE H3k27ac" = "wgEncodeBroadHistone",
                                                     "ENCODE H3k4me1" = "wgEncodeBroadHistone_H3k4me1",
                                                     "ENCODE H3k4me3" = "wgEncodeBroadHistone_H3k4me3",
                                                     "Annotation" = "annotOncoFinemap",
                                                     "Gene" = "Gene",
                                                     "Caption" = "Caption"),
                                                   selected = c("Manhattan", "Recombination")),
                                # Caption settings --------------------------------------------------------
                                #checkboxInput("IncludeCaption", label = "Include caption", value = TRUE),
                                uiOutput("ui_captionText"),
                                hr(),
                                splitLayout(
                                  actionButton("resetInput", "Reset inputs",
                                               icon = icon("ambulance"),
                                               style = "background-color:#C9DD03"),
                                  
                                  actionButton("goToFinalPlot", "3.Final Plot",
                                               icon = icon("area-chart"),
                                               style = "background-color:#85E7FF"))
               ), # END  conditionalPanel("input.plotTypePanelSelected=='Manhattan'"
               
               # Network settings --------------------------------------------------------
               conditionalPanel("input.plotTypePanelSelected=='LD-Network'",
                                hr(),
                                checkboxInput("hitsOnly", label = "Show hit SNPs only", value = TRUE),
                                selectInput("networkLayout", label = "Network layout",
                                            choices = list("nicely" = "layout_nicely",
                                                           "circle" = "layout_in_circle",
                                                           "star" = "layout_as_star",
                                                           "components" = "layout_components",
                                                           "grid" = "layout_on_grid",
                                                           "sphere" = "layout_on_sphere",
                                                           "random" = "layout_randomly",
                                                           "dh" = "layout_with_dh",
                                                           "drl" = "layout_with_drl",
                                                           "fr" = "layout_with_fr",
                                                           "gem" = "layout_with_gem",
                                                           "graphopt" = "layout_with_graphopt",
                                                           "kk" = "layout_with_kk",
                                                           "lgl" = "layout_with_lgl",
                                                           "mds" = "layout_with_mds",
                                                           "sugiyama" = "layout_with_sugiyama"),
                                            selected = "layout_nicely"))
             ), # END sidebarPanel
             # ~~~ mainPanel -----------------------------------------------
             mainPanel(
               tabsetPanel(id = "plotTypePanelSelected",
                           # Manhattan -------------------------------------
                           tabPanel("Manhattan",
                                    h4("Manhattan"),
                                    hr(),
                                    # Plot title zoomed region, link to UCSC
                                    htmlOutput("ui_plotTitle", align = "center"),
                                    # Conditional plots
                                    conditionalPanel("input.ShowHideTracks.indexOf('Chromosome')>-1",
                                                     plotOutput("PlotChromosome",width=800,height=70)),
                                    
                                    # Manhattan with p-values
                                    conditionalPanel("input.ShowHideManhattanPvalues.indexOf('Manhattan')>-1",
                                                     plotOutput("PlotManhattanPvalues",width=800,height=500)),
                                    # Manhattan with PostProbs
                                    conditionalPanel("input.ShowHideManhattanPostProbs.indexOf('Manhattan')>-1",
                                                     plotOutput("PlotManhattanPostProbs",width=800,height=250)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('SNPType')>-1",
                                                     plotOutput("PlotSNPType",width=800,height=70)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('LD')>-1",
                                                     plotOutput("PlotSNPLD",width=800,height=110)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeBroadHistone')>-1",
                                                     plotOutput("PlotwgEncodeBroadHistone",width=800,height=90)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeBroadHistone_H3k4me1')>-1",
                                                     plotOutput("PlotwgEncodeBroadHistone_H3k4me1",width=800,height=90)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('wgEncodeBroadHistone_H3k4me3')>-1",
                                                     plotOutput("PlotwgEncodeBroadHistone_H3k4me3",width=800,height=90)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('annotOncoFinemap')>-1",
                                                     plotOutput("PlotAnnotOncoFinemap",width=800,height=120)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('Gene')>-1",
                                                     plotOutput("PlotGene",width=800,height=350)),
                                    conditionalPanel("input.ShowHideTracks.indexOf('Caption')>-1",
                                                     plotOutput("PlotCaption",width=800,height=100))
                                    
                           ), # END tabPanel("Manhattan"
                           # LD-Heatmap ---------------------------------------
                           tabPanel("LD-Heatmap",
                                    h4("LD-Heatmap"),
                                    hr(),
                                    plotOutput("plotSNP_LDheatmap", width = 800, height = 650)
                                    ),
                           # LD-Network ---------------------------------------
                           tabPanel("LD-Network",
                                    h4("LD-Network"),
                                    h5("Note: this work is still experimental, and will be improved more in future versions."),
                                    hr(),
                                    visNetworkOutput("plotSNP_LDnetwork", width = "100%", height = "800px")
                                    )
                           
               ) # END tabsetPanel(id = "plotTypePanelSelected",
               
             ) # END mainPanel
    ), #tabPanel - "2.Plot Settings"
    # 3.Final Plot ------------------------------------------------------------  
    tabPanel("3.Final Plot",
             # ~~~ sidebarPanel ----------------------------------------------
             sidebarPanel(
               # Choose Title for merged plot
               uiOutput("ui_downloadPlotTitle"),
               
               radioButtons("PlotTheme","Plot theme",
                            choices = list("Shades - Colour picker" = 1,
                                           "Classic - Borders only" = 2,
                                           "None - No colours no borders" = 3),
                            selected = 1),
               
               conditionalPanel("input.PlotTheme == 1",
                                splitLayout(
                                  cellArgs = list(style = "overflow: visible;"),
                                  colourInput("PlotThemeColour1",
                                              "Plot theme shade 1", "#C2C2C2"),
                                  uiOutput("ui_PlotThemeColour2")
                                  # colourInput("PlotThemeColour2",
                                  #             "Plot theme shade 2", "#E5E5E5")
                                )
               ),
               
               # Choose download filename.
               uiOutput("ui_downloadPlotFileName"),
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
                              "Reset inputs", icon = icon("ambulance"),
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
             # ~~~ mainPanel --------------------------------------------------
             mainPanel(
               #plot merge height is dynamic, based on seleceted tracks
               uiOutput("ui_plotMergeUI")
             )
             
    ), #tabPanel - "Final Plot"
    # 4.Help ------------------------------------------------------------------
    navbarMenu(title = "Help",
               #title = NULL,
               #icon = icon("ambulance"),
               icon = icon("info"),
               tabPanel("About",
                        h4("About"),
                        hr(),
                        includeMarkdown("README.md")),
               # tabPanel("Publication",
               #          h4("Publication"),
               #          hr(),
               #          #PDF
               #          img(src = "BioinformaticsPaper_p1.PNG"),
               #          img(src = "BioinformaticsPaper_p2.PNG")),
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
                          downloadButton(outputId = "downloadLDFile",
                                         label = "Download LD file")
                          
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
               tabPanel("Prostate raw data",
                        h4("Prostate raw data"),
                        hr(),
                        splitLayout(
                        img(src = "logoICR.png"),
                        img(src = "PRACTICAL_logo_xs.png")),
                        hr(),
                        includeMarkdown("Markdown/PracticalData.md")
                        ),
               tabPanel("FAQ",
                        h4("Frequently Asked Quesitons"),
                        hr(),
                        includeMarkdown("Markdown/FAQ.md"))
    )#navbarMenu - Help
    )#navbarPage
  )#shinyUI



# TESTING -----------------------------------------------------------------

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



#Legend Floating --------------------------------------------------------
#                absolutePanel(id = "Legend", class = "panel panel-default", fixed = TRUE,
#                              draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
#                              width = 200, height = "auto",
#                              h4("Legend"),
#                              helpText("Coming soon..."),
#                              style = "opacity: 0.75")



#,

# Debug: testing objects ------------------
#h4("testing vars"),
#DT::dataTableOutput("testVars")


#        HTML(
#          ".checkbox-inline { 
#             margin-left: 0px;
#             margin-right: 10px;
#   }
#  .checkbox-inline+.checkbox-inline {
#             margin-left: 0px;
#             margin-right: 10px;
#   }
# "
#        )
