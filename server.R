# Workspace ---------------------------------------------------------------
#require(dplyr)
require(shiny)
require(data.table)
# require(sqldf)
# require(tcltk)
# require(ggplot2)
# require(ggbio)
# require(GenomicFeatures)
# require(RColorBrewer)
# require(reshape2)
# require(grid) #arrow function for ggplot annotation
# require(shinyIncubator) #progress bar for plots

#data for plotting chromosomes
#data(hg19IdeogramCyto, package = "biovizBase")


#load("plotGenes.RData")


# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {
  # Connect database ------------------------------------------------------
  #con <- dbConnect(SQLite(), dbname="meta.sqlite")
  
  # Read reference files ----------------------------------------------------
#   regions <- data.table(
#     dbGetQuery(con,paste("select chr,RegionID, 
#                          min(position) Xstart,
#                          max(position) Xend,
#                          max(P_value_Log) P_value_Log_Max
#                          from meta 
#                          group by RegionID 
#                          order by RegionID")))
  regions <- fread("regions.csv")
#   regionStart <- reactive({regions[ RegionID==input$RegionID, Xstart]})
#   regionEnd <- reactive({regions[ RegionID==input$RegionID, Xend]})
#   regionP_value_Log_Max <- reactive({
#     maxY <- regions[ RegionID==input$RegionID, P_value_Log_Max]
#     return(max(10,ceiling((maxY+1)/5)*5))
#   })
  
  #Plot caption
  #captionText <- read.csv("regions.csv",as.is=TRUE)
  
  # Read plot dat -----------------------------------------------------------
#   dat_region_SQL <- reactive({
#     data.table(dbGetQuery(con,
#                           paste0("select * 
#                                     from meta 
#                                     where RegionID='",input$RegionID,"'")))})
#   
#   dat_region <- reactive({
#     dat_region_SQL()[ P_value_Log >= input$Yrange, ]
#   })
#     
  # Data prep eQTL CIS tumour ------------------------------------------------
#   dat_eQTL <- data.table(dbGetQuery(con,"select * from eQTL"))
#   
#   cis_tumour <-
#     dat_eQTL %>%
#     mutate(Direction=ifelse(t.stat>0,"Up",
#                             ifelse(t.stat<0,"Down",NA)))
#   #normalise t.stat to -1-0-1 range
#   cis_tumour_norm <- 
#     rbind(
#       #scale positive t.stats to 0-1
#       cis_tumour %>% 
#         filter(Direction=="Up") %>%
#         mutate(TStat=t.stat/max(t.stat)),
#       #scale negative t.stats to -1-0
#       cis_tumour %>% 
#         filter(Direction=="Down") %>%
#         mutate(
#           temp=abs(t.stat),
#           TStat=-(temp/max(temp))) %>%
#         select(-temp)) %>%
#     #sort by chr and position
#     arrange(chr,pos)
#   #subset CIS for selected region
#   cis_region <- reactive({
#     cis_tumour_norm %>% 
#       filter(chr==input$Chr &
#                pos>=regionStart() &
#                pos<=regionEnd())})
#   
  # ROI ---------------------------------------------------------------------
  #set Region of Interest - to subset gg_merge
#   roi <- reactive({
#     GRanges(seqnames=ifelse(input$Chr==23,"chrX",
#                             paste0("chr",input$Chr)),
#             IRanges(start=input$Xrange[1],
#                     end=input$Xrange[2],
#                     names=input$RegionID))})
#   
  #region SNPs
#   SNP_type <- read.csv("SNP_Type.csv",as.is=TRUE)
#   regionSNPs <- reactive({
#     hitSNPs <- dat_region()[ HitSNP==1,snp_name]
#     #convert to 1KG names-1KG SNP match
#     hitSNPs <- SNP_type[match(hitSNPs,SNP_type$SNP),"Match1KG"]
#     return(hitSNPs)})
  
  
  # LD data -----------------------------------------------------------------
  #get 1KG LD data
#   LD_1KG <- reactive({
#     #Choose colours for LD
#     nsnp <- length(regionSNPs())
#     mycols <- brewer.pal(ifelse(nsnp<3,3,nsnp),"Set1")
#     #if number of HitSNPs is >9 then recycle colours
#     mycols <- if(nsnp>9){cbind(mycols,1:nsnp)[,1]}else{mycols}
#     #create pallete LD 0 to 100
#     colLD <- lapply(mycols,function(i)colorRampPalette(c("grey95",i))(100))
#     #assign colours to LD 
#     d_LD_1KG<-base::do.call(rbind,
#                             lapply(regionSNPs(),function(snp){
#                               d <- dbGetQuery(con,paste0("select * from LD_1KG where SNP_A='",snp,
#                                                          "' and R2 >= ",
#                                                          input$LDrange))
#                               LDColIndex <- ifelse(round(d$R2,2)==0,1,round(d$R2,2)*100)
#                               LDColIndex <- ifelse(LDColIndex>100,100,LDColIndex)
#                               d$LDSNP <- snp
#                               d$LDSmoothCol <- colLD[[match(snp,regionSNPs())]][100]
#                               d$LDCol <- colLD[[match(snp,regionSNPs())]][LDColIndex]
#                               return(d)
#                             }))
#     # add pvalues for Y value on the plot
#     d_LD_1KG <- 
#       base::merge(d_LD_1KG,dat_region(),by.x="BP_B",by.y="position",all.x=TRUE)
#     d_LD_1KG$R2_Adj <- regionP_value_Log_Max()*d_LD_1KG$R2
#     return(data.table(d_LD_1KG))
#   })
#   
  
  # Recomb Data -------------------------------------------------------------
#   regionRecomb <- reactive({
#     d_regionRecomb <- 
#       dbGetQuery(con,paste0("select *
#                             from recomb_1KG
#                             where chr=",input$Chr," and 
#                             position between ",
#                             regions[RegionID==input$RegionID,Xstart],
#                             " and ",
#                             regions[RegionID==input$RegionID,Xend]))
#     #adjust R2 to max Y value = P_value_Log_Max
#     d_regionRecomb$recomb_adj <- 
#       d_regionRecomb$recomb*regionP_value_Log_Max()/100
#     #return data for ggplot
#     return(d_regionRecomb)})
#   # LNCAP SNPs - colour arrows ----------------------------------------------
#   #subset LNCAP
#   LNCAP <- reactive({dbGetQuery(con,paste("select * from LNCAP where chr=",
#                                           input$Chr,
#                                           "and position between",
#                                           input$Xrange[1],
#                                           "and",
#                                           input$Xrange[2]))})
#   
  # LNCAP SNPs - colour arrows ----------------------------------------------
  #subset LNCAP
  #   LNCAP_dense <- reactive({
  #     temp_bw <- (regionEnd()-regionStart())/100
  #     d1 <- density(LNCAP()$position, bw=temp_bw, n=temp_bw)
  #     d1 <- data.table(POS=d1$x,LNCAP=d1$y)
  #     d2 <- data.table(POS=c(input$Xrange[1],input$Xrange[2]))
  #     d1[d2][between(pos, start, end)]
  #     temp[temp_region,nomatch=0]
  #     temp <- temp[ temp$POS >= input$Xrange[1] &
  #                     temp$POS <= input$Xrange[2],]
  #     return(temp)
  #     })
  #   
  
  
  # Dynamic UI --------------------------------------------------------------
  #RegionID list based on input$Chr
  output$RegionID <- renderUI({
    selectInput(inputId="RegionID", label= h5("RegionID"), 
                choices=regions[ Chr==input$Chr,RegionID],
                selected="2_4")})
  
  #zoom to region X axis BP
#   output$Xrange <-
#     renderUI({
#       sliderInput("Xrange", h5("Region: Start-End"),
#                   min = regions[ RegionID==input$RegionID,Xstart]-300000,
#                   max = regions[ RegionID==input$RegionID,Xend]+300000,
#                   value = c(regions[ RegionID==input$RegionID,Xstart],
#                             regions[ RegionID==input$RegionID,Xend]),
#                   step = 50000)})
#   #download file default name is RegionID
#   output$downloadPlotFileName <- renderUI({
#     textInput(
#       inputId = "downloadPlotFileName",
#       label = h5("Download file name"),
#       value = input$RegionID)})
#   
  # Output ------------------------------------------------------------------  
  
  # Tab - Hyperlinks --------------------------------------------------------
  # Hyperlinks to UCSC region and NCBI SNPs
#   output$Hyperlink <- renderTable(function(){
#     data.frame(
#       Hyperlink=
#         c(
#           #ROI link to UCSC
#           paste0("<a href='http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr",
#                  ifelse(input$Chr==23,"X",input$Chr),":",
#                  round(input$Xrange[1]),"-",round(input$Xrange[2]),
#                  "' target='_blank'>",
#                  #Chr??:start-end
#                  "chr",input$Chr,":",input$Xrange[1],"-",input$Xrange[2],
#                  "</a>"),
#           
#           #Hit SNPs link to NCBI
#           sapply(regionSNPs(),
#                  function(snp)paste0("<a href='http://www.ncbi.nlm.nih.gov/snp/?term=",
#                                      snp,"' target='_blank'>",
#                                      snp,"</a>"))),
#       row.names = c("UCSC",paste0("SNP",1:length(regionSNPs())))
#     )
#   },
#   sanitize.text.function = function(x) x)
#   
#   #caption for plots
#   captionPlot <- reactive({captionText[ captionText$RegionID == input$RegionID,"Caption"]})
#   output$PlotCaptionStatic <- renderText({captionPlot()})
#   output$PlotCaptionInteractive <- renderText({captionPlot()})
  
  # Tab - Static Plot -------------------------------------------------------
  output$PlotStatic <- renderImage({
    list(src = paste0("StaticPlots/",input$RegionID,".jpeg"))},
    deleteFile = FALSE)
  
  # Tab - Help --------------------------------------------------------------
  #Help image
  output$Help <- renderImage({
    list(src = "StaticPlots/Help.jpg")
  }, deleteFile = FALSE)
  
  # Tab - Interactive plot --------------------------------------------------
#   #Interactive Plot
#   plotObject <- reactive({source("source_ggplot.R",local=TRUE)})
#   output$PlotInteractive <- renderPlot({
#     #show progress message
#     withProgress(session, {
#       setProgress(message = "Please wait...")
#       print(plotObject())
#       setProgress(message = "Please wait... Done!")
#     })})
  # Tab - Data --------------------------------------------------------------
#   output$DataSNP <- renderDataTable({
#     #dat_region()
#     xData <- 
#       dat_region()[ position >= input$Xrange[1] &
#                       position <= input$Xrange[2],
#                     c("snp_name","position","Allele1","Allele2","Freq1",
#                          "Effect","P_value_Log","HitSNP","SNPType"),with=FALSE]
#     setnames(xData,c("SNP","POS","A1","A2","Frq",
#                      "Effect","-Log10(Pvalue)","Hit","SNPType"))
#     xData$Hit <- ifelse(xData$Hit==1,"Yes","No")
#     xData$SNPType <- ifelse(xData$SNPType==1,"Typed","Imputed")
#     return(xData)
#   },options=list(iDisplayLength=10))
#   
#   output$DataSNPLD <- renderDataTable({
#     yData <- LD_1KG()[,c("SNP_A","BP_A","SNP_B","BP_B","R2"),with=FALSE]
#   },options=list(iDisplayLength=10))
#   
#   
  
  # Output to a file --------------------------------------------------------
  # Get the selected download file type.
#   downloadPlotType <- reactive({input$downloadPlotType})
#   
#   observe({
#     plotType    <- input$downloadPlotType
#     plotTypePDF <- plotType %in% c("pdf","svg")
#     plotUnit    <- ifelse(plotTypePDF, "inches", "pixels")
#     plotUnitDefHight <- ifelse(plotTypePDF, 10, 800)
#     plotUnitDefWidth <- ifelse(plotTypePDF, 10, 800)
#     
#     updateNumericInput(
#       session,
#       inputId = "downloadPlotHeight",
#       label = sprintf("Height (%s)", plotUnit),
#       value = plotUnitDefHight)
#     
#     updateNumericInput(
#       session,
#       inputId = "downloadPlotWidth",
#       label = sprintf("Width (%s)", plotUnit),
#       value = plotUnitDefWidth)
#   })
#   
#   
#   # Get the download dimensions.
#   downloadPlotHeight <- reactive({
#     input$downloadPlotHeight
#   })
#   
#   downloadPlotWidth <- reactive({
#     input$downloadPlotWidth
#   })
#   
#   # Get the download file name.
#   downloadPlotFileName <- reactive({
#     input$downloadPlotFileName
#   })
#   
#   # Include a downloadable file of the plot in the output list.
#   output$downloadPlot <- downloadHandler(
#     filename = function() {
#       paste(downloadPlotFileName(), downloadPlotType(), sep=".")   
#     },
#     # The argument content below takes filename as a function
#     # and returns what's printed to it.
#     content = function(con) {
#       # Gets the name of the function to use from the 
#       # downloadFileType reactive element. Example:
#       # returns function pdf() if downloadFileType == "pdf".
#       plotFunction <- match.fun(downloadPlotType())
#       plotFunction(con, width = downloadPlotWidth(), height = downloadPlotHeight())
#       print(plotObject())
#       dev.off(which=dev.cur())
#     }
#   )
  
  
  
  
  #Show plot data
  #   output$metaDataTable <-
  #     renderDataTable({
  #       #if no LD seleceted then show all data else only in LD with selected SNP
  #       if(input$LDselect =="None"){
  #         d <- dat()[,c("snp_name","position","Freq1",
  #                       "P_value_Log","SNPType","R2")]}else{
  #                         d <- dat()[dat()$R2>=0.5,
  #                                    c("snp_name","position","Freq1",
  #                                      "P_value_Log","SNPType","R2")]
  #                       }#if else
  #       #add LNCAP annot
  #       d$LNCAP <- 
  #         ifelse(d$position %in% 
  #                  dat_LNCAP()[ dat_LNCAP()[,input$LNCAP]==1,"position"],
  #                input$LNCAP,"-")
  #       return(d)},
  #       options=list(
  #         #color sorted column
  #         bSortClasses=TRUE,
  #         #limit row number to display to 5
  #         aLengthMenu=c(10,50,100),iDisplayLength=10))
  
})#END shinyServer


