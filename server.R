# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt


# About -------------------------------------------------------------------
# Server file for shiny

# Workspace ---------------------------------------------------------------

# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {  
  
  source("source/UDF.R",local = TRUE)
  
  #prostate region names
  regions <- read.csv("Data/ProstateData/regions.csv") %>% 
    mutate(REGIONBED=paste0(CHR,"_",START,"_",END))
  
  #wgEncodeBroadHistone bigwig data description
  wgEncodeBroadHistoneFileDesc <- 
    read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
             stringsAsFactors = FALSE)
  

  # Data level 1 - Raw Input ------------------------------------------------
  datStats <- reactive({
    switch(input$dataType,
           Prostate = {
             fread(paste0("Data/ProstateData/LE/",input$RegionID,"_assoc.txt"),
                   header=TRUE, data.table=FALSE)
             },
           Custom = {
             #input file check
             validate(need(input$FileStats != "", "Please upload Association file"))
             
             inFile <- input$FileStats
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE) 
           },
           Example = {
             fread("Data/CustomDataExample/Association.txt",
                   header=TRUE, data.table=FALSE)
           })
  })
  
  datLD <- reactive({
    switch(input$dataType,
           Prostate = {
             #input file check
             validate(need(input$RegionID != "", "Please select RegionID"))
             fread(paste0("Data/ProstateData/LE/",input$RegionID,"_LD.txt"),
                    header=TRUE, data.table=FALSE)
             },
           Custom = {
             #input file check
             #validate(need(input$FileLD != "", "Please upload LD file"))
             
             # If the LD file is missing then create a dummy LD input, 
             # with top SNP LD at 0.01
             inFile <- input$FileLD
             if(is.null(inFile)){
               datStats() %>% 
                 transmute(CHR_A=RegionChrN(),
                           BP_A=datStats() %>% arrange(P) %>% head(1) %>% .$BP,
                           SNP_A=datStats() %>% arrange(P) %>% head(1) %>% .$SNP,
                           CHR_B=CHR_A,
                           BP_B=BP,
                           SNP_B=SNP,
                           R2=0.01)
               }else {fread(inFile$datapath, header=TRUE, data.table=FALSE) }
           },
           Example = {
             fread("Data/CustomDataExample/LD.txt",
                   header=TRUE, data.table=FALSE) 
           })
  })
  
  datLNCAP <- reactive({
    fread("Data/ProstateLNCAP/ProstateLNCAP.txt",
          colClasses = c("character","numeric","character")) })
  
  datEQTL <- reactive({
    switch(input$dataType,
           Prostate = {
             #input file check
             validate(need(input$RegionID != "", "Please select RegionID"))
             fread(paste0("Data/ProstateData/LE/",input$RegionID,"_EQTL.txt"),
                   header=TRUE, data.table=FALSE)
             },
           Custom = {
             #input file check
             validate(need(input$FileEQTL != "", "Please upload eQTL file"))
             
             inFile <- input$FileEQTL
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE)
           },
           Example = {
             fread("Data/CustomDataExample/eQTL.txt",
                   header=TRUE, data.table=FALSE) 
           })
  })

  datLDlink <- reactive({
    #input file check
    validate(need(input$FileLDlink != "", "Please upload LDlink file"))
    
    inFile <- input$FileLDlink
    if(is.null(inFile)){return(NULL)}
    fread(inFile$datapath, header=TRUE, data.table=FALSE)
    })
  
  datLDlinkProcess <- reactive({
    #x <- fread("../CIDR_Data/CIDR_UTAH/chr8_rs17446916.txt", header=TRUE, data.table=FALSE)
    #x %>% head %>%  
    
    datLDlink() %>% 
    ## select only relevant columns
    dplyr::select(SNP_B = RS_Number,Coord,R2) %>% 
      ## add index SNP for LD comparison
      mutate(
        SNP_B = ifelse(SNP_B==".",Coord,SNP_B),
        SNP_A = SNP_B[c(1)]) %>% 
      ## separate 'Coord' into CHR_B and BP_B
      separate(Coord, into = c("CHR_B", "BP_B"), sep = ":") %>% 
      mutate(CHR_B=gsub("chr", replacement = "", x = CHR_B),
             ## add BP for index SNP
             CHR_A = CHR_B[c(1)],
             BP_A = BP_B[c(1)]) %>%
      ## Reorder columns
      dplyr::select(c(6,7,5,2,3,1,4)) %>% 
      arrange(BP_B)
    })
  
  # Define ROI --------------------------------------------------------------
  RegionFlank <- reactive({
    if(input$Flank==""){rf <- 10000}else{rf <- input$Flank}
    as.numeric(rf)
  })
  RegionChr <- reactive({ datStats()$CHR[1] })
  RegionChrN <- reactive({
    x <- gsub("chr","",RegionChr())
    as.numeric(ifelse(x=="X","23",x))
    })
  RegionStart <- reactive({ 
    x <- datStats() %>% 
      filter(CHR==RegionChr() &
               !is.na(BP)) %>% 
      .$BP %>% min 
    #round down to 10K bp
    as.integer(max(0,round((x-RegionFlank())/RegionFlank()) * RegionFlank()))
    #as.integer(max(0,round((x-10000)/10000) * 10000))
    })
  RegionEnd <- reactive({ 
    x <- datStats() %>% 
      filter(CHR==RegionChr() &
               !is.na(BP)) %>% 
      .$BP %>% max 
    #round up to 10K bp
    as.integer(round((x+RegionFlank())/RegionFlank()) * RegionFlank())
    #as.integer(round((x+10000)/10000) * 10000)
    })
  
  RegionHits <- reactive({ 
    d <- datLD() %>% .$SNP_A %>% unique %>% sort
    d[1:min(length(d),15)] })
  
  RegionHitsSelected <- reactive({input$HitSNPs})
  
  #Genomic ranges to subset bigwig data - wgEncodeBroadHistone
  RegionGR <- 
    reactive({
      GRanges(seqnames = RegionChr(),
              IRanges(start = RegionStart(),
                      end = RegionEnd())) 
      })
  
  # Data level 2 - ROI ------------------------------------------------------
  ROIdatStats <- reactive({ datStats() %>% 
      filter(BP>=RegionStart() &
               BP<=RegionEnd()) %>% 
      mutate(PLog=-log10(P)) })
  ROIdatLD <- reactive({ datLD() %>% 
      filter(CHR_B==RegionChrN() &
               BP_B>=RegionStart() &
               BP_B<=RegionEnd()) })
  ROIdatLNCAP <- reactive({ datLNCAP() %>% 
      filter(CHR==RegionChr() &
               BP>=RegionStart() &
               BP<=RegionEnd()) })
  ROIdatEQTL <- reactive({ datEQTL() %>% 
      filter(CHR==RegionChr() &
               START>=RegionStart() &
               END<=RegionEnd()) })
  
  ROIPLogMax <- reactive({
    #ylim Max for plotting Manhattan
    maxY <- max(ROIdatStats()$PLog,na.rm = TRUE)
    return(max(10,ceiling((maxY+1)/5)*5))
  })
  
  ROIdatGeneticMap <- reactive({ 
    mychr <- ifelse(RegionChr()=="chrX","chr23",RegionChr())
    tabixRegion <- paste0(mychr,":",
                          RegionStart(),"-",
                          RegionEnd())
    x <- tabix.read.table("Data/GeneticMap1KG/GeneticMap1KG.txt.gz",tabixRegion)
    colnames(x) <- c("CHR","BP","RECOMB")
    
    return(x %>% 
             transmute(BP=BP,
                       Recomb=RECOMB,
                       RecombAdj=Recomb*ROIPLogMax()/100))
  })
  
  ROIdatwgEncodeRegDnaseClustered <- reactive({ 
    #define region to subset
    tabixRegion <- paste0(RegionChr(),":",
                          RegionStart(),"-",
                          RegionEnd())
    #subset using seqmineR::tabix.read.table
    x <- tabix.read.table("Data/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.txt.gz",tabixRegion)
    colnames(x) <- c("CHR","START","END","SCORE")
    return(x)
  })
  
  # Define Zoom Start End ---------------------------------------------------
  zoomStart <- reactive({ input$BPrange[1] })
  zoomEnd <- reactive({ input$BPrange[2] })
  
  
  # Data level 3 - Plot data ------------------------------------------------
  plotDatStats <- reactive({
    ROIdatStats() %>% 
      filter(BP>=zoomStart() &
               BP<=zoomEnd() &
               PLog >= input$FilterMinPlog)})
  
  plotDatLD <- reactive({
    #subset LD based on ui input
    LD <- ROIdatLD() %>% 
      filter(
        R2 >= input$FilterMinLD &
          SNP_A %in% RegionHitsSelected())
    
    # http://stackoverflow.com/a/8197703/680068
    # hcl(h=seq(15, 375, length=6), l=65, c=100)[1:5]
    # colors match first 5 colours of default ggplot2
    colourLD5 <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
    #create pallete LD 0 to 100
    colLD <- lapply(colourLD5,function(i)colorRampPalette(c("grey95",i))(100))
    
    #assign colours to LD 
    d_LD <- 
      base::do.call(
        rbind,
        lapply(RegionHitsSelected(), function(snp){
          d <- LD %>% filter(SNP_A == snp)
          #LD round minimum is 1 max 100
          LDColIndex <- ifelse(round(d$R2,2)==0,1,round(d$R2,2)*100)
          LDColIndex <- ifelse(LDColIndex>100,100,LDColIndex)
          d$LDSNP <- snp
          d$LDSmoothCol <- colLD[[match(snp,RegionHitsSelected())]][100]
          d$LDCol <- colLD[[match(snp,RegionHitsSelected())]][LDColIndex]
          return(d)
        }))
    
    # to add smooth LD Y value used for Pvalue and LD 0-1
    # add pvalues for Y value on the plot
    d_LD <- 
      base::merge.data.frame(
        plotDatStats(),
        d_LD[,c("BP_B","R2","LDSNP","LDSmoothCol","LDCol")],
        by.x="BP",by.y="BP_B",all.x=TRUE) %>% 
      mutate(R2_Adj=ROIPLogMax()*R2)
    
    return(d_LD)
  })
  
  #subset of recombination rates for zoomed region
  plotDatGeneticMap <- reactive({ ROIdatGeneticMap() %>% 
      filter(BP>=zoomStart() &
               BP<=zoomEnd()) })
  
  #get granges collapsed genes for ggplot+ggbio
  plotDatGene <- reactive({ 
    udf_GeneSymbol(chrom=RegionChr(),
                   chromStart=zoomStart(),
                   chromEnd=zoomEnd()) })
  
  #number of genes in zoomed region
  plotDatGeneN <- reactive({ 
    res <- try({
      length(unique(plotDatGene()@elementMetadata$gene_id))}, silent=TRUE)
    if(class(res)=="try-error"){res <- 1}
    return(res)
    })
  
  #wgEncodeBroadHistone data 7 big wig data
  plotDatwgEncodeBroadHistone <- reactive({
    #get file list matching Encode description file
    bigWigFiles <- 
      intersect(wgEncodeBroadHistoneFileDesc$File,
                list.files("Data/wgEncodeBroadHistone/","*.bigWig"))
    
    #check if bigWig files are downloaded, if plot warning message.
    if(length(bigWigFiles)>0){
      #subset region 
      wgEncodeBroadHistone <- 
        rbind_all(
          #seven bigwig files
          lapply(bigWigFiles, function(i){
            #i=1  wgEncodeBroadHistoneFileDesc$File[1]
            as.data.frame(
              import(paste0("Data/wgEncodeBroadHistone/",i),
                     which=RegionGR())) %>% 
              filter(score >= 5) %>% 
              transmute(BP=start,
                        ENCODE=wgEncodeBroadHistoneFileDesc[
                          which(wgEncodeBroadHistoneFileDesc$File==i),
                          "Name"],
                        SCORE=round(ifelse(score >= 100, 100, score),0))
          }))
      
      #merge on BP, to plot overlappling
      res <- data.frame(BP=unique(wgEncodeBroadHistone$BP))
      output <- 
        rbind_all(
          lapply(unique(wgEncodeBroadHistone$ENCODE), function(i){
            d <- merge(res,wgEncodeBroadHistone %>% filter(ENCODE==i),by="BP",all.x=TRUE)
            d$ENCODE <- i
            d$SCORE[ is.na(d$SCORE) ] <- 0
            return(d)
          })
        )
      }else{
      #if there are no bigwig files return blank data.frame
      output <- data.frame()
    }
    return(output)
    })
  
  #output of ENCODE data
  #wgEncodeBroadHistone data
  output$SummarywgEncodeBroadHistone <- 
    DT::renderDataTable(plotDatwgEncodeBroadHistone() %>% arrange(BP),
                        options=list(searching=FALSE,searchable=FALSE))
  #wgEncodeRegDnaseClustered
  output$SummarywgEncodeRegDnaseClustered <-
    DT::renderDataTable(ROIdatwgEncodeRegDnaseClustered() %>% arrange(START),
                        options=list(searching=FALSE,searchable=FALSE))
  
  # Output ------------------------------------------------------------------
  # Output Summary ----------------------------------------------------------
  output$SummaryStats <- DT::renderDataTable( 
    datStats() %>% arrange(P) %>% 
      #if SNP name has rs number then convert to a link to NCBI
      mutate(SNP=ifelse(substr(SNP,1,2)=="rs",
                        paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                               gsub("rs","",SNP),'" target="_blank">',SNP,'</a>'),
                        SNP)),
    # FALSE to parse as a link
    escape=FALSE)
  
  output$SummaryLD <- DT::renderDataTable({ 
    datLD() %>% 
      arrange(BP_B) %>% 
      mutate(SNP_A=ifelse(substr(SNP_A,1,2)=="rs",
                          paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                 gsub("rs","",SNP_A),'" target="_blank">',SNP_A,'</a>'),
                          SNP_A),
             SNP_B=ifelse(substr(SNP_B,1,2)=="rs",
                          paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                 gsub("rs","",SNP_B),'" target="_blank">',SNP_B,'</a>'),
                          SNP_B))
    },escape=FALSE)
  
  
  output$SummaryLNCAP <- 
    DT::renderDataTable(ROIdatLNCAP() %>% arrange(BP),
                        options=list(searching=FALSE,searchable=FALSE))
  output$SummaryEQTL <- 
    DT::renderDataTable(datEQTL() %>% arrange(START),
                        options=list(searching=FALSE,searchable=FALSE))
  
  output$SummaryRegion <- 
    renderUI(a(paste0(RegionChr(),':',RegionStart(),'-',RegionEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',RegionStart(),'-',RegionEnd()),target="_blank"))
  
  output$SummaryHits <- DT::renderDataTable( 
    data.table(SNP = RegionHits()) %>% 
      #if SNP name has rs number then convert to a link to NCBI
      mutate(SNP=ifelse(substr(SNP,1,2)=="rs",
                        paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                               gsub("rs","",SNP),'" target="_blank">',SNP,'</a>'),
                        SNP)),
    escape=FALSE,options=list(paging=FALSE,searching=FALSE,searchable=FALSE))

#   output$SummaryFileNrowNcol <- DT::renderDataTable( 
#     data.table(InputFile=c("Association","LD","LNCAP","eQTL"),
#                RowCount=c(nrow(datStats()),nrow(datLD()),
#                       nrow(datLNCAP()),nrow(datEQTL())),
#                ColumnCount=c(ncol(datStats()),ncol(datLD()),
#                       ncol(datLNCAP()),ncol(datEQTL()))),
#     options=list(paging=FALSE,searching=FALSE,searchable=FALSE))
  
  output$SummaryLDlink <- DT::renderDataTable(datLDlink())
  output$SummaryLDlinkProcess <- 
    DT::renderDataTable({
      datLDlinkProcess() %>% 
        arrange(BP_B) %>% 
        mutate(SNP_A=ifelse(substr(SNP_A,1,2)=="rs",
                            paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                   gsub("rs","",SNP_A),'" target="_blank">',SNP_A,'</a>'),
                            SNP_A),
               SNP_B=ifelse(substr(SNP_B,1,2)=="rs",
                            paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                   gsub("rs","",SNP_B),'" target="_blank">',SNP_B,'</a>'),
                            SNP_B))
    },escape=FALSE)
  
  # Plot --------------------------------------------------------------------
  #plot title
  output$plotTitle <- 
     renderUI(a(paste0(RegionChr(),':',zoomStart(),'-',zoomEnd()),
                href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                            RegionChr(),'%3A',zoomStart(),'-',zoomEnd()), target="_blank"))

  #Plot Chr ideogram
  plotObjChromosome <- reactive({source("source/Chromosome.R",local=TRUE)})
  output$PlotChromosome <- renderPlot({print(plotObjChromosome())})
  #Manhattan track
  plotObjManhattan <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    source("source/Manhattan.R",local=TRUE)})
  output$PlotManhattan <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjManhattan())})
  
  #SNPType track
  plotObjSNPType <- reactive({source("source/SNPType.R",local=TRUE)})
  output$PlotSNPType <- renderPlot({print(plotObjSNPType())})
  #SNP LD track
  plotObjSNPLD <- reactive({source("source/LD.R",local=TRUE)})
  output$PlotSNPLD <- renderPlot({print(plotObjSNPLD())})
  
  #wgEncodeBroadHistone 7 bigwig data track
  plotObjwgEncodeBroadHistone <- reactive({source("source/wgEncodeBroadHistone.R",local=TRUE)})
  output$PlotwgEncodeBroadHistone <- renderPlot({print(plotObjwgEncodeBroadHistone())})
  
  #wgEncodeRegDnaseClustered 1 bed file
  plotObjwgEncodeRegDnaseClustered <- reactive({source("source/wgEncodeRegDnaseClustered.R",local=TRUE)})
  output$PlotwgEncodeRegDnaseClustered <- renderPlot({print(plotObjwgEncodeRegDnaseClustered())})
  
  
  #LNCAP Smooth track
  plotObjLNCAP <- reactive({source("source/LNCAP.R",local=TRUE)})
  output$PlotLNCAP <- renderPlot({print(plotObjLNCAP())})
  #eQTL bar track
  plotObjEQTL <- reactive({source("source/eQTL.R",local=TRUE)})
  output$PlotEQTL <- renderPlot({print(plotObjEQTL())})
  
  #Gene track
  plotObjGene <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    source("source/Gene.R",local=TRUE)})
  output$PlotGene <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjGene())})
  # Example plot ------------------------------------------------------------
  output$ExamplePlotJPEG <- renderImage({
    return(list(
        src = "www/chr17_36020000_36140000.jpeg",
        contentType = "image/jpeg",
        alt = "ExamplePlotOutput"
      ))
    }, deleteFile = FALSE)
  # Plot Merge --------------------------------------------------------------
  #Dynamic size for tracks
  RegionHitsCount <- reactive({ length(RegionHitsSelected()) })
  RegionGeneCount <- reactive({ plotDatGeneN() })
  RegionSNPTypeCount <- reactive({ length(unique(plotDatStats()$TYPED)) })
  
  #Default size per track
  trackSize <- reactive({ 
    data.frame(Track=c("Chromosome","Manhattan","LD","SNPType",
                       "wgEncodeBroadHistone","wgEncodeRegDnaseClustered",
                       "LNCAP","eQTL","Gene"),
               Size=c(100,400,
                      RegionHitsCount()*20,
                      RegionSNPTypeCount()*20,
                      60, # wgEncodeBroadHistone
                      20, # wgEncodeRegDnaseClustered
                      20, # lncap
                      30, # eqtl
                      RegionGeneCount()*30)) })

  #Create subset based on selected tracks
  trackHeights <- reactive({
    trackSize() %>% filter(Track %in% input$ShowHideTracks) %>% .$Size })
  trackColours <- reactive({ 
    res <- cbind(trackHeights(),c('grey90','grey80'))[,2][1:length(trackHeights())] 
    #chromosome background must be white
    if("Chromosome" %in% input$ShowHideTracks){res[1] <- "white"}
    return(res)
  })
  
  #output plot
  # See Dynamic UI section
  
  # Output to a file --------------------------------------------------------
  # Get the selected download file type.
  downloadPlotType <- reactive({input$downloadPlotType})
  
  observe({
    plotType    <- input$downloadPlotType
    plotTypePDF <- plotType %in% c("pdf","svg")
    plotUnit    <- ifelse(plotTypePDF, "inches", "pixels")
    plotUnitDefHeight <- ifelse(plotTypePDF, 12, 1200)
    plotUnitDefWidth <- ifelse(plotTypePDF, 10, 1000)
    
    updateNumericInput(
      session,
      inputId = "downloadPlotHeight",
      label = sprintf("Height (%s)", plotUnit),
      value = plotUnitDefHeight)
    
    updateNumericInput(
      session,
      inputId = "downloadPlotWidth",
      label = sprintf("Width (%s)", plotUnit),
      value = plotUnitDefWidth)
  })
  
  
  # Get the download dimensions.
  downloadPlotHeight <- reactive({
    input$downloadPlotHeight
  })
  
  downloadPlotWidth <- reactive({
    input$downloadPlotWidth
  })
  
  # Get the download file name.
  downloadPlotFileName <- reactive({
    input$downloadPlotFileName
  })
  
  # Get the download resolution for jpeg/tiff
  downloadPlotResolution <- reactive({
    input$downloadPlotResolution
  })
  
  # Include a downloadable file of the plot in the output list.
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(downloadPlotFileName(), downloadPlotType(), sep=".")   
    },
    # The argument content below takes filename as a function
    # and returns what's printed to it.
#     content = function(con) {
#       # Gets the name of the function to use from the 
#       # downloadFileType reactive element. Example:
#       # returns function pdf() if downloadFileType == "pdf".
#       plotFunction <- match.fun(downloadPlotType())
#       plotFunction(con, width = downloadPlotWidth(), height = downloadPlotHeight())
#       print(plotObjMerge())
#       dev.off(which=dev.cur())
#     }
  
    content = function(con) {
      # Gets the name of the function to use from the 
      # downloadFileType reactive element. Example:
      # returns function pdf() if downloadFileType == "pdf".
      plotFunction <- match.fun(downloadPlotType())
      
      switch(input$downloadPlotType,
             pdf = {plotFunction(con, 
                                 width = downloadPlotWidth(),
                                 height = downloadPlotHeight(),
                                 useDingbats = FALSE)},
             svg = {plotFunction(con, 
                                 width = downloadPlotWidth(),
                                 height = downloadPlotHeight())},
             jpeg = {plotFunction(con, 
                                  width = downloadPlotWidth(),
                                  height = downloadPlotHeight(),
                                  type="cairo",
                                  res=downloadPlotResolution(),
                                  quality = 100)},
             tiff = {plotFunction(con, 
                                  width = downloadPlotWidth(),
                                  height = downloadPlotHeight(),
                                  type="cairo",
                                  res=downloadPlotResolution())}
             
             )
      
      print(plotObjMerge())
      dev.off(which=dev.cur())
      }) #downloadHandler downloadPlot
  
  #download LDlink processed data as flat file tab separated
  output$downloadLDFile <- downloadHandler(
    filename = function() { paste0(datLDlinkProcess()[1,"SNP_A"],"_LD.txt") },
    content = function(file) {
      write.table(datLDlinkProcess(), file,
                  sep="\t",row.names=FALSE, quote=FALSE)
    })
  
  # Dynamic UI --------------------------------------------------------------
  #Zoom to region X axis BP
  output$BPrange <-
    renderUI({
      sliderInput("BPrange", h5("Region: Start-End"),
                  min = RegionStart(),
                  max = RegionEnd(),
                  value = c(RegionStart(),RegionEnd()),
                  step = 20000)})
  #Select hit SNPs
  output$HitSNPs <- 
    renderUI({
      checkboxGroupInput("HitSNPs", h4("Hit SNPs:"),
                         RegionHits(),
                         #select max of 5 SNPs
                         selected = 1:min(5,length(RegionHits())))
      })
  
  output$RegionID <- renderUI({
    selectInput(inputId="RegionID", label= h5("RegionID"), 
                choices=regions %>% filter(CHR==input$Chr) %>% .$REGIONBED,
                #pre selected a region to demonstrate as an example - HNF1B gene
                selected="chr17_35947000_36204000")})
  
  
  #download file name - default: chr_start_end
  output$downloadPlotFileName <- renderUI({
    textInput(
      inputId = "downloadPlotFileName",
      label = h4("Download file name"),
      value = paste(RegionChr(),zoomStart(),zoomEnd(),sep="_"))})
  #download plot title - default: chr_start_end
  output$downloadPlotTitle <- renderUI({
    textInput(
      inputId = "downloadPlotTitle",
      label = h4("Plot title"),
      value = paste(RegionChr(),zoomStart(),zoomEnd(),sep="_"))})
  
  # Merged Final plot -------------------------------------------------------  
  #merged plot with dynamic plot height
  plotObjMerge <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    source("Source/MergePlots.R",local=TRUE)

    })
  output$plotMerge <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjMerge())
    
    })
  
  #plot merge height is dynamic, based on seleceted tracks
  output$plotMergeUI <- renderUI({
    plotOutput("plotMerge",width=800,height=sum(trackHeights()))
    })
  
  # Observe update ----------------------------------------------------------
  # maximum of 5 SNPs can be selected to display LD, minimum 1 must be ticked.
  observe({
    if(length(input$HitSNPs) > 5){
      updateCheckboxGroupInput(session, "HitSNPs", selected= head(input$HitSNPs,5))}
    if(length(input$HitSNPs) < 1){
      updateCheckboxGroupInput(
        session, "HitSNPs",
        selected= RegionHits()[1:min(5,length(RegionHits()))])}
    if(length(input$ShowHideTracks) < 1){
      updateCheckboxGroupInput(session, "ShowHideTracks", selected= "Manhattan")}
  })
  observeEvent(input$RegionZoom,({
    if(!input$RegionZoom %in% c("chr:start-end","")){
      newStartEnd <- as.numeric(unlist(strsplit(input$RegionZoom,":|-"))[2:3])
      updateSliderInput(session, "BPrange",
                        value = newStartEnd)}
    }))
  #Reset plot options
  observeEvent(input$resetInput,({
    updateSliderInput(session,"FilterMinPlog", value=0)
    updateSliderInput(session,"FilterMinLD", value=0)
    updateSliderInput(session,"BPrange", value = c(RegionStart(),RegionEnd()))
    updateTextInput(session,"RegionZoom",value="chr:start-end")
    updateSelectInput(session,"Flank",
                      choices = list("10KB"=10000,
                                     "50KB"=50000,
                                     "1MB"=100000,
                                     "2MB"=200000), 
                      selected = 10000)
    updateCheckboxGroupInput(session,"HitSNPs",
                             choices = RegionHits(),
                             selected = 
                               RegionHits()[1:min(5,length(RegionHits()))])
    
    updateCheckboxGroupInput(session,"ShowHideTracks",
                             choices = c("Chromosome"="Chromosome",
                                         "Manhattan"="Manhattan",
                                         "LDSmooth"="LDSmooth",
                                         "LD"="LD",
                                         "SNPType"="SNPType",
                                         "wgEncodeBroadHistone"="wgEncodeBroadHistone",
                                         "wgEncodeRegDnaseClustered"="wgEncodeRegDnaseClustered",
                                         "LNCaP Prostate"="LNCAP",
                                         "eQTL"="eQTL",
                                         "Gene"="Gene"),
                             selected=c("Manhattan")
                             #selected=c("Manhattan","LD","LDSmooth")
                             )
    })) #END observeEvent resetInput
  
})#END shinyServer

