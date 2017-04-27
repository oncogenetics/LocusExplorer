# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Server file for shiny

# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {  
  
  # Data level 1 - Raw Input ------------------------------------------------
  datAssoc <- reactive({
    switch(input$dataType,
           iCOGS = {
             fread(paste0("Data/ProstateData/iCOGS/plotData/", input$RegionID, "_assoc.txt"),
                   header = TRUE, data.table = FALSE)},
           OncoArrayFineMapping = {
             fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/", input$RegionID, "_assoc.txt"),
                   header = TRUE, data.table = FALSE)},
           OncoArrayMeta = {
             fread(paste0("Data/ProstateData/OncoArrayMeta/plotData/", input$RegionID, "_assoc.txt"),
                   header = TRUE, data.table = FALSE)},
           Custom = {
             #input file check
             validate(need(input$FileStats != "", "Please upload Association file"))
             
             inFile <- input$FileStats
             
             if(is.null(inFile)){return(NULL)} else {
               fread(inFile$datapath, header = TRUE, data.table = FALSE) 
             }
             
           },
           Example = {
             fread("Data/CustomDataExample/Association.txt",
                   header = TRUE, data.table = FALSE)
           })
  }) # END datAssoc
  
  datLD <- reactive({
    switch(input$dataType,
           iCOGS = {
             fread(paste0("Data/ProstateData/iCOGS/plotData/",input$RegionID,"_LD.txt"),
                   header = TRUE, data.table = FALSE)},
           OncoArrayFineMapping = {
             fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/",input$RegionID,"_LD.txt"),
                   header = TRUE, data.table = FALSE)},
           OncoArrayMeta = {
             fread(paste0("Data/ProstateData/OncoArrayMeta/plotData/",input$RegionID,"_LD.txt"),
                   header = TRUE, data.table = FALSE)},
           Custom = {
             #input file check
             #validate(need(input$FileLD != "", "Please upload LD file"))
             
             # If the LD file is missing then create a dummy LD input, 
             # with top SNP LD at 0.01
             inFile <- input$FileLD
             if(is.null(inFile)){
               datAssoc() %>% 
                 transmute(CHR_A = RegionChrN(),
                           BP_A = datAssoc() %>% arrange(P) %>% head(1) %>% .$BP,
                           SNP_A = datAssoc() %>% arrange(P) %>% head(1) %>% .$SNP,
                           CHR_B = CHR_A,
                           BP_B = BP,
                           SNP_B = SNP,
                           R2 = 0.01)
             }else{fread(inFile$datapath, header = TRUE, data.table = FALSE) }
           },
           Example = {
             fread("Data/CustomDataExample/LD.txt",
                   header = TRUE, data.table = FALSE)}
    )# END switch
  })# END datLD
  
  datBedGraph <- reactive({
    switch(input$dataType,
           iCOGS = {
             res <-
               fread(paste0("Data/ProstateData/iCOGS/plotData/", input$RegionID,
                            "_EQTL.txt"), header = FALSE, data.table = FALSE)
             res <- res[,1:4]
             colnames(res) <- c("CHR", "START", "END", "SCORE")
             return(res)
           },
           OncoArrayFineMapping = { annotOncoFinemap },
           OncoArrayMeta = { annotOncoNewHits },
           Custom = {
             #input file check
             validate(need(input$FileBedGraph != "", "Please upload BedGraph file"))
             
             inFile <- input$FileBedGraph
             if(is.null(inFile)){return(NULL)} else {
               res <- fread(inFile$datapath, header = FALSE, data.table = FALSE)
               res <- res[,1:4]
               colnames(res) <- c("CHR", "START", "END", "SCORE")
               return(res)
             }
           },
           Example = {
             res <- fread("Data/CustomDataExample/bedGraph.txt",
                          header = FALSE, data.table = FALSE)
             res <- res[,1:4]
             colnames(res) <- c("CHR", "START", "END", "SCORE")
             return(res)
           })
  }) # END datBedGraph
  
  datLDlink <- reactive({
    #input file check
    validate(need(input$FileLDlink != "", "Please upload LDlink file"))
    
    inFile <- input$FileLDlink
    if(is.null(inFile)){return(NULL)}
    fread(inFile$datapath, header = TRUE, data.table = FALSE)
  })
  
  datLDlinkProcess <- reactive({
    
    datLDlink() %>%
      ## select only relevant columns
      dplyr::select(SNP_B = RS_Number, Coord, R2) %>%
      ## add index SNP for LD comparison
      mutate(
        SNP_B = ifelse(SNP_B == ".", Coord, SNP_B),
        SNP_A = SNP_B[c(1)]) %>%
      ## separate 'Coord' into CHR_B and BP_B
      separate(Coord, into = c("CHR_B", "BP_B"), sep = ":") %>%
      mutate(CHR_B = gsub("chr", "", x = CHR_B),
             ## add BP for index SNP
             CHR_A = CHR_B[c(1)],
             BP_A = BP_B[c(1)]) %>%
      ## Reorder columns
      dplyr::select(c(6, 7, 5, 2, 3, 1, 4)) %>%
      arrange(BP_B)})
  
  
  datStats <-  reactive({
    if(input$dataType == "OncoArrayFineMapping"){
      datAssoc() %>% filter(BP %in% unique(datLD()$BP_A)) %>%
        transmute(SNP,
                  Method = "JAM",
                  BP,
                  Stats = P)
      
      # fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/",input$RegionID,"_stats.txt"),
      #       header = TRUE) %>%
      #   mutate(SNP = ifelse(is.na(SNP), rsid, SNP))
    } else {
      datAssoc() %>% filter(BP %in% unique(datLD()$BP_A)) %>%
        transmute(SNP,
                  Method = "Forward regression",
                  BP,
                  Stats = P)
    }
  })
  
  # Define ROI --------------------------------------------------------------
  RegionFlank <- reactive({
    as.numeric(input$Flank)
  })
  
  RegionChr <- reactive({ datAssoc()$CHR[1] })
  RegionChrN <- reactive({
    x <- gsub("chr","",RegionChr())
    as.numeric(ifelse(x == "X", "23", x))})
  
  RegionStart <- reactive({
    x <- datAssoc() %>%
      filter(CHR == RegionChr() & !is.na(BP)) %>%
      .$BP %>% min
    #round down to 10K bp
    as.integer(max(0, round((x - RegionFlank())/RegionFlank()) * RegionFlank()))
  })
  RegionEnd <- reactive({
    x <- datAssoc() %>%
      filter(CHR == RegionChr() & !is.na(BP)) %>%
      .$BP %>% max
    #round up to 10K bp
    as.integer(round((x + RegionFlank())/RegionFlank()) * RegionFlank())
  })
  
  RegionHits <- reactive({
    d <- datLD() %>% .$SNP_A %>% unique %>% sort
    d[1:min(length(d),150)]
  })
  
  RegionHitsSelected <- reactive({input$HitSNPs})
  
  #Genomic ranges to subset bigwig data - wgEncodeBroadHistone
  RegionGR <-
    reactive({
      GRanges(seqnames = RegionChr(),
              IRanges(start = RegionStart(),
                      end = RegionEnd()))
    })
  
  #   # Data level 2 - ROI ------------------------------------------------------
  ROIdatAssoc <- reactive({ datAssoc() %>%
      filter(BP >= RegionStart() &
               BP <= RegionEnd()) %>%
      mutate(PLog = -log10(P)) })
  ROIdatLD <- reactive({ datLD() %>%
      filter(CHR_B == RegionChrN() &
               BP_B >= RegionStart() &
               BP_B <= RegionEnd()) })
  ROIdatBedGraph <- reactive({
    datBedGraph() %>%
      filter(CHR == RegionChr()) %>%
      #scale to -1 and 1,
      mutate(SCORE=SCORE/max(abs(SCORE),na.rm = TRUE))
  })
  
  ROIPLogMax <- reactive({
    #ylim Max for plotting Manhattan
    maxY <- max(ROIdatAssoc()$PLog,na.rm = TRUE)
    return(max(10,ceiling((maxY+1)/5)*5))
  })
  
  ROIdatGeneticMap <- reactive({
    mychr <- ifelse(RegionChr() == "chrX", "chr23", RegionChr())
    tabixRegion <- paste0(mychr,":",
                          RegionStart(), "-",
                          RegionEnd())
    x <- tabix.read.table("Data/GeneticMap1KG/GeneticMap1KG.txt.gz", tabixRegion)
    colnames(x) <- c("CHR", "BP", "RECOMB")
    
    return(x)
    # return(x %>%
    #          transmute(BP,
    #                    Recomb=RECOMB,
    #                    RecombAdj=Recomb * ROIPLogMax() / 100))
  })
  

  #   # Define Zoom Start End ---------------------------------------------------
  zoomStart <- reactive({ input$BPrange[1] })
  zoomEnd <- reactive({ input$BPrange[2] })
  
  
  #   # Data level 3 - Plot data ------------------------------------------------
  plotDatAssoc <- reactive({
    ROIdatAssoc() %>%
      filter(PLog >= input$FilterMinPlog)})
  
  plotDatLD <- reactive({
    #subset LD based on ui input
    LD <- ROIdatLD() %>%
      filter(
        R2 >= input$FilterMinLD &
          SNP_A %in% RegionHitsSelected())})

  #subset of recombination rates for zoomed region
  plotDatGeneticMap <- reactive({ ROIdatGeneticMap() %>%
      filter(BP >= zoomStart() &
               BP <= zoomEnd()) })
  
  #get granges collapsed genes for ggplot+ggbio
  plotDatGene <- reactive({
    udf_GeneSymbol(chrom = RegionChr(),
                   chromStart = zoomStart(),
                   chromEnd = zoomEnd()) })
  
  #number of genes in zoomed region
  plotDatGeneN <- reactive({
    res <- try({
      length(unique(plotDatGene()@elementMetadata$gene_id))}, silent=TRUE)
    if(class(res)=="try-error"){res <- 1}
    return(res)
  })
  
  # Output ------------------------------------------------------------------
  # Output Summary ----------------------------------------------------------
  
  # testing vars ------------------------------------------------------------
  output$testVars <- renderDataTable({
    data.frame(Data =
                 c("input$FilterMinPlog",
                   "input$FilterMinLD",
                   "RegionStart",
                   "RegionEnd",
                   "zoomStart", #reactive({ input$BPrange[1] })
                   "zoomEnd",   #reactive({ input$BPrange[2] })
                   "input$suggestiveLine, input$genomewideLine",
                   "input$Chr",
                   
                   "RegionHits",
                   "RegionHitsSelected",
                   
                   "datAssoc",
                   "ROIdatAssoc",
                   "plotDatAssoc",
                   "colnames(plotDatAssoc)",
                   
                   "datLD",
                   "ROIdatLD",
                   "plotDatLD",
                   "colnames(plotDatLD)",
                   
                   "ROIdatGeneticMap",
                   "plotDatGeneticMap",
                   "colnames(plotDatGeneticMap)"
                   ),
               Size = c(
                 input$FilterMinPlog,
                 input$FilterMinLD,
                 RegionStart(),
                 RegionEnd(),
                 zoomStart(),
                 zoomEnd(),
                 paste(input$suggestiveLine, input$genomewideLine, sep = ","),
                 input$Chr,
                 
                 paste(RegionHits(), collapse = ","),
                 paste(RegionHitsSelected(), collapse = ","),
                 
                 paste(dim(datAssoc()), collapse = ","),
                 paste(dim(ROIdatAssoc()), collapse = ","),
                 paste(dim(plotDatAssoc()), collapse = ","),
                 paste(colnames(plotDatAssoc()), collapse = ","),
                 
                 paste(dim(datLD()), collapse = ","),
                 paste(dim(ROIdatLD()), collapse = ","),
                 paste(dim(plotDatLD()), collapse = ","),
                 paste(colnames(plotDatLD()), collapse = ","),
                 
                 paste(dim(ROIdatGeneticMap()), collapse = ","),
                 paste(dim(plotDatGeneticMap()), collapse = ","),
                 paste(colnames(plotDatGeneticMap()), collapse = ",")
                 )
    )
  })
  
  
  
  
  output$refProstatePaper <- renderUI(
    includeMarkdown(paste0("Data/ProstateData/",input$dataType,
                           "/README.md")))
  
  output$SummaryStats <- DT::renderDataTable({
    datAssoc() %>% arrange(P) %>%
      #if SNP name has rs number then convert to a link to NCBI
      mutate(SNP = ifelse(substr(SNP,1,2) == "rs",
                          paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                 gsub("rs", "", SNP), '" target="_blank">', SNP, '</a>'),
                          SNP))},
    # FALSE to parse as a link
    escape = FALSE)
  
  output$SummaryHitSNPStats <- DT::renderDataTable({
    datStats() %>%
      dplyr::select(SNP, Method, BP, Stats) %>% arrange(Method,BP) %>% 
      mutate(SNP = ifelse(substr(SNP, 1, 2) == "rs",
                          paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                                 gsub("rs", "", SNP), '" target="_blank">', SNP, '</a>'),
                          SNP))
  },
  # FALSE to parse as a link
  escape = FALSE)
  
  
  
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
  }, escape = FALSE)
  
  output$SummaryRegion <-
    renderUI(a(paste0(RegionChr(),':',RegionStart(),'-',RegionEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',RegionStart(),'-',RegionEnd()),target="_blank"))
  

  # Plot: Title ---------------------------------------------------------------
  output$plotTitle <-
    renderUI(a(paste0(RegionChr(),':',zoomStart(),'-',zoomEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',zoomStart(),'-',zoomEnd()), target="_blank"))
  
  # Plot: Chr ideogram --------------------------------------------------------
  plotObjChromosome <- reactive({source("Source/Chromosome.R", local = TRUE)})
  output$PlotChromosome <- renderPlot({print(plotObjChromosome())})
  
  # Plot: Manhattan ------------------------------------------------------------
  plotObjManhattan <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    plotManhattan(assoc = plotDatAssoc(),
                  LD = plotDatLD(),
                  geneticMap = plotDatGeneticMap(),
                  suggestiveLine = input$suggestiveLine,
                  genomewideLine = input$genomewideLine,
                  hits = RegionHitsSelected(),
                  xStart = zoomStart(),
                  xEnd = zoomEnd(),
                  opts = intersect(
                    c("Hits", "LD", "SuggestiveLine", "GenomewideLine", 
                      input$ShowHideTracks), 
                    c("Recombination", "LDSmooth",
                      "Hits", "LD", "SuggestiveLine", "GenomewideLine"))
                  ) + theme_LE()
                  #opts = c("LD", "LDSmooth", "Hits"))
    
    #source("Source/Manhattan.R",local=TRUE)
  })
  output$PlotManhattan <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjManhattan())})
  
  # Plot: SNPType -------------------------------------------------------------
  # SNP is imputed or typed
  plotObjSNPType <- reactive({
    plotSNPtype(assoc = plotDatAssoc(),
                xStart = zoomStart(),
                xEnd = zoomEnd()) + theme_LE()
    
    #source("Source/SNPType.R",local=TRUE)
    })
  output$PlotSNPType <- renderPlot({print(plotObjSNPType())})
  
  # Plot: SNP LD --------------------------------------------------------------
  plotObjSNPLD <- reactive({
    plotLD(data = plotDatLD(),
           hits = RegionHitsSelected(),
           xStart = zoomStart(),
           xEnd = zoomEnd()) + theme_LE()
    #source("Source/LD.R",local=TRUE)
    })
  output$PlotSNPLD <- renderPlot({print(plotObjSNPLD())})
  
  # Plot: wgEncodeBroadHistone 7 bigwig ---------------------------------------
  plotObjwgEncodeBroadHistone <- reactive({
    plotHistone(folder = "Data/wgEncodeBroadHistone/",
                chr = input$Chr,
                xStart = zoomStart(),
                xEnd = zoomEnd()
                ) + theme_LE()
    
    #source("Source/wgEncodeBroadHistone.R",local=TRUE)
    })
  output$PlotwgEncodeBroadHistone <- renderPlot({print(plotObjwgEncodeBroadHistone())})
  

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
                         selected = RegionHits()[1:min(c(5, length(RegionHits())))]
                         
      ) #END checkboxGroupInput
    })
  
  
  output$Chr <- renderUI({
    selectInput(inputId = "Chr", label = h5("Chr"),
                choices = regions %>%
                  filter(DATA == input$dataType) %>%
                  .$CHR %>% unique,
                selected = "chr2"
    )
  })
  
  output$RegionID <- renderUI({
    selectInput(inputId = "RegionID", label= h5("Region ID"),
                choices = regions %>%
                  filter(CHR == input$Chr &
                           DATA == input$dataType) %>%
                  .$REGIONBED,
                #pre selected a region to demonstrate as an example - HNF1B gene
                selected = "chr2_241657087_242920971"
    )})
  
})#END shinyServer

