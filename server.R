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
      fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/",input$RegionID,"_stats.txt"),
            header = TRUE) %>%
        mutate(SNP = ifelse(is.na(SNP), rsid, SNP))
    } else {
      datAssoc() %>% filter(BP %in% unique(datLD()$BP_A)) %>%
        transmute(SNP,
                  Method = "glm",
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

    return(x %>%
             transmute(BP=BP,
                       Recomb=RECOMB,
                       RecombAdj=Recomb * ROIPLogMax() / 100))
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

#   # Define Zoom Start End ---------------------------------------------------
  zoomStart <- reactive({ input$BPrange[1] })
  zoomEnd <- reactive({ input$BPrange[2] })


#   # Data level 3 - Plot data ------------------------------------------------
  plotdatAssoc <- reactive({
    ROIdatAssoc() %>%
      filter(PLog >= input$FilterMinPlog)})

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
        plotdatAssoc(),
        d_LD[,c("BP_B","R2","LDSNP","LDSmoothCol","LDCol")],
        by.x="BP",by.y="BP_B",all=TRUE) %>%
      mutate(R2_Adj=ROIPLogMax()*R2)

    return(d_LD)
  })

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

#   # Output ------------------------------------------------------------------
#   # Output Summary ----------------------------------------------------------
  
  output$refProstatePaper <- renderUI(
    includeMarkdown(paste0("Data/ProstateData/",input$dataType,
                           "/README.md")))
  
  output$SummaryStats <- DT::renderDataTable({
    datAssoc() %>% arrange(P) %>%
      #if SNP name has rs number then convert to a link to NCBI
      mutate(SNP=ifelse(substr(SNP,1,2)=="rs",
                        paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                               gsub("rs","",SNP),'" target="_blank">',SNP,'</a>'),
                        SNP))},
    # FALSE to parse as a link
    escape = FALSE)

  output$SummaryHitSNPStats <- DT::renderDataTable({
   datStats() %>%
               dplyr::select(SNP, Method, BP, Stats) %>% arrange(Method,BP)
    })
    
    

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

   # Plot --------------------------------------------------------------------
   #plot title
  output$plotTitle <-
     renderUI(a(paste0(RegionChr(),':',zoomStart(),'-',zoomEnd()),
                href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                            RegionChr(),'%3A',zoomStart(),'-',zoomEnd()), target="_blank"))

  #Plot Chr ideogram
  plotObjChromosome <- reactive({source("Source/Chromosome.R",local=TRUE)})
  output$PlotChromosome <- renderPlot({print(plotObjChromosome())})
  #Manhattan track
  plotObjManhattan <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    plotManhattan(assoc = plotdatAssoc(), LD = plotDatLD(), geneticMap = plotDatGeneticMap())

    #source("Source/Manhattan.R",local=TRUE)
    })
  output$PlotManhattan <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    print(plotObjManhattan())})

  #SNPType track
  plotObjSNPType <- reactive({source("Source/SNPType.R",local=TRUE)})
  output$PlotSNPType <- renderPlot({print(plotObjSNPType())})
  #SNP LD track
  plotObjSNPLD <- reactive({source("Source/LD.R",local=TRUE)})
  output$PlotSNPLD <- renderPlot({print(plotObjSNPLD())})

  #wgEncodeBroadHistone 7 bigwig data track
  plotObjwgEncodeBroadHistone <- reactive({source("Source/wgEncodeBroadHistone.R",local=TRUE)})
  output$PlotwgEncodeBroadHistone <- renderPlot({print(plotObjwgEncodeBroadHistone())})

  #wgEncodeRegDnaseClustered 1 bed file
  plotObjwgEncodeRegDnaseClustered <- reactive({source("Source/wgEncodeRegDnaseClustered.R",local=TRUE)})
  output$PlotwgEncodeRegDnaseClustered <- renderPlot({print(plotObjwgEncodeRegDnaseClustered())})


  #Select hit SNPs
  output$HitSNPs <-
    renderUI({
      checkboxGroupInput("HitSNPs", h4("Hit SNPs:"),
                         RegionHits(),
                         #select max of 5 SNPs
                         #selected = 1:min(5,length(RegionHits()))

                         #if the input data oncoarray then select Forward hits only
                         selected = if(input$dataType == "OncoArray" &
                                       input$HitSNPsType == "Stepwise Forward"){
                           datStats() %>%
                             filter(Method == input$HitSNPsType) %>%
                             arrange(Stats) %>%
                             head(5) %>%
                             .$SNP

                           } else if(input$dataType == "OncoArray" &
                               input$HitSNPsType == "BVS"){
                            datStats() %>%
                               filter(Method == input$HitSNPsType) %>%
                               arrange(-Stats) %>%
                               head(5) %>%
                               .$SNP
                           } else {1:min(5,length(RegionHits()))}
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

