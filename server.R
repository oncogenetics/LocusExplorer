# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Server file for shiny

# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {  
  
  # Data level 1 - Raw Input ------------------------------------------------
  datAssoc <- reactive({
    req(input$dataType)
    switch(input$dataType,
           OncoArrayFineMapping = {
             req(input$RegionID)
             fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/", input$RegionID, "_assoc.txt"),
                   header = TRUE, data.table = FALSE)},
           # OncoArrayMeta = {
           #   fread(paste0("Data/ProstateData/OncoArrayMeta/plotData/", input$RegionID, "_assoc.txt"),
           #         header = TRUE, data.table = FALSE)},
           Custom = {
             #input file check
             validate(need(input$FileStats != "", "Please upload Association file"))
             
             inFile <- input$FileStats
             req(inFile)
             fread(inFile$datapath, header = TRUE, data.table = FALSE) 
           },
           Example = {
             fread("Data/CustomDataExample/Association.txt",
                   header = TRUE, data.table = FALSE)
           })
  }) # END datAssoc
  
  datLD <- reactive({
    req(input$dataType)
    switch(input$dataType,
           OncoArrayFineMapping = {
             req(input$RegionID)
             fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/",input$RegionID,"_LD.txt"),
                   header = TRUE, data.table = FALSE)},
           # OncoArrayMeta = {
           #   fread(paste0("Data/ProstateData/OncoArrayMeta/plotData/",input$RegionID,"_LD.txt"),
           #         header = TRUE, data.table = FALSE)},
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
  
  datStats <-  reactive({
    req(input$dataType)
    req(input$RegionID)
    
    if(input$dataType == "OncoArrayFineMapping"){
      
      d <- fread(paste0("Data/ProstateData/OncoArrayFineMapping/plotData/",
                        input$RegionID,"_stats.txt")) 
      #   SNP,rsid,BP,PP_best_tag,PP_tag_r2,PostProb,BF,JAM99
      #   rs587636640,1:150283370:C:AA,150283370,1:150283370:C:AA,1,0.0129,2.601395291,0
      
      d$Method <- "JAM"
      d$Stats <- as.numeric(d$BF)
      
      #return  
      d
    } else {
      datAssoc() %>% dplyr::filter(BP %in% unique(datLD()$BP_A)) %>%
        transmute(SNP,
                  Method = "Forward regression",
                  BP,
                  Stats = P)
    }
  })
  
  datAnnot <-   reactive({
    annotOncoFinemap
    # req(input$dataType)
    # switch(input$dataType,
    #        OncoArrayFineMapping = { annotOncoFinemap },
    #        # OncoArrayMeta = { NULL },
    #        Custom = { NULL },
    #        Example = { NULL })
  }) # END datAnnot
  
  datAnnotEQTL <-   reactive({
    annotOncoFinemapEQTL
    # req(input$dataType)
    # switch(input$dataType,
    #        OncoArrayFineMapping = { annotOncoFinemapEQTL },
    #        # OncoArrayMeta = { NULL },
    #        Custom = { NULL },
    #        Example = { NULL })
  }) # END datAnnotEQTL
  
  datLDlink <- reactive({
    #input file check
    validate(need(input$FileLDlink != "", "Please upload LDlink file"))
    
    inFile <- input$FileLDlink
    req(inFile)
    #if(is.null(inFile)){return(NULL)}
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
  
  # Define ROI --------------------------------------------------------------
  RegionFlank <- reactive({ req(input$Flank)
    max(c(1, input$Flank * 1000)) })
  
  RegionChr <- reactive({ 
    req(datAssoc())
    datAssoc()$CHR[1] 
    })
  RegionChrN <- reactive({
    req(RegionChr())
    x <- gsub("chr", "", RegionChr(), fixed = TRUE)
    if(x == "X") x <- "23"
    as.numeric(x)
    })
  
  RegionStart <- reactive({
    req(datAssoc())
    req(RegionChr())
    x <- datAssoc() %>%
      dplyr::filter(CHR == RegionChr() & !is.na(BP)) %>%
      .$BP %>% min
    #round down to 10K bp
    as.integer(max(0, round((x - RegionFlank())/RegionFlank()) * RegionFlank()))
  })
  RegionEnd <- reactive({
    req(datAssoc())
    req(RegionChr())
    x <- datAssoc() %>%
      dplyr::filter(CHR == RegionChr() & !is.na(BP)) %>%
      .$BP %>% max
    #round up to 10K bp
    as.integer(round((x + RegionFlank())/RegionFlank()) * RegionFlank())
  })
  
  RegionHits <- reactive({
    req(datLD())
    d <- datLD() %>% .$SNP_A %>% unique %>% sort
    d[1:min(length(d), 150)] 
  })
  
  RegionHitsSelected <- reactive({
    x <- unique(c(input$HitSNPs, input$otherHits))
    
    if(is.null(x)) x <- "" 
       
    x   
    })
  
  #Genomic ranges to subset bigwig data - wgEncodeBroadHistone
  RegionGR <-
    reactive({
      GRanges(seqnames = RegionChr(),
              IRanges(start = RegionStart(),
                      end = RegionEnd()))
    })
  
  #   # Data level 2 - ROI ------------------------------------------------------
  ROIdatAssoc <- reactive({
    req(datAssoc())
    datAssoc() %>%
      dplyr::filter(BP >= RegionStart() &
               BP <= RegionEnd()) %>%
      mutate(PLog = -log10(P)) })
  ROIdatLD <- reactive({
    req(datLD())
    datLD() %>%
      dplyr::filter(CHR_B == RegionChrN() &
                      BP_B >= RegionStart() & BP_B <= RegionEnd()) })
  # ROIdatBedGraph <- reactive({
  #   datBedGraph() %>%
  #     dplyr::filter(CHR == RegionChr()) %>%
  #     #scale to -1 and 1,
  #     mutate(SCORE = SCORE/max(abs(SCORE), na.rm = TRUE))
  # })
  
  ROIdatAnnot <- reactive({
    req(datAnnot())
    datAnnot() %>%
      dplyr::filter(CHR == RegionChrN() &
                      BP >= RegionStart() & BP <= RegionEnd())
  })
  
  ROIdatAnnotEQTL <- reactive({
    req(datAnnotEQTL())
    datAnnotEQTL() %>%
      dplyr::filter(CHR == RegionChrN() &
                      (
                        (SNPhit_BP >= RegionStart() & SNPhit_BP <= RegionEnd()) |
                          (SNP_BP >= RegionStart() & SNP_BP <= RegionEnd())
                      )
                    )
  })
  
  ROIPLogMax <- reactive({
    #ylim Max for plotting Manhattan
    maxY <- max(ROIdatAssoc()$PLog, na.rm = TRUE)
    max(10, ceiling((maxY+1)/5)*5)
  })
  
  ROIdatGeneticMap <- reactive({
    req(GeneticMap1KG)
    GeneticMap1KG[ CHR == RegionChrN() & BP >= RegionStart() & BP <= RegionEnd()]
  })
  
  
  # Define Zoom Start End ---------------------------------------------------
  zoomStart <- reactive({ input$BPrange[1] })
  zoomEnd <- reactive({ input$BPrange[2] })
  
  
  # Data level 3 - Plot data ------------------------------------------------
  plotDatAssoc <- reactive({
    d <- ROIdatAssoc()
    d <- d[ d$PLog >= input$FilterMinPlog, ]
    d[ order(d$BP), ] 
    })
  
  plotDatLD <- reactive({
    #subset LD based on ui input
    LD <- ROIdatLD() %>%
      dplyr::filter(
        R2 >= input$FilterMinLD &
          SNP_A %in% RegionHitsSelected())})
  
  #subset of recombination rates for zoomed region
  plotDatGeneticMap <- reactive({ ROIdatGeneticMap() %>%
      dplyr::filter(between(BP, zoomStart(), zoomEnd())) })
  
  #get granges collapsed genes for ggplot+ggbio
  plotDatGene <- reactive({
    GeneSymbol(chrom = RegionChr(),
               chromStart = zoomStart(),
               chromEnd = zoomEnd()) })
  
  #number of genes in zoomed region
  # if not genes, set to 1, to plot a message: no genes.
  plotDatGeneN <- reactive({
    res <- try({
      length(unique(plotDatGene()@elementMetadata$gene_id))}, silent = TRUE)
    if(class(res) == "try-error"){res <- 1}
    #return
    res })
  
  # OncoArray finemapping Annotaion and TCGA EQTL
  plotDatAnnot <- reactive({ annotOncoFinemap[ CHR == RegionChrN(), ] })
  plotDatAnnotEQTL <- reactive({
    # c(CHR, BP, TYPE1, TYPE2, COLOUR_HEX, TYPE2N) 
    annotOncoFinemapEQTL %>%
      filter(CHR == RegionChrN() &
               between(SNP_BP, zoomStart(), zoomEnd())) %>% 
      #filter(SNPhit == SNP) %>% 
      #filter(SNP_BP %in% datLD()[datLD()$SNP_A %in% RegionHitsSelected(), "BP_A"]) %>% 
      transmute(CHR,
                BP = SNP_BP,
                TYPE1 = "EQTL",
                TYPE2 = GENE,
                COLOUR_HEX = "#F00034", #ICRcolours("BrighRed")
                TYPE2N = 4) %>%
      unique()
    })
  
    
    
  # Output ------------------------------------------------------------------
  # Output Summary ----------------------------------------------------------
  
  output$SummaryStats <- DT::renderDataTable({
    datAssoc() %>% arrange(P) %>%
      #if SNP name has rs number then convert to a link to NCBI
      mutate(SNP = hyperlink(SNP))},
    # FALSE to parse as a link
    escape = FALSE)
  
  output$SummaryHitSNPStats <- DT::renderDataTable({
    datStats() %>%
      dplyr::select(SNP, Method, BP, Stats) %>%
      arrange(Method, BP) %>% 
      mutate(
        # to output Inf as "Inf" convert to character, relevant GitHub issue:
        # https://github.com/jeroen/jsonlite/issues/94
        # https://stackoverflow.com/q/38771807/680068
        Stats = as.character(Stats),
        # output SNP names as links to NCBI
        SNP = hyperlink(SNP))
  },
  # FALSE to parse as a link
  escape = FALSE)
  
  
  output$SummaryLD <- DT::renderDataTable({
    req(input$dataType)
    datLD() %>%
      arrange(BP_B) %>%
      mutate(SNP_A = hyperlink(SNP_A),
             SNP_B = hyperlink(SNP_B))
  }, escape = FALSE)
  
  output$SummaryROIdatAnnot <- DT::renderDataTable({
    merge(datAssoc()[, c("SNP", "BP")], ROIdatAnnot(), by = "BP") %>% 
      transmute(SNP = hyperlink(SNP),
                BP,
                TYPE1, TYPE2) %>% 
      arrange(BP)
    
    }, escape = FALSE)
  
  output$SummaryROIdatAnnotEQTL <- DT::renderDataTable({
    ROIdatAnnotEQTL() %>% 
      mutate(SNPhit = hyperlink(SNPhit),
             SNP  = hyperlink(SNP),
             R2,
             GENE = hyperlink(GENE, type = "geneSymbol"))
    }, escape = FALSE)
  
  # output$SummaryBedGraph <- 
  #   DT::renderDataTable({datBedGraph() %>% arrange(START)},
  #                       options = list(searching = FALSE, searchable = FALSE))
  
  output$ui_SummaryRegion <-
    renderUI(a(paste0(RegionChr(),':',RegionStart(),'-',RegionEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',RegionStart(),'-',RegionEnd()),target="_blank"))
  
  
  # ~~~ TAB: Manhattan --------------------------------------------------------
  
  # Plot: Title ---------------------------------------------------------------
  output$ui_plotTitle <-
    renderUI(a(paste0(RegionChr(),':',zoomStart(),'-',zoomEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',zoomStart(),'-',zoomEnd()), target="_blank"))
  
  # Plot: Chr ideogram --------------------------------------------------------
  plotObjChromosome <- reactive({
    Ideogram(genome = "hg19",
             subchr = RegionChr(),
             #highlight ROI - region
             color = "#F00034", fill = "#F00034",
             zoom.region = c(RegionStart(), RegionEnd())) +
      ylab(label = element_blank())
    })
  output$PlotChromosome <- renderPlot({print(plotObjChromosome())})
  
  # Plot: Manhattan Pvalues ----------------------------------------------------
  plotObjManhattanPvalues <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    plotHits <- plotDatAssoc()[ plotDatAssoc()$SNP %in% RegionHitsSelected(), "SNP" ]
    
    plotManhattan(assoc = plotDatAssoc(),
                  LD = plotDatLD(),
                  geneticMap = plotDatGeneticMap(),
                  suggestiveLine = input$suggestiveLine,
                  genomewideLine = input$genomewideLine,
                  hits = plotHits,
                  xStart = zoomStart(),
                  xEnd = zoomEnd(),
                  opts = intersect(
                    c("Hits", "LD", "SuggestiveLine", "GenomewideLine", 
                      input$ShowHideManhattanPvalues), 
                    c("Recombination", "LDSmooth",
                      "Hits", "LD", "SuggestiveLine", "GenomewideLine", "Effect")),
                  pad = TRUE,
                  postprob = FALSE) +
      theme_LE() 
      
  })
  output$PlotManhattanPvalues <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjManhattanPvalues())})
  
  
  # Plot: Manhattan PostProbs --------------------------------------------------
  plotObjManhattanPostProbs <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # c("SNP","BP","P","TYPED") 
    datStats <- datStats()
    datStats$TYPED <- ifelse(datStats$SNP %in% datAssoc()[datAssoc()$TYPED == 2, "SNP"], 2, 1)
    plotDat <- data.frame(SNP = datStats$SNP,
                          BP = datStats$BP,
                          P = datStats$PostProb,
                          TYPED = datStats$TYPED,
                          BF = datStats$BF)
    plotDat <- plotDat[ order(plotDat$BP), ]
    # add BF to snp names, e.g: rs1234 (12.34)
    plotHits <- plotDat[ plotDat$SNP %in% RegionHitsSelected(), "SNP" ]
    plotHitsName <- paste0(plotHits, " (",
                           formatC(plotDat[ plotDat$SNP %in% RegionHitsSelected(), "BF"],
                                   digits = 2, format = "f"),
                          ")")
    
    # x <- c(0.1111, -Inf, Inf, 23439924.22222, -0.00000001, 0.000002)
    # format(round(x, 2), nsmall = 2)
    # formatC(x, digits = 2, format = "f")
    
    plotManhattan(assoc = plotDat,
                  LD = plotDatLD(),
                  geneticMap = plotDatGeneticMap(),
                  suggestiveLine = 0,
                  genomewideLine = 0,
                  hits = plotHits,
                  #hitsName = plotHitsName,
                  xStart = zoomStart(),
                  xEnd = zoomEnd(),
                  opts = intersect(
                    c("Hits", "LD", "SuggestiveLine", "GenomewideLine", 
                      input$ShowHideManhattanPostProbs), 
                    c("Recombination", "LDSmooth",
                      "Hits", "LD", "SuggestiveLine", "GenomewideLine")),
                  pad = TRUE,
                  postprob = TRUE) +
      theme_LE()
  })
  output$PlotManhattanPostProbs <- renderPlot({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Almost there...", value = 80)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    print(plotObjManhattanPostProbs())})
  
  # Plot: SNPType -------------------------------------------------------------
  # SNP is imputed or typed
  plotObjSNPType <- reactive({
    plotSNPtype(assoc = plotDatAssoc(),
                xStart = zoomStart(),
                xEnd = zoomEnd(),
                pad = TRUE) +
      theme_LE()
    
  })
  output$PlotSNPType <- renderPlot({print(plotObjSNPType())})
  
  # Plot: SNP LD --------------------------------------------------------------
  plotObjSNPLD <- reactive({
    plotLD(data = plotDatLD(),
           hits = RegionHitsSelected(),
           xStart = zoomStart(),
           xEnd = zoomEnd(),
           pad = TRUE) +
      theme_LE()
  })
  output$PlotSNPLD <- renderPlot({print(plotObjSNPLD())})
  
  # Plot: wgEncodeBroadHistone 7 bigwig ---------------------------------------
  plotObjwgEncodeBroadHistone <- reactive({
    
    gg <- try({
      plotHistone(folder = "Data/wgEncodeBroadHistone/",
                  chr = RegionChr(),
                  xStart = zoomStart(),
                  xEnd = zoomEnd(),
                  pad = TRUE
      ) + theme_LE()}, silent = TRUE)
    
    if(class(gg) == "try-error"){ 
      gg <- plotBlank(zoomStart(), zoomEnd(),
                      yLabel = expression(ENCODE[]),
                      textAnnot = "Error: Histone bigwig files are missing.") +
        theme_LE() #+
        #ylab(expression(ENCODE[]))
      }
    
    #return 
    gg
    })
  
  output$PlotwgEncodeBroadHistone <- renderPlot({print(plotObjwgEncodeBroadHistone())})
  
  # Plot: BedGraph ------------------------------------------------------------
  # plotObjBedGraph <- reactive({
  #   plotBlank(xStart = 1, xEnd = 10, yLabel = "test") + theme_LE()
  #   
  # })
  # output$PlotBedGraph <- renderPlot({print(plotObjBedGraph())})
  
  # Plot: Annot & eQTL OncoFinemap -------------------------------------------
  plotObjAnnotOncoFinemap <- reactive({
    plotDat <- rbind(plotDatAnnot(), plotDatAnnotEQTL())
  
    
    gg <- plotAnnot(plotDat,
                    chrom = RegionChrN(),
                    xStart = zoomStart(), xEnd = zoomEnd(),
                    collapse = TRUE) +
      theme_LE()
    # add gene names if there eqtl
    # if(nrow(plotDatAnnotEQTL()) > 0) {
    #   gg <- gg +
    #     geom_text_repel(data = plotDatAnnotEQTL(),
    #                     aes(x = BP, y = 4, label = TYPE2), col = "black")}
    #return ggplot
    gg
    
  })
  output$PlotAnnotOncoFinemap <- renderPlot({print(plotObjAnnotOncoFinemap())})
  
  
  
  # Plot: Gene ---------------------------------------------------------------
  # returns ggplot and count of genes
  plotObjGene <- reactive({
    plotGene(chrom = RegionChr(),
             chromStart = zoomStart(), chromEnd = zoomEnd(),
             hits = unique(plotDatAnnotEQTL()$TYPE2),
             vline = unique(plotDatLD()$BP_A),
             pad = TRUE)})
  # ggplot object to plot gene track
  plotObjGenePlot <- reactive({ 
    plotObjGene()$genePlot + theme_LE()})
  
  # numeric count of genes, passed for merging tracks, more genes vertical space.
  # used for cowplot::plot_grid(), rel_heights
  # plotObjGene()$geneCnt
  
  output$PlotGene <- renderPlot({print(plotObjGenePlot())})
  
  # ~~~ TAB: LD-Heatmap ------------------------------------------------------
  plotObjSNP_LDheatmap <- reactive({
    plotLDmatrix(data = plotDatLD(), hits = RegionHitsSelected())
  })
  
  output$plotSNP_LDheatmap <- renderPlot({print(plotObjSNP_LDheatmap())})
  
  # ~~~ TAB: LD-Network -----------------------------------------------------
  
  # 1. Nodes ------------------------------------------------------------------
  nodes <- reactive({
    nodeHits <-  data.frame(id = RegionHitsSelected(),
                            label = RegionHitsSelected(),
                            group = "Hits")
    
    nodeTags <- 
      rbind.data.frame(
        plotDatLD() %>% 
          transmute(
            id = SNP_A,
            label = SNP_A,
            group = "Tags"),
        plotDatLD() %>% 
          transmute(
            id = SNP_B,
            label = SNP_B,
            group = "Tags")) %>% 
      filter(!(id %in% nodeHits$id)) %>% 
      unique 
    # %>% 
    #   filter(
    #     # filter on BF
    #     id %in% metaRegionClean[ MaxBF > input$filterMinBF, SNP]
    #   )
    
    #merge hit nodes and other tag nodes
    if(nrow(nodeTags) > 0) {
      res <- rbind.data.frame(nodeHits, nodeTags)} else {
        res <- nodeHits}
    
    #keep only hits?
    if(input$hitsOnly){
      res <- res %>%
        filter(id %in% RegionHitsSelected())
    }
    
    # remove nodes that have no links
    res <- res %>%
      filter(id %in% RegionHitsSelected() |
               id %in% c(links()$from, links()$to))
    
    #add meta info: pvalue, infoscore, allele, maf
    # res <- merge(res,
    #              metaRegionClean[, list(SNPid, MaxBF, title)],
    #              by = "SNPid", all.x = TRUE)
    
    
    # # set the size of nodes
    # if(input$nodeSizeBF){
    #   res <- res %>% 
    #     mutate(size = as.numeric(
    #       as.character(cut(MaxBF,
    #                        breaks = c(-1, 0, 1,  3,  5, 10, 100, Inf),
    #                        labels =    c(1, 4, 10, 12, 15, 20, 25)))))
    #   #return
    #   res %>%
    #     arrange(group) %>% 
    #     select(id, label, group, size, title) %>% unique
    # } else {
    #   #return
    #   res %>%
    #     arrange(group) %>% 
    #     select(id, label, group, title) %>% unique
    #   
    # }
    
    #return
    res %>%
      arrange(group) %>% 
      dplyr::select(id, label, group) %>%
      unique
    
  })
  
  # 2. Links ---------------------------------------------------------------------
  links <- reactive({
    if(input$hitsOnly){
      res <-
        plotDatLD() %>%
        transmute(from = SNP_A,
                  to = SNP_B,
                  value = R2)%>%
        filter(from != to) %>%
        filter(from %in% RegionHitsSelected() &
                 to %in% RegionHitsSelected())
    } else {
      res <-
        plotDatLD() %>%
        transmute(from = SNP_A,
                  to = SNP_B,
                  value = R2) %>%
        filter(from != to &
                 from %in% RegionHitsSelected())
    }
    
    #res$color <- if_else(res$value < 0.1, "gray", "lightblue")
    res$color <- "lightblue"
    #res$value <- round(res$value, 2) * 100
    res$title <- paste0("R2:", round(res$value, 2))
    
    res %>%
      filter(value >= input$FilterMinLD)
    
  })
  # 3. networkHits ----------------------------------------------------------
  output$plotSNP_LDnetwork <- renderVisNetwork({
    visNetwork(nodes(), links(),
               main = "LD network",
               submain = paste0("Filter R2>", input$FilterMinLD)
               ) %>%
      visOptions(highlightNearest = TRUE) %>%
      visLegend(useGroups = FALSE,
                addNodes = data.frame(label = c("Hit SNPs", "Tag SNPs"),
                                      shape = c("circle", "circle"),
                                      color = c("#97C2FC", "yellow")),
                addEdges = data.frame(label = "LD",
                                      color = "lightblue",
                                      arrows = "none", width = 2),
                zoom = FALSE) %>% 
      visInteraction(navigationButtons = TRUE) %>% 
      visIgraphLayout(layout = input$networkLayout, randomSeed = 12) %>% 
      visExport(name = paste(RegionChr(), zoomStart(), zoomEnd(), "network", sep = "_"),
                label = "Save as PNG",
                style = "background-color:#85E7FF")
  })
  
  # Dynamic UI --------------------------------------------------------------
  output$ui_refProstatePaper <- renderUI(
    if(input$dataType %in% c("OncoArrayFineMapping",
                             "OncoArrayMeta","iCOGS")){
      includeMarkdown(
        paste0("Data/ProstateData/", input$dataType, "/README.md"))
      } else if(input$dataType == "Custom") {
        includeMarkdown("Data/CustomDataExample/README.md")
        }
    )
  
  #Zoom to region X axis BP
  output$ui_BPrange <-
    renderUI({
      req(RegionStart())
      req(RegionEnd())
      sliderInput("BPrange", h5("Use sliders to zoom in to required region."),
                  min = RegionStart(),
                  max = RegionEnd(),
                  value = c(RegionStart(),RegionEnd()),
                  step = 20000)})
  
  #Select hit SNPs
  output$ui_HitSNPs <-
    renderUI({
      list(tags$div(align = 'left', 
                    class = 'multicol', 
                    checkboxGroupInput(inputId  = 'HitSNPs', 
                                       label    = "Hit SNPs:", 
                                       choices  = RegionHits(),
                                       selected = 
                                         if(input$dataType == "OncoArrayFineMapping") {
                                           datStats() %>% dplyr::filter(JAM99 == 1) %>% .$SNP
                                           # SNP,rsid,BP,PP_best_tag,PP_tag_r2,PostProb,BF,JAM99
                                           # rs6724057,chr2_242113751_A_G,242113751,chr2_242113751_A_G,1,0.0507,25.58232942,0
                                           # rs77559646,chr2_242135265_A_G,242135265,chr2_242135265_A_G,1,0.9998,Inf,1
                                         } else {
                                           RegionHits()[1:min(c(5, length(RegionHits())))]},
                                       inline   = FALSE))) 
      
      
      # checkboxGroupInput("HitSNPs", h4("Hit SNPs:"),
      #                    RegionHits(),
      #                    #select max of 5 SNPs
      #                    selected = RegionHits()[1:min(c(5, length(RegionHits())))],
      #                    inline = TRUE
      
      #) #END checkboxGroupInput
    })
  output$ui_otherHits <- renderUI({
    selectizeInput("otherHits", "Other SNPs:",
                   choices = plotDatAssoc()$SNP,
                   multiple = TRUE)
  })
  
  output$ui_Chr <- renderUI({
    selectInput(inputId = "Chr", label = h5("Chr"),
                choices = regions %>%
                  dplyr::filter(DATA == input$dataType) %>%
                  .$CHR %>% unique,
                selected = "chr2"
    )
  })
  
  output$ui_RegionID <- renderUI({
    req(input$Chr)
    req(input$dataType)
    selectInput(inputId = "RegionID", label= h5("Region ID"),
                choices = regions %>%
                  dplyr::filter(CHR == input$Chr &
                           DATA == input$dataType) %>%
                  .$REGIONBED,
                #pre selected a region to demonstrate as an example - HNF1B gene
                selected = "chr2_241657087_242920971"
    )})
  
  #download file name - default: chr_start_end
  output$ui_downloadPlotFileName <- renderUI({
    textInput(
      inputId = "downloadPlotFileName",
      label = "Download file name",
      value = paste(RegionChr(), zoomStart(), zoomEnd(), sep = "_"))})
  #download plot title - default: chr:start-end
  output$ui_downloadPlotTitle <- renderUI({
    textInput(
      inputId = "downloadPlotTitle",
      label = "Plot title",
      value = paste0(RegionChr(), ':', zoomStart(), '-', zoomEnd(), sep=""))})
  
  
  output$ui_PlotThemeColour2 <- renderUI({
    req(input$PlotThemeColour1)
    x <- input$PlotThemeColour1
    colourInput("PlotThemeColour2",
                "Plot theme shade 2",
                #"#E5E5E5"
                shadeTintColour(input$PlotThemeColour1, 10)
                )})
  
  # Merged Final plot -------------------------------------------------------  
  #merged plot with dynamic plot height
  plotList <- reactive({
    plots <- list(
      if("Chromosome" %in% input$ShowHideTracks) plotObjChromosome() else NA,
        #scale_y_continuous(expand = c(0.05, 0)) 
      
      if("Manhattan" %in% input$ShowHideManhattanPvalues) plotObjManhattanPvalues() else NA,
      if("Manhattan" %in% input$ShowHideManhattanPostProbs) plotObjManhattanPostProbs() else NA,
      if("SNPType" %in% input$ShowHideTracks) plotObjSNPType() + 
        scale_y_continuous(breaks = c(0.5, 1.5),
                           labels = strPadLeft(c("Imputed", "Typed")),
                           limits = c(-0.90, 2),
                           name = expression(SNP[])) else NA,
      if("LD" %in% input$ShowHideTracks) plotObjSNPLD() else NA,
      if("wgEncodeBroadHistone" %in% input$ShowHideTracks) plotObjwgEncodeBroadHistone() else NA,
      if("annotOncoFinemap" %in% input$ShowHideTracks) plotObjAnnotOncoFinemap() else NA,
      if("Gene" %in% input$ShowHideTracks) plotObjGenePlot() else NA
    )
    
    trackHeights <- c(
      100, # Chromosome
      300, # ManhattanPvalues
      120, # ManhattanPostProbs
      50,  # SNPType
      20 * length(input$HitSNPs), # LD
      60,  # wgEncodeBroadHistone
      80,  # annotOncoFinemap,
      min(c(240, 30 * plotObjGene()$geneCnt)) # Gene
    )
    
    ixKeep <- !is.na(plots)
    # return list of plots and heights
    list(Tracks = plots[ ixKeep ],
         Heights = trackHeights[ ixKeep ])
  })
  
  trackColours <- reactive({ 
    myCols <- 
      if(input$PlotTheme == "1"){
        c(input$PlotThemeColour1, input$PlotThemeColour2)
        } else { c("#FFFFFF", "#FFFFFF") }

    res <- cbind(plotList()$Heights, myCols)[, 2][ 1:length(plotList()$Heights) ]  
    #chromosome background must be white
    if("Chromosome" %in% input$ShowHideTracks){res <- head(c("white", res), -1)}
    return(res)
  })
    
  plotObjMerge <- reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Plotting please wait...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())

    tracks(
      plotList()$Tracks,
      heights = plotList()$Heights,
      padding = unit(-0.2, "lines"),
      track.plot.color = trackColours(),
      title = if(!is.null(input$downloadPlotTitle) &
                 nchar(input$downloadPlotTitle) > 0) { input$downloadPlotTitle
        } else { NULL }
      ) +
      if(input$PlotTheme == "2") {
        theme(panel.border = element_rect(fill = NA, color = 'grey80')) }
      

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
  output$ui_plotMergeUI <- renderUI({
    plotOutput("plotMerge", width = 800, height = sum(plotList()$Heights))
    #plotOutput("plotMerge", width = 800, height = 1000)
  })
  
  # Output to a file --------------------------------------------------------
  # Get the selected download file type.
  downloadPlotType <- reactive({input$downloadPlotType})
  
  observe({
    plotType    <- input$downloadPlotType
    plotTypePDF <- plotType %in% c("pdf","svg")
    plotUnit    <- ifelse(plotTypePDF, "inches", "pixels")
    plotUnitDefHeight <- ifelse(plotTypePDF, 12, 1200)
    plotUnitDefWidth <- ifelse(plotTypePDF, 10, 1000)
    #plotUnitDefMin <- ifelse(plotTypePDF, 5, 800)
    #plotUnitDefMax <- ifelse(plotTypePDF, 50, 10000)
    plotUnitDefStep <- ifelse(plotTypePDF, 1, 100)
    
    updateSliderInput( session,
                       inputId = "downloadPlotHeight",
                       label = sprintf("Height (%s)", plotUnit),
                       value = plotUnitDefHeight,
                       step = plotUnitDefStep)
    
    updateSliderInput(
      session,
      inputId = "downloadPlotWidth",
      label = sprintf("Width (%s)", plotUnit),
      value = plotUnitDefWidth,
      step = plotUnitDefStep)
  })
  
  
  downloadPlotHeight <- reactive({ input$downloadPlotHeight })
  downloadPlotWidth <- reactive({ input$downloadPlotWidth })
  downloadPlotFileName <- reactive({ input$downloadPlotFileName })
  downloadPlotRes <- reactive({ input$downloadPlotRes })
  downloadPlotPaper <- reactive({ input$downloadPlotPaper })
  downloadPlotType <- reactive({ input$downloadPlotType })
  downloadPlotPointSize <- reactive({ input$downloadPlotPointSize })
  downloadPlotQuality <- reactive({ input$downloadPlotQuality })
  downloadPlotTypeJPEG <- reactive({ input$downloadPlotTypeJPEG })
  downloadPlotTypeCompression <- reactive({ input$downloadPlotCompression })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(downloadPlotFileName(), downloadPlotType(), sep = ".")   
    },
    # The argument content below takes filename as a function
    # and returns what's printed to it.
    content = function(con) {
      # Gets the name of the function to use from the 
      # downloadFileType reactive element. Example:
      # returns function pdf() if downloadFileType == "pdf".
      plotFunction <- match.fun(downloadPlotType())
      
      switch(input$downloadPlotType,
             pdf = {plotFunction(con, 
                                 width = downloadPlotWidth(),
                                 height = downloadPlotHeight(),
                                 pointsize = downloadPlotPointSize(),
                                 paper = downloadPlotPaper(),
                                 useDingbats = FALSE)},
             svg = {plotFunction(con, 
                                 width = downloadPlotWidth(),
                                 height = downloadPlotHeight(),
                                 pointsize = downloadPlotPointSize())},
             jpeg = {plotFunction(con, 
                                  width = downloadPlotWidth(),
                                  height = downloadPlotHeight(),
                                  pointsize = downloadPlotPointSize(),
                                  res = downloadPlotRes(),
                                  quality = downloadPlotQuality(),
                                  type = downloadPlotTypeJPEG())},
             tiff = {plotFunction(con, 
                                  width = downloadPlotWidth(),
                                  height = downloadPlotHeight(),
                                  pointsize = downloadPlotPointSize(),
                                  res = downloadPlotRes(),
                                  type = downloadPlotTypeJPEG(),
                                  compression = downloadPlotTypeCompression())}
             
      )
      
      print(plotObjMerge())
      dev.off(which=dev.cur())
    }) #downloadHandler downloadPlot
  
  # Observe update ----------------------------------------------------------
  # maximum of 5 SNPs can be selected to display LD, minimum 1 must be ticked.
  observe({
    # if(length(input$HitSNPs) < 1){
    #   updateCheckboxGroupInput(
    #     session, "HitSNPs",
    #     #selected = RegionHits()[1:min(5, length(RegionHits()))]
    #     selected = RegionHits()[1]
    #     )}
    if(length(c(input$ShowHideTracks,
                input$ShowHideManhattanPvalues,
                input$ShowHideManhattanPostProbs)) < 1){
      updateCheckboxGroupInput(session, "ShowHideManhattanPvalues",
                               selected = "Manhattan")}
    if(input$dataType %in% c("Custom", "Example")){

      updateCheckboxGroupInput(session, "ShowHideManhattanPostProbs",
                               label = "Manhattan: PostProbs",
                               choices = c("Manhattan" = "Manhattan",
                                           "Recombination" = "Recombination",
                                           "LD smooth" = "LDSmooth"),
                               selected = NULL)
    }
  })
  # 
  # #manual chr:start-end zoom values
  # observeEvent(input$RegionZoom,({
  #   if(!input$RegionZoom %in% c("chr:start-end","")){
  #     newStartEnd <- as.numeric(unlist(strsplit(input$RegionZoom,":|-"))[2:3])
  #     updateSliderInput(session, "BPrange",
  #                       value = newStartEnd)}
  # }))
  
  #Reset plot options - 2.Plot Settings
  observeEvent(input$resetInput,({
    updateSliderInput(session,"FilterMinPlog", value=0)
    updateSliderInput(session,"FilterMinLD", value=0.2)
    
    updateNumericInput(session,"suggestiveLine", value=5)
    updateNumericInput(session,"genomewideLine", value=8)
    
    updateCheckboxInput(session,"adjustLabels", value=TRUE)
    updateSliderInput(session,"repFact", value=5)
    
    updateSliderInput(session,"BPrange", value = c(RegionStart(),RegionEnd()))
    updateTextInput(session,"RegionZoom",value="chr:start-end")
    updateSelectInput(session,"Flank",
                      selected = 10000)
    updateCheckboxGroupInput(session,"HitSNPs",
                             choices = RegionHits(),
                             selected =
                               RegionHits()[1:min(5,length(RegionHits()))])
    
    updateCheckboxGroupInput(session,"ShowHideTracks",
                             selected = c("Manhattan", "Recombination"))
  })
  ) #END observeEvent resetInput
  
  #reset plot output file settings - 3.Final Plot
  observeEvent(input$resetDownloadPlotSettings |
                 input$downloadPlotAdvancedSettings,({
                   updateSelectInput(session, "downloadPlotType", selected = "jpeg")
                   updateNumericInput(session, "downloadPlotWidth", value = 1000)
                   updateNumericInput(session, "downloadPlotHeight", value = 1200)
                   updateSliderInput(session, "downloadPlotPointSize", value = 12)
                   updateSelectInput(session, "downloadPlotPaper", selected = "special")
                   updateSliderInput(session,"downloadPlotRes", value = 100)
                   updateSliderInput(session, "downloadPlotQuality", value = 100)
                   updateSelectInput(session, "downloadPlotTypeJPEG", selected = "cairo")
                   updateSelectInput(session, "downloadPlotTypeCompression", selected = "lzw")
                 })
  ) #END observeEvent resetDownloadPlotSettings
  
  #If manhattan track is not selected then Recomb and LDSmooth track unticked.
  observeEvent(input$ShowHideTracks,({
    selectedTracks <- input$ShowHideTracks
    if(!"Manhattan" %in% selectedTracks){
      selectedTracks <- setdiff(selectedTracks, c("Recombination", "LDSmooth"))
      
      updateCheckboxGroupInput(session, "ShowHideTracks",
                               selected = selectedTracks)
    } #END if
  })) #END observeEvent ShowHideTracks
  
  #Action buttons to switch between nav bars 1.Input 2.Settings 3.Final
  observeEvent(input$goToPlotSettings,{
    updateNavbarPage(session, "navBarPageID", selected = "2.Plot Settings")
  })
  observeEvent(input$goToFinalPlot,{
    updateNavbarPage(session, "navBarPageID", selected = "3.Final Plot")
  })
  
  
})#END shinyServer



# TESTING -----------------------------------------------------------------


# datBedGraph <- reactive({
#   switch(input$dataType,
#          OncoArrayFineMapping = { annotOncoFinemapEQTL },
#          OncoArrayMeta = { annotOncoNewHits },
#          Custom = {
#            #input file check
#            validate(need(input$FileBedGraph != "", "Please upload BedGraph file"))
# 
#            inFile <- input$FileBedGraph
#            req(inFile)
#            #if(is.null(inFile)){return(NULL)} else {
#            res <- fread(inFile$datapath, header = FALSE, data.table = FALSE)
#            res <- res[,1:4]
#            colnames(res) <- c("CHR", "START", "END", "SCORE")
#            return(res)
# 
#          },
#          Example = {
#            res <- fread("Data/CustomDataExample/bedGraph.txt",
#                         header = FALSE, data.table = FALSE)
#            res <- res[,1:4]
#            colnames(res) <- c("CHR", "START", "END", "SCORE")
#            return(res)
#          })
# }) # END datBedGraph


# return(x %>%
#          transmute(BP,
#                    Recomb=RECOMB,
#                    RecombAdj=Recomb * ROIPLogMax() / 100))



# ROIdatGeneticMap <- reactive({
#   req(GeneticMap1KG)
#   GeneticMap1KG[ CHR == RegionChrN() & BP >= RegionStart() & BP <= RegionEnd()]
#   
#   # mychr <- ifelse(RegionChr() == "chrX", "chr23", RegionChr())
#   # tabixRegion <- paste0(mychr,":",
#   #                       RegionStart(), "-",
#   #                       RegionEnd())
#   # x <- tabix.read.table("Data/GeneticMap1KG/GeneticMap1KG.txt.gz", tabixRegion)
#   # colnames(x) <- c("CHR", "BP", "RECOMB")
#   # 
#   # #return
#   # x
# })


# Testing vars ------------------------------------------------------------
# output$testVars <- DT::renderDataTable({
#   data.frame(Data =
#                c("input$FilterMinPlog",
#                  "input$FilterMinLD",
#                  "RegionFlank",
#                  "RegionStart",
#                  "RegionEnd",
#                  "zoomStart", #reactive({ input$BPrange[1] })
#                  "zoomEnd",   #reactive({ input$BPrange[2] })
#                  "input$suggestiveLine, input$genomewideLine",
#                  "input$Chr",
#                  
#                  "RegionHits",
#                  "RegionHitsSelected",
#                  
#                  "datAssoc",
#                  "ROIdatAssoc",
#                  "plotDatAssoc",
#                  "colnames(plotDatAssoc)",
#                  
#                  "datLD",
#                  "ROIdatLD",
#                  "plotDatLD",
#                  "colnames(plotDatLD)",
#                  
#                  "datStats",
#                  "colnames(datStats)",
#                  
#                  "ROIdatGeneticMap",
#                  "plotDatGeneticMap",
#                  "colnames(plotDatGeneticMap)"
#                ),
#              Size = c(
#                input$FilterMinPlog,
#                input$FilterMinLD,
#                input$Flank,
#                RegionStart(),
#                RegionEnd(),
#                zoomStart(),
#                zoomEnd(),
#                paste(input$suggestiveLine, input$genomewideLine, sep = ","),
#                input$Chr,
#                
#                paste(RegionHits(), collapse = ","),
#                paste(RegionHitsSelected(), collapse = ","),
#                
#                paste(dim(datAssoc()), collapse = ","),
#                paste(dim(ROIdatAssoc()), collapse = ","),
#                paste(dim(plotDatAssoc()), collapse = ","),
#                paste(colnames(plotDatAssoc()), collapse = ","),
#                
#                paste(dim(datLD()), collapse = ","),
#                paste(dim(ROIdatLD()), collapse = ","),
#                paste(dim(plotDatLD()), collapse = ","),
#                paste(colnames(plotDatLD()), collapse = ","),
#                
#                paste(dim(datStats()), collapse = ","),
#                paste(colnames(datStats()), collapse = ","),
#                
#                paste(dim(ROIdatGeneticMap()), collapse = ","),
#                paste(dim(plotDatGeneticMap()), collapse = ","),
#                paste(colnames(plotDatGeneticMap()), collapse = ",")
#              )
#   ) %>% datatable(options = list(pageLength = 100))
#   
#   
# })

