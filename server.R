# Workspace ---------------------------------------------------------------


# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {  
  
  source("source/UDF.R",local = TRUE)
  
  # Data level 1 - Raw Input ------------------------------------------------
  datStats <- reactive({
    switch(input$dataType,
           Prostate = {},
           Custom = {
             #input file check
             validate(need(input$FileStats != "", "Please upload file"))
             
             inFile <- input$FileStats
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE) 
           },
           Example = {
             fread("Data/CustomDataExample/stats.txt",
                   header=TRUE, data.table=FALSE)
           })
  })
  
  datLD <- reactive({
    switch(input$dataType,
           Prostate = {},
           Custom = {
             #input file check
             validate(need(input$FileLD != "", "Please upload file"))
             
             inFile <- input$FileLD
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE)
           },
           Example = {
             fread("Data/CustomDataExample/LD.txt",
                   header=TRUE, data.table=FALSE) 
           })
  })
  
  datLNCAP <- reactive({
    switch(input$dataType,
           Prostate = {},
           Custom = {
             #input file check
             validate(need(input$FileLNCAP != "", "Please upload file"))
             
             inFile <- input$FileLNCAP
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE) 
           },
           Example = {
             fread("Data/CustomDataExample/LNCAP.txt",
                   header=TRUE, data.table=FALSE) 
           })
  })
  
  datEQTL <- reactive({
    switch(input$dataType,
           Prostate = {},
           Custom = {
             #input file check
             validate(need(input$FileEQTL != "", "Please upload file"))
             
             inFile <- input$FileEQTL
             if(is.null(inFile)){return(NULL)}
             fread(inFile$datapath, header=TRUE, data.table=FALSE)
           },
           Example = {
             fread("Data/CustomDataExample/eQTL.txt",
                   header=TRUE, data.table=FALSE) 
           })
  })
  
  # Define ROI --------------------------------------------------------------
  RegionChr <- reactive({ datStats()$CHR[1] })
  RegionStart <- reactive({ 
    x <- datStats() %>% 
      filter(CHR==RegionChr() &
               !is.na(BP)) %>% 
      .$BP %>% min 
    #round down to 10K bp
    as.integer(max(0,round((x-10000)/10000) * 10000))
    })
  RegionEnd <- reactive({ 
    x <- datStats() %>% 
      filter(CHR==RegionChr() &
               !is.na(BP)) %>% 
      .$BP %>% max 
    #round up to 10K bp
    as.integer(round((x+10000)/10000) * 10000)
    })
  RegionHits <- reactive({ 
    #CHR_A	BP_A	SNP_A	CHR_B	BP_B	SNP_B	R2
    #maximum LD SNPs is 5
    d <- datLD() %>% .$SNP_A %>% unique %>% sort
    d[1:min(length(d),5)] })
  
  RegionHitsSelected <- reactive({input$HitSNPs})
  
  # Genetic Map 1KG ---------------------------------------------------------
  datGeneticMap <- reactive({
    fread(paste0("Data/GeneticMap1KG/genetic_map_",
                 RegionChr(),
                 "_combined_b37.txt"),
          header=TRUE, sep=" ", data.table = FALSE)  %>% 
      transmute(BP=position,
                Recomb=`COMBINED_rate(cM/Mb)`,
                RecombAdj=Recomb*ROIPLogMax()/100)})
  
  # Data level 2 - ROI ------------------------------------------------------
  ROIdatStats <- reactive({ datStats() %>% 
      filter(BP>=RegionStart() &
               BP<=RegionEnd()) %>% 
      mutate(PLog=-log10(P)) })
  ROIdatLD <- reactive({ datLD() %>% 
      filter(CHR_B==as.numeric(gsub("chr","",RegionChr())) &
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
  ROIdatGeneticMap <- reactive({ datGeneticMap() %>% 
      filter(BP>=RegionStart() &
               BP<=RegionEnd()) })
  ROIPLogMax <- reactive({
    #ylim Max for plotting Manhattan
    maxY <- max(ROIdatStats()$PLog,na.rm = TRUE)
    return(max(10,ceiling((maxY+1)/5)*5))
  })
  
  # Define Zoom Start End ---------------------------------------------------
  zoomStart <- reactive({ input$BPrange[1] })
  zoomEnd <- reactive({ input$BPrange[2] })
  
  
  # Data level 3 - Plot data ------------------------------------------------
  # http://stackoverflow.com/a/8197703/680068
  # hcl(h=seq(15, 375, length=6), l=65, c=100)[1:5]
  # colors match first 5 colours of default ggplot2
  colourLD5 <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
  #create pallete LD 0 to 100
  colLD <- lapply(colourLD5,function(i)colorRampPalette(c("grey95",i))(100))
  
  plotDatManhattan <- reactive({
    Stats <- 
      ROIdatStats() %>% 
      filter(BP>=zoomStart() &
               BP<=zoomEnd() &
               PLog >= input$FilterMinPlog)
    
    
    LD <- ROIdatLD() %>% 
      filter(
        R2 >= input$FilterMinLD &
          SNP_A %in% RegionHitsSelected())
    
    #assign colours to LD 
    d_LD <- base::do.call(rbind,
                          lapply(RegionHitsSelected(), function(snp){
                            #snp="rs13410475"
                            d <- LD %>% filter(SNP_A == snp)
                            #LD round minimum is 1 max 100
                            LDColIndex <- ifelse(round(d$R2,2)==0,1,round(d$R2,2)*100)
                            LDColIndex <- ifelse(LDColIndex>100,100,LDColIndex)
                            d$LDSNP <- snp
                            d$LDSmoothCol <- colLD[[match(snp,RegionHits())]][100]
                            d$LDCol <- colLD[[match(snp,RegionHits())]][LDColIndex]
                            return(d)
                          }))
    # to add smooth LD Y value used for Pvalue and LD 0-1
    # add pvalues for Y value on the plot
    d_LD <- 
      base::merge.data.frame(
        Stats,
        d_LD[,c("BP_B","R2","LDSNP","LDSmoothCol","LDCol")],
        by.x="BP",by.y="BP_B",all.x=TRUE) %>% 
      mutate(R2_Adj=ROIPLogMax()*R2)
    
    
    return(d_LD)
  })
  
  #subset of recombination rates for zoomed region
  plotDatGeneticMap <- reactive({ datGeneticMap() %>% 
      filter(BP>=zoomStart() &
               BP<=zoomEnd()) })
  
  #get granges collapsed genes for ggplot+ggbio
  plotDatGene <- reactive({ 
    udf_GeneSymbol(chrom=RegionChr(),
                   chromStart=zoomStart(),
                   chromEnd=zoomEnd()) })
  
  #number of genes in zoomed region
  plotDatGeneN <- reactive({ 
    length(unique(plotDatGene()@elementMetadata$gene_id)) })
  
  
  # Output ------------------------------------------------------------------
  # Output Summary ----------------------------------------------------------
  output$SummaryStats <- renderDataTable({ datStats() %>% arrange(BP) }) 
  output$SummaryLD <- renderDataTable({ datLD() %>% group_by(SNP_A) %>% arrange(BP_B) })
  output$SummaryLNCAP <- renderDataTable({ datLNCAP() %>% arrange(BP) })
  output$SummaryEQTL <- renderDataTable({ datEQTL() %>% arrange(START) })
  
  output$SummaryRegion <- renderTable({ 
    data.frame(Chr=RegionChr(),
               Start=RegionStart(),
               End=RegionEnd(),
               Size=RegionEnd()-RegionStart(),
               Hits=paste(RegionHits(),collapse=", "))})
  
  output$SummaryFileNrowNcol <- renderTable({ 
    data.frame(InputFile=c("Association","LD","LNCAP","eQTL"),
               Nrow=c(nrow(datStats()),nrow(datLD()),
                      nrow(datLNCAP()),nrow(datEQTL())),
               Ncol=c(ncol(datStats()),ncol(datLD()),
                      ncol(datLNCAP()),ncol(datEQTL()))
    )})
  
  output$SummaryHeadPlotStats <- renderTable({head(plotDatStats())})
  output$SummaryHeadPlotDatManhattan <- renderDataTable({plotDatManhattan()})
  output$SummaryDimPlotDatManhattan <- renderTable({as.data.frame(dim(plotDatManhattan()))})
  #output$SummaryZoom <- renderText({paste(zoomStart(),zoomEnd(),sep="-")})
  
  # Plot --------------------------------------------------------------------
  #Plot Chr ideogram
  plotObjChromosome <- reactive({source("source/Chromosome.R",local=TRUE)})
  output$PlotChromosome <- renderPlot({print(plotObjChromosome())})
  #Manhattan track
  plotObjManhattan <- reactive({source("source/Manhattan.R",local=TRUE)})
  output$PlotManhattan <- renderPlot({print(plotObjManhattan())})
  #SNPType track
  plotObjSNPType <- reactive({source("source/SNPType.R",local=TRUE)})
  output$PlotSNPType <- renderPlot({print(plotObjSNPType())})
  #SNP LD track
  plotObjSNPLD <- reactive({source("source/LD.R",local=TRUE)})
  output$PlotSNPLD <- renderPlot({print(plotObjSNPLD())})
  #LNCAP Smooth track
  plotObjLNCAP <- reactive({source("source/LNCAP.R",local=TRUE)})
  output$PlotLNCAP <- renderPlot({print(plotObjLNCAP())})
  #eQTL bar track
  plotObjEQTL <- reactive({source("source/eQTL.R",local=TRUE)})
  output$PlotEQTL <- renderPlot({print(plotObjEQTL())})
  #Gene track
  plotObjGene <- reactive({source("source/Genes.R",local=TRUE)})
  output$PlotGene <- renderPlot({print(plotObjGene())})
  
  
  # Dynamic UI --------------------------------------------------------------
  #Zoom to region X axis BP
  output$BPrange <-
    renderUI({
      sliderInput("BPrange", h5("Region: Start-End"),
                  min = RegionStart(),
                  max = RegionEnd(),
                  value = c(RegionStart(),RegionEnd()),
                  step = 50000)})
  
  #Select hit SNPs
  output$HitSNPs <- 
    renderUI({
      checkboxGroupInput("HitSNPs", "Hit SNPs",
                         RegionHits(), selected = RegionHits())
    })
  observe({
    # maximum of 5 SNPs can be selected to display LD, minimum 1 must be ticked.
    if(length(input$HitSNPs) > 5){
      updateCheckboxGroupInput(session, "HitSNPs", selected= tail(input$HitSNPs,5))}
    if(length(input$HitSNPs) < 1){
      updateCheckboxGroupInput(session, "HitSNPs", selected= RegionHits()[1])}
  })
  
  
  
})#END shinyServer

