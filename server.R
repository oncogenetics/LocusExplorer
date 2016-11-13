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
             
#              TESTING validation...
#                customData <- fread(inFile$datapath, header=TRUE, data.table=FALSE) 
#                #check if the input has required columns
#                if(isTRUE(all.equal(colnames(customData),
#                                    c("CHR","SNP","BP","P","TYPED")))
#                   ){return(NULL)} else {
#                     return(customData) }

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
  
#   # datLNCAP <- reactive({
#   #   fread("Data/ProstateLNCAP/ProstateLNCAP.txt",
#   #         colClasses = c("character","numeric","character")) })
#   
#   datAnnotSmooth <- reactive ({
#     annotSmooth <- fread("Data/Annotation/annotSmooth.txt")
#     annotSmooth$TYPEN <-
#       as.numeric(
#         factor(annotSmooth$TYPE,
#                levels = c(
#                  "DNaseI","Conserved",
#                  #ChromHMM
#                  "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
#                  "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter")))
#     return(annotSmooth)
#   })

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

#   datAnnotEQTL <- reactive({
#     # Annotation data from Ed and Eze, merged, multiple eQTL per SNPs (many to one)
#     d <- fread("Data/Annotation/paper_regions_All_SNPs_annotated_Eze_4_Tokhir.tsv")
#     d$Gene <- 
#       sapply(d$eqtl.Eze.gene,
#              function(i)unlist(strsplit(i, split = "|", fixed = TRUE))[1])
#     return(d)
#   })
#     
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
 
#   datHaploreg <-  reactive({
#      
#       datHaploreg <- fread("Data/Haploreg/20160323.txt", sep = "\t", data.table = FALSE)
#       datHaploreg[datHaploreg=="." | datHaploreg==""] <- NA
#       #x <- RegionChrN()
#       datHaploreg %>% 
#         transmute(
#           CHR=chr,
#           SNP=rsID,
#           SiPhy_cons,
#           EnhancerHistoneMarks=as.integer(!is.na(Chromatin_Marks)),
#           DNAse=as.integer(!is.na(DNAse)),
#           ProteinBound=as.integer(!is.na(Proteins)),
#           MotifsChanged=as.integer(!is.na(Motifs)),
#           GERP_cons=GERP_cons,
#           eQTL=as.integer(is.na(grasp)+is.na(eQTL) < 2),
#           functional=substr(dbSNP_functional_annotation,1,3)) %>% 
#         gather(annot,value,SiPhy_cons:functional) %>% 
#         filter(value!=0) %>% 
#         mutate(annot=ifelse(value==1,as.character(annot),paste0("Func_",value))) %>% 
#         dplyr::select(CHR,SNP,annot) %>% 
#         filter(annot %in% input$haploreg &
#                  CHR == RegionChrN())
#   
#       })
#   
#       
  # Define ROI --------------------------------------------------------------
  RegionFlank <- reactive({
    #if(input$Flank == "0"){rf <- 1} else {rf <- input$Flank}
    #as.numeric(rf)
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
    #as.integer(max(0,round((x-10000)/10000) * 10000))
    })
  RegionEnd <- reactive({
    x <- datAssoc() %>%
      filter(CHR == RegionChr() & !is.na(BP)) %>%
      .$BP %>% max
    #round up to 10K bp
    as.integer(round((x + RegionFlank())/RegionFlank()) * RegionFlank())
    #as.integer(round((x+10000)/10000) * 10000)
    })

  RegionHits <- reactive({
    d <- datLD() %>% .$SNP_A %>% unique %>% sort
    #d[1:min(length(d),15)]
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
  # ROIdatLNCAP <- reactive({ datLNCAP() %>%
  #     filter(CHR==RegionChr() &
  #              BP>=RegionStart() &
  #              BP<=RegionEnd()) })
  # ROIdatAnnotSmooth <- reactive({ datAnnotSmooth() %>%
  #     filter(CHR == RegionChr() &
  #              START >= RegionStart() &
  #              END <= RegionEnd()) })
  ROIdatBedGraph <- reactive({
    datBedGraph() %>%
      filter(CHR==RegionChr()
             # &
             #   START>=RegionStart() &
             #   END<=RegionEnd()
             ) %>%
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

#   ROIdatAnnotEQTL <- reactive({
#     datAnnotEQTL() %>% 
#       dplyr::filter(Chr == RegionChrN() &
#                Start >= RegionStart() &
#                End <= RegionEnd() &
#                MargBF > 100 &
#                selected == 1) %>% 
#       dplyr::transmute(SNP = name,
#                 GENE = Gene,
#                 x = Start,
#                 xend = gene.start,
#                 eQTL_TCGA = 1,
#                 yend = 0.35) %>% 
#       base::unique()
#   })
#   
#   ROIdatAnnotEQTL_GeneLabel <- reactive({
#     # data is used to add Xaxis labels - Gene names, also to add
#     #  vertical lines on gene plot to mark genes with a line.
#     ROIdatAnnotEQTL() %>% 
#       dplyr::transmute(x = xend, eQTL_TCGA = 0.35, label = GENE) %>%
#       base::unique()
#     })
#     
#   ROIdatSankeyAnnotation <- reactive({
#     #step1 - method to hits
#     dd1 <- datHitSNPStats()
#     step1 <-
#       rbind(
#         dd1 %>%
#           dplyr::filter(Method == "Stepwise Forward") %>%
#           dplyr::transmute(From = "Stepwise",
#                     To = SNP,
#                     W = round(rescale(-log10(Stats), to = c(10,100)))
#                     #W=10
#                     ),
#         dd1 %>%
#           dplyr::filter(Method == "BVS" &
#                           Stats > as.numeric(input$FilterMinBVS_BF)) %>%
#           dplyr::transmute(From = "BVS",
#                     To = SNP,
#                     W = round(rescale(Stats, to = c(10,100)))
#                     #W=10
#                     )
#       )#rbind
#     
#     #step2 - hits to Annotation
#     dd2 <- datAnnotEQTL()
#     step2 <-
#       dd2 %>%
#       dplyr::filter(
#         Chr == RegionChrN() &
#           Start >= RegionStart() &
#           End <= RegionEnd() &
#           best_tag %in% dd1$rsid & 
#           MargBF > as.numeric(input$FilterMinBVS_BF)) %>%
#                       #MargBF > 100) %>%
#       dplyr::select(c(11, 24, 17:23, 25:37)) %>%
#       base::unique() %>%
#       gather(Type, score, -c(1:2)) %>%
#       #dplyr::filter(best_tag %in% step1$To & score != 0) %>%
#       dplyr::filter(score != 0) %>%
#       dplyr::transmute(From = best_tag,
#                        To = Type,
#                        W = round(rescale(RawScore, to = c(10,100))))
#     
#     # update SNP names to match Stats
#     step2 <- merge(step2, dd1[, c("rsid", "SNP")],
#                    by.x = "From", by.y = "rsid" ) %>% 
#       mutate(From = SNP) %>% 
#       dplyr::select(From, To, W) %>% 
#       base:: unique()
#     
#     step3 <- 
#       step1[ ! step1$To %in% step2$From, ] %>% 
#       dplyr::transmute(From = To,
#                        To = "None",
#                        W = 10)
#     
#     #return: if W is inf then max
#     res <- base:: unique(rbindlist(list(step1, step2, step3)))
#     #res <- step2
#     res$W <- ifelse(is.infinite(res$W), 100, res$W)
#     return(res)
#   })
# 
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

  #output$tempSummaryplotDatLD <- renderDataTable(plotDatLD())

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


#   #subset eQTL TCGA for zoom
#   plotDatAnnotEQTL <- reactive({
#     ROIdatAnnotEQTL() %>% 
#       filter(x >= zoomStart() &
#                x <= zoomEnd() &
#                xend >= zoomStart() &
#                xend <= zoomEnd()
#                )
#     })
#   plotDatAnnotEQTL_GeneLabel <- reactive({
#     ROIdatAnnotEQTL_GeneLabel() %>% 
#       filter(x >= zoomStart() &
#                x <= zoomEnd()) 
#     })
#   
#   
#   #wgEncodeBroadHistone data 7 big wig data
#   plotDatwgEncodeBroadHistone <- reactive({
#     #get file list matching Encode description file
#     bigWigFiles <- 
#       intersect(wgEncodeBroadHistoneFileDesc$File,
#                 list.files("Data/wgEncodeBroadHistone/","*.bigWig"))
#     
#     #check if bigWig files are downloaded, if plot warning message.
#     if(length(bigWigFiles)>0){
#       #subset region 
#       wgEncodeBroadHistone <- 
#         rbind_all(
#           #seven bigwig files
#           lapply(bigWigFiles, function(i){
#             #i=1  wgEncodeBroadHistoneFileDesc$File[1]
#             as.data.frame(
#               import(paste0("Data/wgEncodeBroadHistone/",i),
#                      which=RegionGR())) %>% 
#               filter(score >= 5) %>% 
#               transmute(BP=start,
#                         ENCODE=wgEncodeBroadHistoneFileDesc[
#                           which(wgEncodeBroadHistoneFileDesc$File==i),
#                           "Name"],
#                         SCORE=round(ifelse(score >= 100, 100, score),0))
#           }))
#       
#       #merge on BP, to plot overlappling
#       res <- data.frame(BP=unique(wgEncodeBroadHistone$BP))
#       output <- 
#         rbind_all(
#           lapply(unique(wgEncodeBroadHistone$ENCODE), function(i){
#             d <- merge(res,wgEncodeBroadHistone %>% filter(ENCODE==i),by="BP",all.x=TRUE)
#             d$ENCODE <- i
#             d$SCORE[ is.na(d$SCORE) ] <- 0
#             return(d)
#           })
#         )
#       }else{
#       #if there are no bigwig files return blank data.frame
#       output <- data.frame()
#     }
#     return(output)
#     })
#   
#   #output of ENCODE data
#   #wgEncodeBroadHistone data
#   output$SummarywgEncodeBroadHistone <- 
#     DT::renderDataTable(plotDatwgEncodeBroadHistone() %>% arrange(BP),
#                         options=list(searching=FALSE,searchable=FALSE))
#   #wgEncodeRegDnaseClustered
#   output$SummarywgEncodeRegDnaseClustered <-
#     DT::renderDataTable(ROIdatwgEncodeRegDnaseClustered() %>% arrange(START),
#                         options=list(searching=FALSE,searchable=FALSE))
#   
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
#   
#   
#   # output$SummaryLNCAP <- 
#   #   DT::renderDataTable(ROIdatLNCAP() %>% arrange(BP),
#   #                       options=list(searching=FALSE,searchable=FALSE))
#   
#   output$SummaryAnnotSmooth <- 
#     DT::renderDataTable(ROIdatAnnotSmooth() %>% arrange(START),
#                         options=list(searching=FALSE,searchable=FALSE))
#   
#   output$SummaryBedGraph <- 
#     DT::renderDataTable(datBedGraph() %>% arrange(START),
#                         options=list(searching=FALSE,searchable=FALSE))
#   
#   output$SummaryROIdatAnnotEQTL <- 
#     DT::renderDataTable(ROIdatAnnotEQTL() %>% arrange(x),
#                         options=list(searching=FALSE,searchable=FALSE))

  output$SummaryRegion <-
    renderUI(a(paste0(RegionChr(),':',RegionStart(),'-',RegionEnd()),
               href=paste0('http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',
                           RegionChr(),'%3A',RegionStart(),'-',RegionEnd()),target="_blank"))

#   output$SummaryHits <- DT::renderDataTable( 
#     data.table(SNP = RegionHits()) %>% 
#       #if SNP name has rs number then convert to a link to NCBI
#       mutate(SNP=ifelse(substr(SNP,1,2)=="rs",
#                         paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
#                                gsub("rs","",SNP),'" target="_blank">',SNP,'</a>'),
#                         SNP)),
#     escape=FALSE,options=list(paging=FALSE,searching=FALSE,searchable=FALSE))
# 
# #   output$SummaryFileNrowNcol <- DT::renderDataTable( 
# #     data.table(InputFile=c("Association","LD","LNCAP","BED"),
# #                RowCount=c(nrow(datAssoc()),nrow(datLD()),
# #                       nrow(datLNCAP()),nrow(datBED())),
# #                ColumnCount=c(ncol(datAssoc()),ncol(datLD()),
# #                       ncol(datLNCAP()),ncol(datBED()))),
# #     options=list(paging=FALSE,searching=FALSE,searchable=FALSE))
#   
#   #output$tempPlotDatLD <- DT::renderDataTable(ROIdatLD())
#   
#   
#   output$SummaryLDlink <- DT::renderDataTable(datLDlink())
#   output$SummaryLDlinkProcess <- 
#     DT::renderDataTable({
#       datLDlinkProcess() %>% 
#         arrange(BP_B) %>% 
#         mutate(SNP_A=ifelse(substr(SNP_A,1,2)=="rs",
#                             paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
#                                    gsub("rs","",SNP_A),'" target="_blank">',SNP_A,'</a>'),
#                             SNP_A),
#                SNP_B=ifelse(substr(SNP_B,1,2)=="rs",
#                             paste0('<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
#                                    gsub("rs","",SNP_B),'" target="_blank">',SNP_B,'</a>'),
#                             SNP_B))
#     },escape=FALSE)
#   
#   output$SummaryHaploreg <- DT::renderDataTable(datHaploreg())
#   
#   
#   output$SummaryROIdatSankeyAnnotation <- DT::renderDataTable(ROIdatSankeyAnnotation())
#   
#   output$SummarydatAnnotEQTL <- DT::renderDataTable(datAnnotEQTL() %>%
#                                                       filter(best_tag %in% datHitSNPStats()$rsid), 
#                                                     options = list(pageLength = 100))
#       
#   # Plot --------------------------------------------------------------------
#   #plot title
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

    source("Source/Manhattan.R",local=TRUE)})
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


#   #LNCAP Smooth track
#   #plotObjLNCAP <- reactive({source("Source/LNCAP.R",local=TRUE)})
#   #output$PlotLNCAP <- renderPlot({print(plotObjLNCAP())})
#   
#   #annotSmooth track
#   plotObjAnnotSmooth <- reactive({source("Source/annotSmooth.R",local = TRUE)})
#   output$PlotAnnotSmooth <- renderPlot({print(plotObjAnnotSmooth())})
#   
#   #BedGraph bar track
#   plotObjBedGraph <- reactive({source("Source/BedGraph.R",local=TRUE)})
#   output$PlotBedGraph <- renderPlot({print(plotObjBedGraph())})
#   
#   #Gene track
#   plotObjGene <- reactive({
#     # Create a Progress object
#     progress <- shiny::Progress$new()
#     progress$set(message = "Plotting please wait...", value = 0)
#     # Close the progress when this reactive exits (even if there's an error)
#     on.exit(progress$close())
#     
#     source("Source/Gene.R",local=TRUE)})
#   output$PlotGene <- renderPlot({
#     # Create a Progress object
#     progress <- shiny::Progress$new()
#     progress$set(message = "Almost there...", value = 80)
#     # Close the progress when this reactive exits (even if there's an error)
#     on.exit(progress$close())
#     
#     print(plotObjGene())})
#   
#   
#   # ROIdatSankeyAnnotation
#   output$PlotSankeyAnnotation <- renderGvis({
#     dat <- ROIdatSankeyAnnotation()
#     gvisSankey(dat, from = "From",
#                to = "To", weight = "W",
#                options=list(
#                  height = 800,
#                  width = 800,
#                  #sankey="{link:{colorMode: 'gradient'}}"))
#                  sankey="{
#                  node:{nodePadding: 20,
#                  colors: [
#                  '#C9DD03',
#                  '#FFD602',
#                  '#F9A100',
#                  '#EE7EA6',
#                  '#A71930',
#                  '#616365',
#                  '#726E20',
#                  '#6E273D',
#                  '#F00034',
#                  '#ADAFAF',
#                  '#003D4C'
#                  ]},
#                  link:{colorMode: 'source',
#                  colors: [
#                  '#C9DD03',
#                  '#FFD602',
#                  '#F9A100',
#                  '#EE7EA6',
#                  '#A71930',
#                  '#616365',
#                  '#726E20',
#                  '#6E273D',
#                  '#F00034',
#                  '#ADAFAF',
#                  '#003D4C'
#                  ]
#                  }
#   }"))
# 
#       
#       }) #END: PlotSankeyAnnotation
#   
#   
#   
#   
#   # Example plot ------------------------------------------------------------
#   output$ExamplePlotJPEG <- renderImage({
#     return(list(
#         src = "www/Figure1.jpe",
#         contentType = "image/jpeg",
#         alt = "ExamplePlotOutput"
#       ))
#     }, deleteFile = FALSE)
#   
#   # Graph plots for Oncoarray -----------------------------------------------
#   # graph plot of hits and methods
#   output$plotGraphStats <- renderForceNetwork({
#     # Data
#     datGraphStats <- datHitSNPStats()
#     datGraphAssoc <- datAssoc()
#     
#     MisNodes <- data.table(name=unique(c(datGraphStats$Method,datGraphStats$SNP)),
#                            id=0:(length(unique(c(datGraphStats$Method,datGraphStats$SNP)))-1))
#     
#     MisLinks <- merge(datGraphStats,MisNodes,by.x="Method",by.y="name",all.x=TRUE)
#     MisLinks <- merge(MisLinks,MisNodes,by.x="SNP",by.y="name",all.x=TRUE)
#     
#     
#     MisLinks <- MisLinks %>% 
#       mutate(source=id.x,
#              target=id.y,
#              value=1) %>% 
#       arrange(source,target)
#     
#     MisNodes$group <- 
#       ifelse(MisNodes$name %in% c("BVS","ElasticNet","Stepwise Forward"),
#              "Method",
#              ifelse(MisNodes$name %in% 
#                       datHaploreg()[datHaploreg()$annot %in% input$haploreg,"SNP"],
#                     "Haploreg",
#                     "Other"))
#     
#     #get Pvalues, use to set the size of nodes
#     MisNodes <- merge(MisNodes,datGraphAssoc[,c("SNP","P")],
#                       by.x="name",by.y="SNP",all.x=TRUE)
#     
#     MisNodes$size <- ifelse(is.na(MisNodes$P),1,
#                             ifelse(MisNodes$P<0.0001,50,1))
#     MisNodes <- MisNodes %>% arrange(id)
#     
#     #plot
#     forceNetwork(Links=MisLinks, Nodes=MisNodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  Group="group",opacity=0.9,zoom=TRUE,legend=TRUE,
#                  bounded=TRUE,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()"))  
#     
#   }) # END plotGraphStats
#   
#   
#   # graph plot of hits and other SNPs that are connected to at least 2 hits.
#   output$plotGraphLD_All <- renderForceNetwork({
#     # Data
#     datGraph <- datLD() %>% 
#       filter(SNP_A!=SNP_B & 
#                R2 >= input$FilterMinLD) %>% 
#       group_by(SNP_B) %>% 
#       filter(n()>1)
#     
#     datGraphAssoc <- datAssoc()
#     
#     # Nodes and Links
#     MisNodes <- data.table(name=unique(c(datGraph$SNP_A,datGraph$SNP_B)),
#                            #id must be zero based
#                            id=0:(length(unique(c(datGraph$SNP_A,datGraph$SNP_B)))-1))
#     MisLinks <- merge(datGraph,MisNodes,by.x="SNP_A",by.y="name",all.x=TRUE)
#     MisLinks <- merge(MisLinks,MisNodes,by.x="SNP_B",by.y="name",all.x=TRUE)
#     MisLinks <- MisLinks %>% 
#       transmute(source=ifelse(id.x<=id.y,id.x,id.y),
#                 target=ifelse(id.y>id.x,id.y,id.x),
#                 value=R2) %>% unique
#     #group nodes
#     MisNodes <- MisNodes %>% 
#       mutate(group=ifelse(name %in% datHitSNPStats()$SNP,"Hits","Other")) 
#     
#     #get Pvalues, use to set the size of nodes
#     MisNodes <- merge(MisNodes,datGraphAssoc[,c("SNP","P")],by.x="name",by.y="SNP",all.x=TRUE)
#     MisNodes$size <- ifelse(is.na(MisNodes$P),1,
#                             ifelse(MisNodes$P<0.0001,50,1))
#     MisNodes <- MisNodes %>% arrange(id)
# 
#     #plot
#     forceNetwork(Links=MisLinks, Nodes=MisNodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()"))
#     
#     }) #END plotGraphLD_All
#   
#   
#   # graph plot of for Stepwise hits
#   output$plotGraphLD_StepwiseForward <- renderForceNetwork({
#     res <- udf_NodesAndLinks("Stepwise Forward")
#     
#     #plot
#     forceNetwork(Links=res$MisLinks, Nodes=res$MisNodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  radiusCalculation = JS("d.nodesize"),
#                  Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
#                  charge=-300,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()"))
#   }) #END plotGraphLD_StepwiseForward
#   
# #   output$plotGraphLD_StepwiseBackward <- renderForceNetwork({
# #     res <- udf_NodesAndLinks("Stepwise Backward")
# #     #plot
# #     forceNetwork(Links=res$MisLinks, Nodes=res$MisNodes,
# #                  Source="source", Target="target",
# #                  Value="value",NodeID="name",Nodesize="size",
# #                  radiusCalculation = JS("d.nodesize"),
# #                  Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
# #                  charge=-300,
# #                  opacityNoHover = 1,
# #                  colourScale = JS("d3.scale.category10()"))
# #   }) #END plotGraphLD_StepwiseBackward
#   
#   output$plotGraphLD_ElasticNet <- renderForceNetwork({
#     res <- udf_NodesAndLinks("ElasticNet")
#     #plot
#     forceNetwork(Links=res$MisLinks, Nodes=res$MisNodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  radiusCalculation = JS("d.nodesize"),
#                  Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
#                  charge=-300,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()"))
#   }) #END plotGraphLD_ElasticNet
#   
#   output$plotGraphLD_BVS <- renderForceNetwork({
#     res <- udf_NodesAndLinks("BVS")
#     #plot
#     forceNetwork(Links=res$MisLinks, Nodes=res$MisNodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  radiusCalculation = JS("d.nodesize"),
#                  Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
#                  charge=-300,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()"))
#   }) #END plotGraphLD_BVS
#   
#   
#   # graph plot of for Stepwise hits
#   # output$plotGraphLD_ElasticNet <- renderForceNetwork({
#   #   # Data
#   #   #c("BVS","ElasticNet","Stepwise Backward","Stepwise Forward"),"Method","SNP")   
#   #   
#   #   #ElasticNet hits only
#   #   methodHits <- datHitSNPStats() %>% 
#   #     filter(Method == "ElasticNet") %>% 
#   #     .$SNP %>% unique
#   #   
#   #   #make all combos for hits, in case there is no LD, then we set it to 0
#   #   # so they show up on graph plot even if they are not connected to anything.
#   #   methodHitsCombn <-
#   #     as.data.frame(t(combn(methodHits,2)))
#   #   colnames(methodHitsCombn) <- c("SNP_A","SNP_B")
#   #   methodHitsCombn$SNP_A <- as.character(methodHitsCombn$SNP_A)
#   #   methodHitsCombn$SNP_B <- as.character(methodHitsCombn$SNP_B)
#   #   
#   #   datGraph <- merge(
#   #     datLD()[,c("SNP_A","SNP_B","R2")] %>% 
#   #       filter(SNP_A!=SNP_B & 
#   #                R2 >= input$FilterMinLD &
#   #                (SNP_A %in% methodHits |
#   #                   SNP_B %in% methodHits)) %>% 
#   #       group_by(SNP_A) %>% 
#   #       filter(n()>1),
#   #     methodHitsCombn, 
#   #     all=TRUE)
#   #   #set missing R2 to zero - 0
#   #   datGraph$R2 <- ifelse(is.na(datGraph$R2),0,datGraph$R2)
#   #   
#   #   #get single SNP assoc results to set size of circles
#   #   datGraphAssoc <- datAssoc()
#   #   
#   #   #prepare graph plot data
#   #   # Nodes and Links
#   #   MisNodes <- data.table(name=unique(c(datGraph$SNP_A,datGraph$SNP_B)),
#   #                          #id must be zero based
#   #                          id=0:(length(unique(c(datGraph$SNP_A,datGraph$SNP_B)))-1))
#   #   MisLinks <- merge(datGraph,MisNodes,by.x="SNP_A",by.y="name",all.x=TRUE)
#   #   MisLinks <- merge(MisLinks,MisNodes,by.x="SNP_B",by.y="name",all.x=TRUE)
#   #   MisLinks <- MisLinks %>% 
#   #     transmute(source=ifelse(id.x<=id.y,id.x,id.y),
#   #               target=ifelse(id.y>id.x,id.y,id.x),
#   #               value=R2) %>% unique
#   #   #group nodes
#   #   MisNodes <- MisNodes %>% 
#   #     mutate(group=ifelse(name %in% methodHits,"Hits","Other")) 
#   #   
#   #   #get Pvalues, use to set the size of nodes
#   #   MisNodes <- merge(MisNodes,datGraphAssoc[,c("SNP","P")],
#   #                     by.x="name",by.y="SNP",all.x=TRUE)
#   #   MisNodes$size <- ifelse(is.na(MisNodes$P),1,
#   #                           ifelse(MisNodes$P<0.0001,50,1))
#   #   MisNodes <- MisNodes %>% arrange(id)
#   #   
#   #   #plot
#   #   forceNetwork(Links=MisLinks, Nodes=MisNodes,
#   #                Source="source", Target="target",
#   #                Value="value",NodeID="name",Nodesize="size",
#   #                Group="group",opacity=0.8,zoom=TRUE,legend=TRUE,
#   #                charge=-300,
#   #                opacityNoHover = 1,
#   #                colourScale = JS("d3.scale.category10()"))
#   #   
#   # }) #END plotGraphLD_Stepwise
#   
#   
#   # graph plot of methods, hits, and other SNPs
#   output$plotGraphLDStats <- renderForceNetwork({
#     datGraphStats <- datHitSNPStats()
# #     datGraphStats <- fread("chr2_241657087_242920971_stats.txt",
# #                            header=TRUE, data.table=FALSE)
#     
#     datGraphAssoc <- datAssoc()
# #     datGraphAssoc <- fread("chr2_241657087_242920971_assoc.txt",
# #                            header=TRUE, data.table=FALSE)
#     
#     #datGraphLD <- datLD()
# #     datGraphLD <- fread("chr2_241657087_242920971_LD.txt",
# #                         header=TRUE, data.table=FALSE)
#     
#     datGraphLD <- datLD() %>% 
#       filter(R2 > input$FilterMinLD &
#                SNP_A!=SNP_B) %>% 
#       group_by(SNP_B) %>% 
#       filter(n()>1)
#     
#     x <- unique(c(datGraphStats$Method,
#                   datGraphStats$SNP,
#                   datGraphLD$SNP_A,
#                   datGraphLD$SNP_B))
#     Nodes <- data.table(name=x,
#                         id=0:(length(x)-1))
#     
#     Nodes$group <- 
#       ifelse(Nodes$name %in% c("BVS","ElasticNet","Stepwise Forward"),
#              "Method",
#              ifelse(Nodes$name %in% datGraphStats$SNP, "Hit SNPs","Other SNPs")) 
#     
#     #     #get Pvalues, use to set the size of nodes
#     #     Nodes <- merge(Nodes,datGraphAssoc[,c("SNP","P")],
#     #                    by.x="name",by.y="SNP",all.x=TRUE)
#     Nodes$size <- 1
#     #     Nodes$size <- ifelse(is.na(Nodes$P),1,
#     #                          ifelse(Nodes$P<0.0001,50,1))
#     #     Nodes <- Nodes %>% arrange(id)
#     
#     #links for stats groups
#     LinksStats <- merge(datGraphStats,Nodes,by.x="Method",by.y="name",all.x=TRUE)
#     LinksStats <- merge(LinksStats,Nodes,by.x="SNP",by.y="name",all.x=TRUE)
#     
#     LinksStats <- LinksStats %>% 
#       mutate(source=id.x,
#              target=id.y,
#              value=1) %>% 
#       arrange(source,target)
#     
#     
#     #links for LD between SNPs
#     LinksLD <- merge(datGraphLD,Nodes,by.x="SNP_A",by.y="name",all.x=TRUE)
#     LinksLD <- merge(LinksLD,Nodes,by.x="SNP_B",by.y="name",all.x=TRUE)
#     
#     LinksLD <- LinksLD %>% 
#       transmute(source=ifelse(id.x<=id.y,id.x,id.y),
#                 target=ifelse(id.y>id.x,id.y,id.x),
#                 value=ceiling(R2*10)) %>% 
#       arrange(source,target) %>% unique
#     
#     #Merge Methods and LD links
#     Links <- rbind(LinksStats[,c("source","target","value")],
#                    LinksLD[,c("source","target","value")]) %>% 
#       arrange(source,target)
#     
#     #plot
#     forceNetwork(Links=Links, Nodes=Nodes,
#                  Source="source", Target="target",
#                  Value="value",NodeID="name",Nodesize="size",
#                  Group="group",opacity=0.9,zoom=TRUE,legend=TRUE,
#                  opacityNoHover = 1,
#                  colourScale = JS("d3.scale.category10()")) 
#   }) #END plotGraphLDStats
#   
#   # Arc plots -----------------------------------------------------------------
#   output$plotArcLD_All <- renderPlot({
#     res <- udf_arcPlotData(method=c("Stepwise Forward","ElasticNet","BVS"))
#     #plot
#     arcplot(res$lab, lwd.arcs = 8*E(res$glab)$weight,
#             col.nodes = res$cols, bg.nodes = res$cols, show.nodes = TRUE,
#             cex.labels = 1, horizontal = TRUE)
#     }) # END plotArcLD_All
#   
#   #plot Arc LD for hit SNPs
#   output$plotArcLD_StepwiseForward <- renderPlot({
#     res <- udf_arcPlotData(method="Stepwise Forward")
#     #plot
#     arcplot(res$lab, lwd.arcs = 8*E(res$glab)$weight,
#             col.nodes = res$cols, bg.nodes = res$cols, show.nodes = TRUE,
#             cex.labels = 1, horizontal = TRUE)
#     }) # END plotArcLD_StepwiseForward
#   
# #   output$plotArcLD_StepwiseBackward <- renderPlot({
# #     res <- udf_arcPlotData(method="Stepwise Backward")
# #     #plot
# #     arcplot(res$lab, lwd.arcs = 8*E(res$glab)$weight,
# #             col.nodes = res$cols, bg.nodes = res$cols, show.nodes = TRUE,
# #             cex.labels = 1, horizontal = TRUE)
# #     }) # END plotArcLD_StepwiseForward
#   
#   output$plotArcLD_ElasticNet <- renderPlot({
#     res <- udf_arcPlotData(method="ElasticNet")
#     #plot
#     arcplot(res$lab, lwd.arcs = 8*E(res$glab)$weight,
#             col.nodes = res$cols, bg.nodes = res$cols, show.nodes = TRUE,
#             cex.labels = 1, horizontal = TRUE)
#   }) # END plotArcLD_ElasticNet
#   
#   output$plotArcLD_BVS <- renderPlot({
#     res <- udf_arcPlotData(method="BVS")
#     #plot
#     arcplot(res$lab, lwd.arcs = 8*E(res$glab)$weight,
#             col.nodes = res$cols, bg.nodes = res$cols, show.nodes = TRUE,
#             cex.labels = 1, horizontal = TRUE)
#   }) # END plotArcLD_BVS
#   
#   # Manhattan Plot Merge ------------------------------------------------------
#   #Dynamic size for tracks
#   RegionHitsCount <- reactive({ length(RegionHitsSelected()) })
#   RegionGeneCount <- reactive({ plotDatGeneN() })
#   RegionSNPTypeCount <- reactive({ length(unique(plotdatAssoc()$TYPED)) })
#   
#   #Default size per track
#   trackSize <- reactive({ 
#     data.frame(Track=c("Chromosome",
#                        "Manhattan",
#                        "SNPType",
#                        "LD",
#                        "BedGraph",
#                        "wgEncodeBroadHistone",
#                        "wgEncodeRegDnaseClustered",
#                        "annotSmooth",
#                        "Gene"),
#                Size=c(100, #Chromosome
#                       400, #Manhattan
#                       max(30,RegionSNPTypeCount()*20), #SNPType
#                       RegionHitsCount()*20, #LD R^2
#                       100, # BedGraph
#                       40, # wgEncodeBroadHistone
#                       20, # wgEncodeRegDnaseClustered
#                       160, # Annotation
#                       RegionGeneCount()*15) #Gene
#                )})
#   
#   #Create subset based on selected tracks
#   trackHeights <- reactive({
#     trackSize() %>% filter(Track %in% input$ShowHideTracks) %>% .$Size })
#   trackColours <- reactive({ 
#     
#     
#     myCols <- 
#       if(input$PlotTheme=="1"){
#         c(input$PlotThemeColour1,
#           input$PlotThemeColour2) } else { c("#FFFFFF","#FFFFFF") }
#     
# 
#     res <- cbind(trackHeights(),myCols)[,2][1:length(trackHeights())] 
#     #chromosome background must be white
#     if("Chromosome" %in% input$ShowHideTracks){res[1] <- "white"}
#     return(res)
#   })
#   
#   #output plot
#   # See Dynamic UI section
#   
#   # Output to a file --------------------------------------------------------
#   # Get the selected download file type.
#   downloadPlotType <- reactive({input$downloadPlotType})
#   
#   observe({
#     plotType    <- input$downloadPlotType
#     plotTypePDF <- plotType %in% c("pdf","svg")
#     plotUnit    <- ifelse(plotTypePDF, "inches", "pixels")
#     plotUnitDefHeight <- ifelse(plotTypePDF, 12, 1200)
#     plotUnitDefWidth <- ifelse(plotTypePDF, 10, 1000)
#     #plotUnitDefMin <- ifelse(plotTypePDF, 5, 800)
#     #plotUnitDefMax <- ifelse(plotTypePDF, 50, 10000)
#     plotUnitDefStep <- ifelse(plotTypePDF, 1, 100)
#     
#     updateSliderInput( session,
#       inputId = "downloadPlotHeight",
#       label = sprintf("Height (%s)", plotUnit),
#       value = plotUnitDefHeight,
#       step = plotUnitDefStep)
#     
#     updateSliderInput(
#       session,
#       inputId = "downloadPlotWidth",
#       label = sprintf("Width (%s)", plotUnit),
#       value = plotUnitDefWidth,
#       step = plotUnitDefStep)
#   })
#   
#   
#   downloadPlotHeight <- reactive({ input$downloadPlotHeight })
#   downloadPlotWidth <- reactive({ input$downloadPlotWidth })
#   downloadPlotFileName <- reactive({ input$downloadPlotFileName })
#   downloadPlotRes <- reactive({ input$downloadPlotRes })
#   downloadPlotPaper <- reactive({ input$downloadPlotPaper })
#   downloadPlotType <- reactive({ input$downloadPlotType })
#   downloadPlotPointSize <- reactive({ input$downloadPlotPointSize })
#   downloadPlotQuality <- reactive({ input$downloadPlotQuality })
#   downloadPlotTypeJPEG <- reactive({ input$downloadPlotTypeJPEG })
#   downloadPlotTypeCompression <- reactive({ input$downloadPlotCompression })
#   
#   
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
#       
#       switch(input$downloadPlotType,
#              pdf = {plotFunction(con, 
#                                  width = downloadPlotWidth(),
#                                  height = downloadPlotHeight(),
#                                  pointsize = downloadPlotPointSize(),
#                                  paper = downloadPlotPaper(),
#                                  useDingbats = FALSE)},
#              svg = {plotFunction(con, 
#                                  width = downloadPlotWidth(),
#                                  height = downloadPlotHeight(),
#                                  pointsize = downloadPlotPointSize())},
#              jpeg = {plotFunction(con, 
#                                   width = downloadPlotWidth(),
#                                   height = downloadPlotHeight(),
#                                   pointsize = downloadPlotPointSize(),
#                                   res = downloadPlotRes(),
#                                   quality = downloadPlotQuality(),
#                                   type = downloadPlotTypeJPEG())},
#              tiff = {plotFunction(con, 
#                                   width = downloadPlotWidth(),
#                                   height = downloadPlotHeight(),
#                                   pointsize = downloadPlotPointSize(),
#                                   res = downloadPlotRes(),
#                                   type = downloadPlotTypeJPEG(),
#                                   compression = downloadPlotTypeCompression())}
#              
#              )
#       
#       print(plotObjMerge())
#       dev.off(which=dev.cur())
#       }) #downloadHandler downloadPlot
#   
#   #download LDlink processed data as flat file tab separated
#   output$downloadLDFile <- downloadHandler(
#     filename = function() { paste0(datLDlinkProcess()[1,"SNP_A"],"_LD.txt") },
#     content = function(file) {
#       write.table(datLDlinkProcess(), file,
#                   sep="\t",row.names=FALSE, quote=FALSE)
#     })
#   
#   # Dynamic UI --------------------------------------------------------------
#   #Zoom to region X axis BP
#   output$BPrange <-
#     renderUI({
#       sliderInput("BPrange", h5("Region: Start-End"),
#                   min = RegionStart(),
#                   max = RegionEnd(),
#                   value = c(RegionStart(),RegionEnd()),
#                   step = 20000)})
#   #Select hit SNPs
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
                .$CHR %>% unique
              #paste0("chr", c(1:22,"X")),
              #selected = "chr1"
              )
              })
  
  output$RegionID <- renderUI({
    selectInput(inputId = "RegionID", label= h5("Region ID"),
                choices = regions %>%
                  filter(CHR == input$Chr &
                           DATA == input$dataType) %>%
                  .$REGIONBED
                #pre selected a region to demonstrate as an example - HNF1B gene
                #selected = "chr17_35947000_36204000"
                )})
#   
#   
#   #download file name - default: chr_start_end
#   output$downloadPlotFileName <- renderUI({
#     textInput(
#       inputId = "downloadPlotFileName",
#       label = h4("Download file name"),
#       value = paste(RegionChr(),zoomStart(),zoomEnd(),sep="_"))})
#   #download plot title - default: chr:start-end
#   output$downloadPlotTitle <- renderUI({
#     textInput(
#       inputId = "downloadPlotTitle",
#       label = h4("Plot title"),
#       value = paste0(RegionChr(), ':', zoomStart(), '-', zoomEnd(),sep=""))})
#   
#   #Custom BedGraph file track name, default is eQTL
#   output$FileBedGraphName <- renderUI({
#     textInput(inputId = "FileBedGraphName",
#               label = "bedGraph File track name",
#               value = {
#                 inFile <- input$FileBedGraph
#                 if(input$dataType == "Prostate"){"eQTL"} else {
#                   if(is.null(inFile)){"eQTL"} else {
#                     substr(basename(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", 
#                                         inFile$name)),1,15)}}
#                 })#textInput
#     })
# 
# 
#   # Merged Final plot -------------------------------------------------------  
#   #merged plot with dynamic plot height
#   plotObjMerge <- reactive({
#     # Create a Progress object
#     progress <- shiny::Progress$new()
#     progress$set(message = "Plotting please wait...", value = 0)
#     # Close the progress when this reactive exits (even if there's an error)
#     on.exit(progress$close())
#     
#     source("Source/MergePlots.R",local=TRUE)
# 
#     })
#   output$plotMerge <- renderPlot({
#     # Create a Progress object
#     progress <- shiny::Progress$new()
#     progress$set(message = "Almost there...", value = 80)
#     # Close the progress when this reactive exits (even if there's an error)
#     on.exit(progress$close())
#     
#     print(plotObjMerge())
#     
#     })
#   
#   #plot merge height is dynamic, based on seleceted tracks
#   output$plotMergeUI <- renderUI({
#     plotOutput("plotMerge",width=800,height=sum(trackHeights()))
#   })
#     
#   
#     # Observe update ----------------------------------------------------------
#   # maximum of 5 SNPs can be selected to display LD, minimum 1 must be ticked.
#   observe({
#     if(length(input$HitSNPs) > 5){
#       updateCheckboxGroupInput(session, "HitSNPs", selected= head(input$HitSNPs,5))}
#     if(length(input$HitSNPs) < 1){
#       updateCheckboxGroupInput(
#         session, "HitSNPs",
#         selected= RegionHits()[1:min(5,length(RegionHits()))])}
#     if(length(input$ShowHideTracks) < 1){
#       updateCheckboxGroupInput(session, "ShowHideTracks", selected= "Manhattan")}
#   })
#   
#   #manual chr:start-end zoom values
#   observeEvent(input$RegionZoom,({
#     if(!input$RegionZoom %in% c("chr:start-end","")){
#       newStartEnd <- as.numeric(unlist(strsplit(input$RegionZoom,":|-"))[2:3])
#       updateSliderInput(session, "BPrange",
#                         value = newStartEnd)}
#     }))
#   
#   #Reset plot options - 2.Plot Settings
#   observeEvent(input$resetInput,({
#     updateSliderInput(session,"FilterMinPlog", value=0)
#     updateSliderInput(session,"FilterMinLD", value=0.2)
#     
#     updateNumericInput(session,"suggestiveLine", value=5)
#     updateNumericInput(session,"genomewideLine", value=8)
#     
#     updateCheckboxInput(session,"adjustLabels", value=TRUE)
#     updateSliderInput(session,"repFact", value=5)
#     
#     updateSliderInput(session,"BPrange", value = c(RegionStart(),RegionEnd()))
#     updateTextInput(session,"RegionZoom",value="chr:start-end")
#     updateSelectInput(session,"Flank",
#                       selected = 10000)
#     updateCheckboxGroupInput(session,"HitSNPs",
#                              choices = RegionHits(),
#                              selected = 
#                                RegionHits()[1:min(5,length(RegionHits()))])
#     
#     updateCheckboxGroupInput(session,"ShowHideTracks",
#                              selected=c("Manhattan","Recombination"))
#     })
#     ) #END observeEvent resetInput
#   
#   #reset plot output file settings - 3.Final Plot
#   observeEvent(input$resetDownloadPlotSettings |
#                  input$downloadPlotAdvancedSettings,({
#     updateSelectInput(session, "downloadPlotType", selected = "jpeg")
#     updateNumericInput(session,"downloadPlotWidth", value = 1000)
#     updateNumericInput(session,"downloadPlotHeight", value = 1200)
#     updateSliderInput(session,"downloadPlotPointSize", value = 12)
#     updateSelectInput(session,"downloadPlotPaper", selected = "special")
#     updateSliderInput(session,"downloadPlotRes", value = 100)
#     updateSliderInput(session, "downloadPlotQuality", value = 100)
#     updateSelectInput(session, "downloadPlotTypeJPEG", selected = "cairo")
#     updateSelectInput(session, "downloadPlotTypeCompression", selected = "lzw")
#     }) 
#     ) #END observeEvent resetDownloadPlotSettings
# 
#   #If manhattan track is not selected then Recomb and LDSmooth track unticked.
#   observeEvent(input$ShowHideTracks,({
#     selectedTracks <- input$ShowHideTracks
#     if(!"Manhattan"%in% selectedTracks){
#       selectedTracks <- setdiff(selectedTracks,c("Recombination","LDSmooth"))
#       
#       updateCheckboxGroupInput(session,"ShowHideTracks",
#                                selected=c(selectedTracks))
#       } #END if
#     })) #END observeEvent ShowHideTracks
#   
#   #Action buttons to switch between nav bars 1.Input 2.Settings 3.Final
#   observeEvent(input$goToPlotSettings,{
#     updateNavbarPage(session, "navBarPageID", selected = "2.Plot Settings")
#   })
#   observeEvent(input$goToFinalPlot,{
#     updateNavbarPage(session, "navBarPageID", selected = "3.Final Plot")
#   })
#   
#   
#   
  
  
  
})#END shinyServer

