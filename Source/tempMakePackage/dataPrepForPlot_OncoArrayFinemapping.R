# 15/02/2017
# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# plot for Oncoarray selected Finemapping regions for paper

# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(oncofunco) # custom library
#library(cowplot)
library(rtracklayer) # bigwig file subset
#library(ggbio) # tracks

#setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/fix_manhattan_function")
setwd("N:/Translational Cancer Genetics Team/Bioinformatics/Development/R_ShinyApps/LocusExplorer/Source/tempMakePackage")
source("plotFunctions.R")

# Data --------------------------------------------------------------------
#folderLEData <- "C:/Users/tdadaev/Desktop/LocusExplorerPractical_v2/Data"
folderLEData <- "N:/Translational Cancer Genetics Team/Bioinformatics/Development/R_ShinyApps/LocusExplorer/Data"

#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv(paste0(folderLEData,"/wgEncodeBroadHistone/wgEncodeBroadHistone.csv"),
           stringsAsFactors = FALSE)

# Annotation data
#annotOncoFinemap
#annotOncoMeta
annot <- fread(paste0(folderLEData, "/Annotation/annot_OncoArrayFineMapping.txt"))
annot$TYPEN <-
  as.numeric(
    factor(annot$TYPE,
           levels = c(
             "DNaseI","Conserved",
             #ChromHMM
             "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
             "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter")))

# TCGA data
tcgaEQTL <- fread(paste0(folderLEData, "/Annotation/eQTL_OncoArrayFineMapping.tsv"))
#tcgaEQTL <- fread("annot/eQTL_OncoArrayMeta.tsv")
datEQTL <-
  tcgaEQTL %>% 
  filter( geneName != "" &
            !is.na(geneStart) &
            !is.na(BP)) %>% 
  dplyr::transmute(CHR,
                   SNP = rsid,
                   GENE = geneName,
                   x = BP,
                   xend = geneStart,
                   eQTL_TCGA = 1,
                   yend = 0.35) %>%
  base::unique()

regions <- fread("regions_Final.csv")
regionPaper <- regions[ grepl("Paper", Comment), rn ]
# [1] 15 28 38 52 72 80 85
regionJapanese <- regions[ grepl("Japanese", Comment) &
                             grepl("keep", Comment), rn ]

selectedRegions <- regionPaper
#selectedRegions <- setdiff(regionJapanese, regionPaper)

# Loop and plot -----------------------------------------------------------
for(i in selectedRegions) {
  #define region vars
  # i = selectedRegions[1]
  # i=26 #TET2
  # i=49 #RAD23B Clara
  # i=68 #RAD51B Clara
  # i=15 #ANO7/FARP2 Zsofia talk
  # i=80 #KLK        Zsofia talk
  
  # for paper simple region
  # i=38 #RFX chr6, onco_merged_chr6_unphased_116666036_117710052

  
  regionTitle = regions[i, regionNameClean]
  regionName = regions[i, regionNameClean]
  regionChrN = regions[i, chr]
  regionChr = paste0("chr", regionChrN)
  regionStart = regions[i, regionStartSplit]
  regionEnd = regions[i, regionEndSplit]
  regionStartZoom = regions[i, zoomStart]
  regionEndZoom = regions[i, zoomEnd]
  
  # ANO7 eqtl
  #regionStartZoom = 242050000
  #regionEndZoom = 242627000
  
  #KLK eqtl
  #regionStartZoom = 51300000
  #regionEndZoom = 51538500 #zoom1
  #regionEndZoom = 51400000 #zoom2,3
  #regionStartZoom = 51325000 #zoom3
  
  
  #Genomic ranges to subset bigwig data - wgEncodeBroadHistone
  regionGR <-
    GRanges(seqnames = regionChr,
            IRanges(start = regionStart,
                    end = regionEnd))
  
  regionGRzoom <- 
    GRanges(seqnames = regionChr,
            IRanges(start = regionStartZoom,
                    end = regionEndZoom))
  
  fileAssoc <- paste0(folderLEData, "/ProstateData/OncoArrayFineMapping/plotData/",regionName, "_assoc.txt")
  fileLD <- paste0(folderLEData, "/ProstateData/OncoArrayFineMapping/plotData/",regionName, "_LD.txt")
  fileStat <- paste0(folderLEData, "/ProstateData/OncoArrayFineMapping/plotData/",regionName, "_stats.txt")
  
  #data required for manhattan
  assoc <- fread(fileAssoc)
  assoc$PLog <- -log10(assoc$P)
  #zoom
  assoc <- assoc[ BP >= regionStartZoom & BP <= regionEndZoom, ]
  
  # Hits
  stats <- fread(fileStat)
  hits <- list(
    Stepwise = unique(stats[ Method == "Stepwise", SNP]),
    JAMmeta = intersect(
      stats[ Method == "JAM_metaS1_credSet", SNP],
      stats[ Method == "JAM_metaS2_credSet", SNP]),
    JAMonco = intersect(
      stats[ Method == "JAM_oncoS1_credSet", SNP],
      stats[ Method == "JAM_oncoS2_credSet", SNP]),
    StepwiseTag = unique(stats[ grepl("^S.*Tag", Method), SNP]),
    JAMmetaTag = unique(stats[ grepl("^JAM_m.*Tag", Method), SNP]),
    JAMoncoTag = unique(stats[ grepl("^JAM_o.*Tag", Method), SNP]))
  
  #get mean JAM BF for meta, to be desplayed on plot together with SNP name.
  hitsMetaBF <- 
    stats %>%
    filter(rsid %in% hits$JAMmeta & grepl("JAM_metaS._credSet$", Method)) %>% 
    mutate(Method1 = gsub("S1|S2", "", Method)) %>% 
    group_by(rsid) %>% 
      summarise(BF = round(mean(Stats))) %>% 
    mutate(rsid_BF = paste0(rsid, ":BF=", BF ))
  
  
  #region top SNP pvalue log10
  ROIPLogMax <- round(max(-log10(assoc$P)) + 5, -1)
  
  #Recombination map
  geneticMap <- seqminer::tabix.read.table(
    paste0(folderLEData,"/GeneticMap1KG/GeneticMap1KG.txt.gz"),
    paste0(ifelse(regionChr == "chrX", "chr23", regionChr),
           ":", regionStart, "-", regionEnd)) %>%
    transmute(BP = V2,
              RECOMB = V3,
              RECOMB_ADJ = RECOMB * ROIPLogMax / 100)
  # zoom
  geneticMap <- setDT(geneticMap)[BP >= regionStartZoom & BP <= regionEndZoom, ]
  
  #LD data
  LD <- fread(fileLD)
  if(nrow(LD) == 0){
    LD <- data.table(
      CHR_A = regionChrN,
      BP_A = assoc[ order(P), BP][1],
      SNP_A = assoc[ order(P), SNP][1],
      CHR_B = regionChrN,
      BP_B = assoc[ order(P), BP][1],
      SNP_B = assoc[ order(P), SNP][1],
      R2 = 1)}
  #zoom
  LD <- LD[ BP_B >= regionStartZoom & BP_B <= regionEndZoom,]
  
  #FilterMinLD, Onco data used to calc R2, so we have too much data.
  LD <- LD[ R2 >= 0.05]
  
  # bigwig data -------------------------------------------------------------
  #wgEncodeBroadHistone data 7 big wig data
  #subset region
  wgEncodeBroadHistone <-
    bind_rows(
      #seven bigwig files
      lapply(wgEncodeBroadHistoneFileDesc$File, function(i){
        #i=1  wgEncodeBroadHistoneFileDesc$File[1]
        as.data.frame(
          import(paste0(folderLEData,"/wgEncodeBroadHistone/",i),
                 which = regionGRzoom)) %>%
          filter(score >= 5) %>%
          transmute(BP = start,
                    ENCODE = wgEncodeBroadHistoneFileDesc[
                      which(wgEncodeBroadHistoneFileDesc$File == i),
                      "Name"],
                    SCORE = round(ifelse(score >= 100, 100, score),0))
      }))
  #merge on BP, to plot overlappling
  tempBP <- data.frame(BP = unique(wgEncodeBroadHistone$BP))
  wgEncodeBroadHistone <-
    bind_rows(
      lapply(unique(wgEncodeBroadHistone$ENCODE), function(i){
        d <- merge(tempBP,
                   wgEncodeBroadHistone %>%
                     filter(ENCODE == i),
                   by = "BP", all.x = TRUE)
        d$ENCODE <- i
        d$SCORE[ is.na(d$SCORE) ] <- 0
        return(d)}))
  rm(tempBP)
  

  # Manhattan data mDat list object -----------------------------------------
  #Manhattan input data - list
  mDat <- list(
    assoc = assoc,
    geneticMap = geneticMap,
    stats = stats,
    LD = LD,
    suggestiveLine = 5,
    genomewideLine = 8,
    zoomStart = regionStartZoom,
    zoomEnd =   regionEndZoom,
    hits = unique(LD$SNP_A),
    hitsLabel = TRUE)
  # Plot: Locus Explorer ------------------------------------------------------
  # loop through hit types
  #for(hitType in names(hits)[1:2]){
  
  # JAM hits only  
  hitType = names(hits)[2]
  mDat$hits <- hits[[hitType]]
  #mDat$hits <- c(hits[["JAM"]], hits[["JAMtag"]])
  rm(plotManhattan)
  
  # ggManhattan <- 
  #   plotManhattan(data = mDat,
  #                 # opts = c("Recombination","LD","LDSmooth","SuggestiveLine",
  #                 #          "GenomewideLine","Hits"))
  #                 opts = c("Recombination","LD","SuggestiveLine",
  #                          "GenomewideLine","Hits"))
  
  
  
  
  
  ggManhattan <-
    plotManhattan1(assoc = assoc, LD = LD, geneticMap = geneticMap,
                xStart = regionStartZoom, xEnd = regionEndZoom,
                hits = hits$JAMmeta,
                hitsName = hitsMetaBF$rsid_BF[match(hits$JAMmeta, hitsMetaBF$rsid)],
                title = regionTitle, pad = FALSE,
                opts = c("Recombination","LD",
                         "GenomewideLine","Hits")) + theme_LE()
  
  
  #hits <- hits$JAMmeta
  plotDat <- assoc[ assoc$SNP %in% hits$JAMmeta, ]
  plotDat$label <- hitsName[match(plotDat$SNP, hits)]
  
  
  ggplot(plotDat, aes(BP, PLog, label = label)) +
  geom_point() +
  geom_text_repel()
  
  #ggLD <- plotLD(data = mDat)
  
  #ggLDarc1 <- plotLDarc(data = mDat, minR2 = 0.4)
  #ggLDarc2 <- plotLDarc(data = mDat, minR2 = 0.1, hitsOnly = TRUE)
  
  #ggSNPtype <- plotSNPtype(mDat)
  ggSNPtype <-
    plotSNPtype(assoc = assoc, 
                xStart = regionStartZoom, xEnd = regionEndZoom, pad = FALSE) +
    thene_LE()
  
  ggGeneList <- plotGene(chrom = regionChr,
                         chromStart = regionStartZoom,
                         chromEnd = regionEndZoom,
                         vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  ggGene <- ggGeneList$genePlot
  ggGeneCnt <- ggGeneList$geneCnt
  
  # ggAnnot <- plotAnnot(data = annot, chrom = regionChr,
  #                      chromStart = regionStartZoom,
  #                      chromEnd = regionEndZoom,
  #                      vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  
  # eqtl as a box ------------------------------------------------------------
  # merged with annot collapsed
  #    CHR    START      END      FILE            TYPE  COLOUR_RGB COLOUR_HEX TYPEN
  # 1: chr10 45559800 45583200  ChromHMM Heterochromatin 211,211,211    #D3D3D3     3
  # 2: chr10 45583200 45583400  ChromHMM            CTCF   255,165,0    #FFA500     4
  
  # eqtl subset
  plotDatEQTL <- datEQTL %>% 
    filter(CHR == regionChrN &
             x %in% unique(LD[ SNP_A %in% hits$JAMmeta & R2 > 0.4, BP_B])) %>% 
    group_by(GENE) %>% 
    summarise(START = min(x),
              END = max(xend)) %>% 
    transmute(CHR = regionChr,
              START,
              END,
              FILE = "EQTL",
              TYPE = GENE,
              COLOUR_RGB = "255,0,0",
              COLOUR_HEX = "#FF0000",
              TYPEN = 4)
  
  annot_and_eqtl <- rbind(annot, plotDatEQTL)
  setDT(annot_and_eqtl)
  
  ggAnnotCollapsed <- plotAnnot(data = annot_and_eqtl, chrom = regionChr,
                                chromStart = regionStart,
                                chromEnd = regionEnd,
                                collapse = TRUE,
                                vline = mDat$assoc[ SNP %in% mDat$hits, BP]) +
    coord_cartesian(xlim = c(regionStartZoom, regionEndZoom)) +
    #geom_text(data = plotDatEQTL, aes(x = START, y = 4, label = TYPE))
    geom_text_repel(data = plotDatEQTL, aes(x = START, y = 4, label = TYPE))
  
  #add lines to gene track to match with eqtl
  x <- ggplot_build(ggGene)
  x <- sort(x$layout$panel_ranges[[1]]$y.labels)
  geneNames <- seq(x)
  names(geneNames) <- x
  
  # plotDatEQTL$yend = geneNames[plotDatEQTL$TYPE]
  # ggGene <- ggGene +
  #   geom_segment(data = plotDatEQTL, aes(x = START, xend = END,
  #                                        y = ggGeneCnt, yend = yend), col = "red", alpha = 0.5)
  

  
  # ggAnnotCollapsed <- plotAnnot(data = annot, chrom = regionChr,
  #                               chromStart = regionStartZoom,
  #                               chromEnd = regionEndZoom,
  #                               collapse = TRUE,
  #                               vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  
  # ggEQTL <- plotEQTL(data = datEQTL, chromN = regionChrN,
  #                    chromStart = regionStartZoom, chromEnd = regionEndZoom,
  #                    BP_hits = unique(mDat$LD[ SNP_A %in% mDat$hits, BP_A]),
  #                    BP_hits_tag = unique(mDat$LD[ SNP_A %in% mDat$hits & R2 > 0.4, BP_B]))
  
  ggHistone <- plotHistone(data = wgEncodeBroadHistone,
                           chromStart = regionStartZoom,
                           chromEnd = regionEndZoom)
  
  # Plot: JAM -----------------------------------------------------------------  
  fileJAMmeta <- paste0("JAM_meta/20160830_r2_9_BF3_manualPruning_withoutPvalue/",
                        regionNameSplit,"_jam.csv")
  fileJAMmetaRData <- paste0("JAM_meta/20160830_r2_9_BF3_manualPruning_withoutPvalue/",
                             regionNameSplit,".RData")
  fileJAMonco <- paste0("JAM_onco/20161118_r2_9_BF3_manualPruning_withoutPvalue_oncoFullData/",
                        regionNameSplit,"_fullData_jam.csv")
  fileJAMoncoRData <- paste0("JAM_onco/20161118_r2_9_BF3_manualPruning_withoutPvalue_oncoFullData/",
                             regionNameSplit,"_fullData.RData")
  
  # Jam data
  jamMeta <- fread(fileJAMmeta)
  jamOnco <- fread(fileJAMonco)
  
  # hits
  hitsJAMmeta <- jamMeta %>%
    mutate(rn = as.numeric(gsub("SNP_", "", SNPid))) %>%
    arrange(rn) %>%
    filter(PostProb > 0.01 | TopModel95 == 1)
  hitsJAMonco <- jamOnco %>%
    mutate(rn = as.numeric(gsub("SNP_", "", SNPid))) %>%
    arrange(rn) %>%
    filter(PostProb > 0.01 | TopModel95 == 1)
  
  # JAM binary logistic results
  load(fileJAMmetaRData)
  jamMetaResuts <- binary.results
  load(fileJAMoncoRData)
  jamOncoResuts <- logistic.results
  
  # Plot: PostProbs -----------------------------------------------------------
  ggJAM <- list(
    # META: JAM postprob & BF
    metaPP =
      #data
      jamMeta %>% 
      mutate(rsid_BF = if_else(PostProb > 0.05,
                               paste0(rsid, " | BF:", round(BF)),
                               NA_character_),
             col = if_else(TopModel95 == 1, "red", "grey80")) %>% 
      #plot
      ggplot(aes(BP, PostProb, label = rsid_BF, fill = col, col = col)) +
      geom_point() + 
      geom_text_repel(col = "black", na.rm = TRUE) +
      scale_color_identity() + scale_fill_identity() +
      scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0,1)) +
      theme_minimal() +
      ggtitle("meta PostProb > 0.05"),
    # ONCO: JAM postprob & BF
    oncoPP =
      #data
      jamOnco %>% 
      mutate(rsid_BF = if_else(PostProb > 0.05,
                               paste0(rsid, " | BF:", round(BF)),
                               NA_character_),
             col = if_else(TopModel95 == 1, "red", "grey80")) %>% 
      #plot
      ggplot(aes(BP, PostProb, label = rsid_BF, fill = col, col = col)) +
      geom_point() + 
      geom_text_repel(col = "black", na.rm = TRUE) +
      scale_color_identity() + scale_fill_identity() +
      scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0,1)) +
      theme_minimal() +
      ggtitle("onco PostProb > 0.05"),
    # META: AutoCorr plot
    metaAuto = 
      AutocorrelationPlotGG(jamMetaResuts,
                            covariates.to.include = hitsJAMmeta$SNPid,
                            #plot.title = paste0(regionName, ": meta PostProb > 0.01")) +
                            plot.title = "meta PostProb > 0.01") +
      scale_x_discrete(labels = hitsJAMmeta$rsid, name = NULL),
    oncoAuto = 
      AutocorrelationPlotGG(jamOncoResuts,
                            covariates.to.include = hitsJAMonco$SNPid,
                            #plot.title = paste0(regionName, ": onco PostProb > 0.01")) +
                            plot.title = "onco PostProb > 0.01") +
      scale_x_discrete(labels = hitsJAMonco$rsid, name = NULL)
  ) # END ggJAM list
  
  # Plot: Merge --------------------------------------------------------------
  #merge JAM plots
  # outPlot2 <- 
  #   cowplot::plot_grid(cowplot::plot_grid(ggJAM$metaPP,
  #                                         ggJAM$oncoPP, nrow = 1),
  #                      cowplot::plot_grid(ggJAM$metaAuto,
  #                                         ggJAM$oncoAuto, nrow = 1,
  #                                         align = "h"),
  #                      nrow = 2)
  # 
  
  relHeights = 
    c(250, #manhattan
      30, #SNPType
      #35, #eqtl
      30, #histone
      60, #annotcollapsed
      min(c(ggGeneCnt * 20, 120)) #gene
    )
  
  
  
  outPlot1 <- 
    cowplot::plot_grid(
      ggManhattan + ggtitle(regionName),
      ggSNPtype,
      #ggEQTL,
      ggHistone,
      ggAnnotCollapsed,
      ggGene,
      nrow = 5,
      rel_heights = relHeights,
      align = "v")
  
  pdf(paste0("ANO7_",regionTitle, ".pdf"),width = 8.27, height = 11.69, useDingbats = FALSE)
  #pdf(paste0("ANO7_",regionTitle, ".pdf"),width = 6, height = 8, useDingbats = FALSE)
  print(outPlot1)
  dev.off()
  
  outPlot2 <- 
    cowplot::plot_grid(ggJAM$metaPP,
                       ggJAM$metaAuto,
                       nrow = 2
                       # rel_heights = c(relHeights[1],
                       #                 sum(relHeights[2:length(relHeights)]))
    )
  
  
  outPlotFinal <- cowplot::plot_grid(
    outPlot1,
    outPlot2,
    ncol = 2
  )
  
  # Plot: output -------------------------------------------------------------
  pdf(paste0("plot/", regionTitle, ".pdf"),
      width = 11.69, height = 8.27,
      #width = 12, height = 12,
      useDingbats = FALSE)
  print(outPlotFinal)
  # print(outPlot1)
  # print(outPlot2)
  dev.off()
} # END for(i in regions)
