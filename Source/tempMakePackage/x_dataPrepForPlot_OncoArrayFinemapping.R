# 29/09/2016
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

# if(Sys.info()['sysname'] == "Windows") {
#   windowsFonts(Courier=windowsFont("TT Courier New"))
#   setInternet2(use = NA)}

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
    #hits = hits$JAMmeta,
    hitsLabel = TRUE)
  # zoomStart = min(assoc$BP, na.rm = TRUE),
  # zoomEnd =   max(assoc$BP, na.rm = TRUE),
  # hits = hitJAM,
  # hits = sort(c("rs3110641","rs3110644","rs11649743","rs9901746","rs11263763","rs117836847")),
  # hits = c("rs3110641","rs3110644","rs11649743","rs9901746","rs11263763"),
  # hits = sort(c("rs3110644","rs11649743","rs9901746","rs11263763")),
  # hits = c("rs11649743","yy","rs11263763","xx"),
  
  # LD arc plots ---------------
  # mDat1 <- mDat
  # mDat1$stats <- mDat$stats[SNP != "rs13174377",]
  # mDat1$LD <- mDat$LD[SNP_A != "rs13174377" & 
  #                       SNP_B != "rs13174377",]
  # pdf(paste0("plot/arc_", regionTitle, "_R2_2.pdf"), width = 11.69, height = 8.27)
  # plotLDarc(data = mDat, minR2 = 0.2) +
  #   ggtitle(paste0(regionTitle, "; R2 >= 0.2"))
  # dev.off()
  
  # Plot --------------------------------------------------------------------
  
  pdf(paste0("plot/", regionTitle, ".pdf"),width = 8.27, height = 11.69, useDingbats = FALSE)
  #pdf(paste0("TET2_",regionTitle, ".pdf"),width = 8.27, height = 11.69)
  #pdf(paste0("RAD23B_",regionTitle, ".pdf"),width = 8.27, height = 11.69)
  #pdf(paste0("RAD51B_",regionTitle, ".pdf"),width = 8.27, height = 11.69)
  
  #pdf(paste0("ANO7_",regionTitle, ".pdf"),width = 8.27, height = 11.69, useDingbats = FALSE)
  #jpeg(paste0("ANO7_",regionTitle, ".jpeg"),width = 1000, height = 1400,
  #     res = 100, quality = 100, pointsize = 12)
  
  #pdf(paste0("KLK_",regionTitle, ".pdf"),width = 8.27, height = 11.69, useDingbats = FALSE)
  #svg(paste0("KLK_",regionTitle, ".svg"),width = 8.27, height = 11.69)
  # jpeg(paste0("KLK_",regionTitle, ".jpeg"),width = 1000, height = 1400,
  #      res = 100, quality = 100, pointsize = 12)
  # loop through hit types
  for(hitType in names(hits)[1:2]){
    #hitType = names(hits)[2]
    mDat$hits <- hits[[hitType]]
    #mDat$hits <- c(hits[["JAM"]], hits[["JAMtag"]])
    
    ggManhattan <- 
      plotManhattan(data = mDat,
                    opts = c("Recombination","LD","LDSmooth","SuggestiveLine",
                             "GenomewideLine","Hits"))
    
    ggLD <- plotLD(data = mDat)
    
    ggLDarc1 <- plotLDarc(data = mDat, minR2 = 0.4)
    #ggLDarc2 <- plotLDarc(data = mDat, minR2 = 0.1, hitsOnly = TRUE)
    
    # temp for plotting Arc plots separately
    # pdf(paste0("RAD51B_",regionTitle, "_arc.pdf"),width = 11.69, height = 8.27)
    # pdf(paste0("RAD23B_",regionTitle, "_arc.pdf"),width = 11.69, height = 8.27)
    # plotLDarc(data = mDat, minR2 = 0.2) + ggtitle("RAD23B, R2 >= 0.2")
    # plotLDarc(data = mDat, minR2 = 0.4) + ggtitle("RAD23B, R2 >= 0.4")
    # plotLDarc(data = mDat, minR2 = 0.1, hitsOnly = TRUE) + ggtitle("RAD23B, R2 >= 0.1, hits only")
    # dev.off()
    
    ggSNPtype <- plotSNPtype(mDat)
    
    
    ggGeneList <- plotGene(chrom = regionChr,
                           chromStart = regionStartZoom,
                           chromEnd = regionEndZoom,
                           vline = mDat$assoc[ SNP %in% mDat$hits, BP])
    ggGene <- ggGeneList$genePlot
    ggGeneCnt <- ggGeneList$geneCnt
    
    ggAnnot <- plotAnnot(data = annot, chrom = regionChr,
                         chromStart = regionStartZoom,
                         chromEnd = regionEndZoom,
                         vline = mDat$assoc[ SNP %in% mDat$hits, BP])
    
    ggAnnotCollapsed <- plotAnnot(data = annot, chrom = regionChr,
                                  chromStart = regionStartZoom,
                                  chromEnd = regionEndZoom,
                                  collapse = TRUE,
                                  vline = mDat$assoc[ SNP %in% mDat$hits, BP])
    
    #BP_hits = mDat$assoc[ SNP %in% mDat$hits, BP]
    #BP_hits = mDat$LD[, BP_B]
    ggEQTL <- plotEQTL(data = datEQTL, chromN = regionChrN,
                       chromStart = regionStartZoom, chromEnd = regionEndZoom,
                       BP_hits = unique(mDat$LD[ SNP_A %in% mDat$hits, BP_A]),
                       BP_hits_tag = unique(mDat$LD[ SNP_A %in% mDat$hits & R2 > 0.4, BP_B]))
    
    ggHistone <- plotHistone(data = wgEncodeBroadHistone,
                             chromStart = regionStartZoom,
                             chromEnd = regionEndZoom)
    
    
    # Merged Plot -------------------------------------------------------------
    # using ggbio::tracks
    print(
      tracks(ggManhattan,
             ggLD, 
             #ggLDarc1,
             #ggLDarc2,
             ggSNPtype,
             ggEQTL,
             ggHistone,
             #ggAnnot,
             ggAnnotCollapsed,
             ggGene,
             heights = c(300, #manhattan
                         max(c(15, length(mDat$hits) * 10)), #LD
                         #80,
                         #100,
                         25, #SNPType
                         70, #eqtl
                         25, #histone
                         #110, #annot
                         35, #annotcollapsed
                         min(c(ggGeneCnt * 10, 120)) #gene
             ), #trackHeights(),
             track.plot.color = c("grey70", "grey80","grey70","grey80",
                                  "grey70", "grey80", "grey70"), #trackColours(),
             title = paste0(regionTitle, " Hits: ", hitType),
             xlim = c(regionStartZoom, regionEndZoom)
             #padding = unit(-0.5, "lines")
      ) +
        #input$downloadPlotTitle) +
        # change Y axis text to defualt font
        theme(axis.text.y = element_text(family = ''))
    ) # print tracks
  } # END for(hitType in names(hits))
  dev.off()
} # END for(i in regions)



# TESTING -----------------------------------------------------------------


# test cowplot -----
# cowplot: Switching axis messes up plot_grid alignment
#ggManhattan <- ggdraw(switch_axis_position(ggManhattan, axis = "x"))


# theme(legend.position = "none",
#       panel.background = element_rect(fill = "white"), 
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(), 
#       panel.grid.major.y = element_line(colour = "grey60", linetype = "dotted"),
#       axis.title.x = element_blank(), 
#       axis.line = element_blank(),
#       panel.border = element_blank(), 
#       axis.text.y = element_text(family = "Courier", colour = "grey20"))



# annotOncoFinemap %>% 
#   transmute(CHR = Chr,
#             BP = End,
#             DNaseI =
#               as.numeric((GSM1008595_RWPE1_DNaseI +
#                             GSM1024742_PrEC_DNaseI +
#                             GSM1024743_PrEC_DNaseI +
#                             GSM736565_LNCaP_DNaseI +
#                             GSM736603_LNCaP_DNaseI +
#                             `GSM816634_LNCaP+androgen_DNaseI` +
#                             GSM816637_LNCaP_DNaseI) > 0),
#             Conserved =
#               as.numeric((`GERP++` +
#                             SiPhy_Omega +
#                             SiPhy_Pi +
#                             PhastCons) > 0),
#             Heterochromatin = PrEC_ChromHMM_Heterochromatin,
#             CTCF = PrEC_ChromHMM_CTCF,
#             CTCF_Enhancer = 0,
#             Promoter = PrEC_ChromHMM_Promoter,
#             Enhancer = PrEC_ChromHMM_Enhancer,
#             Poised_Promoter = PrEC_ChromHMM_Poised_Promoter,
#             Transcribed = PrEC_ChromHMM_Transcribed,
#             Repressed = PrEC_ChromHMM_Repressed,
#             CTCF_Promoter = 0
#             ) %>% head
# 
# write.csv(head(annotOncoFinemap), "temp.csv")


# Merge = dont use cowplot, can't align ggbio genes track... :(
# plot_grid(ggManhattan,
#           ggLD,
#           ggSNPtype,
#           ggGene,
#           ncol = 1, rel_heights = c(0.6,0.1,0.1,0.3), align = "hv")


# test tracks backup ----------
# tracks(ggManhattan,
#        ggLD, 
#        ggSNPtype,
#        #ggEQTL,
#        ggHistone,
#        #ggAnnot,
#        ggAnnotCollapsed,
#        ggGene,
#        heights = c(300, #manhattan
#                    20, #LD
#                    30, #SNPType
#                    #50, #eqtl
#                    25, #histone
#                    #110, #annot
#                    40, #annotcollapsed
#                    120 #gene
#        ), #trackHeights(),
#        track.plot.color = c("grey70", "grey80","grey70","grey80",
#                             "grey70", "grey80"), #trackColours(),
#        title = regionTitle,
#        xlim = c(regionStartZoom, regionEndZoom)) +
#   #input$downloadPlotTitle) +
#   # change Y axis text to defualt font
#   theme(axis.text.y = element_text(family = ''))



#regionChr = unlist(strsplit(regionName,"_"))[1]
#regionChrN = ifelse(regionChr == "chrX", 23, as.numeric(sub("chr","",regionChr)))
#regionStart = as.numeric(unlist(strsplit(regionName,"_"))[2])
#regionEnd = as.numeric(unlist(strsplit(regionName,"_"))[3])


# test Haploreg -------------
# #haploreg <- fread(paste0(folderLEData,"/Haploreg/20161027.txt"))
# haploreg[haploreg == "." | haploreg == ""] <- NA
# haploreg <- haploreg %>%
#         transmute(
#           CHR = chr,
#           SNP = rsID,
#           SiPhy_cons,
#           EnhancerHistoneMarks = as.integer(!is.na(Chromatin_Marks)),
#           DNAse = as.integer(!is.na(DNAse)),
#           ProteinBound = as.integer(!is.na(Proteins)),
#           MotifsChanged = as.integer(!is.na(Motifs)),
#           GERP_cons = GERP_cons,
#           eQTL = as.integer(is.na(grasp) + is.na(eQTL) < 2),
#           functional = dbSNP_functional_annotation)substr(dbSNP_functional_annotation, 1, 3)) %>%
#         gather(annot,value,SiPhy_cons:functional) %>%
#         filter(value != 0) %>%
#         mutate(annot = if_else(value == 1,
#                                as.character(annot), paste0("Func_", value))) %>%
#         dplyr::select(CHR, SNP, ANNOT = annot)
# setDT(haploreg)
# # table(haploreg$annot)


#haploreg = haploreg[ CHR == regionChrN, ],
# annot = annot[ CHR == regionChr &
#                  START >= regionStartZoom &
#                  END <= regionEndZoom &
#                  # c("Heterochromatin","CTCF","CTCF+Enhancer",
#                  #   "Promoter","Enhancer","Poised_Promoter",
#                  #   "Transcribed","Repressed","CTCF+Promoter",
#                  #   "DNaseI","Conserved")
#                  TYPE %in%  c("Heterochromatin","CTCF","CTCF+Enhancer",
#                               "Promoter","Enhancer","Poised_Promoter",
#                               "Transcribed","Repressed","CTCF+Promoter")],
# c("Heterochromatin","CTCF","CTCF+Enhancer",
#   "Promoter","Enhancer","Poised_Promoter",
#   "Transcribed","Repressed","CTCF+Promoter",
#   "DNaseI","Conserved")
# TYPE %in%  c("Heterochromatin","CTCF","CTCF+Enhancer",
#              "Promoter","Enhancer","Poised_Promoter",
#              "Transcribed","Repressed","CTCF+Promoter")]






# selected region for main paper for Zsofia - zoom
# regionSelected <- c("chr11_107643000_108644000", # 107800000,108500000
#                     "chr12_12371000_13372000", #12700000,13000000
#                     "chr15_55886000_56887000", #56000000,
#                     "chr12_64513000_65514000") #56800000

# regionSelected <- paste(regions$CHR, regions$START, regions$END, sep = "_")
# regionSelected <- regionSelected[ regions$DATA == "OncoArrayFineMapping"]
# regionExclude <- c("chr2_9611973_11210730",
#                    "chr3_86610674_87967332",
#                    "chr4_20378153_21378153",
#                    "chr5_780028_2395829",
#                    "chr6_29573776_30573776",
#                    "chr6_30618511_32900939",
#                    "chr8_127333841_129040776",
#                    "chr17_46302314_47952263")
#regionSelected <- regionSelected[!regionSelected %in% regionExclude]

# test Regions for paper ---------
# onco_merged_chr2_unphased_241657087_242920971,rs3771570_2q37_(FARP2)
# onco_merged_chr6_unphased_116666036_117710052,rs339331_6q22_(RFX6)
# onco_merged_chr10_unphased_51049496_52049496,rs10993994_10q11_(MSMB)
# onco_merged_chr17_unphased_35547276_36603565,rs11649743_17q12_(HNF1B); rs4430796_17q12_(HNF1B)
# onco_merged_chr19_unphased_50840794_51864623,rs2735839_19q13_(KLK2/3)
# onco_merged_chr21_unphased_42401421_43401421,rs1041449_21q22_(TMPRSS2)
# onco_merged_chr5_unphased_780028_2395829_1,Split:TERT
# regionsPaper <- c("chr2_241657087_242920971",
#                   "chr6_116666036_117710052",
#                   "chr10_51049496_52049496",
#                   "chr17_35547276_36603565",
#                   "chr19_50840794_51864623",
#                   "chr21_42401421_43401421",
#                   "chr5_780028_1600000")



# test Japanese hits ------------
# [1]  9 17 18 28 29 35 38 54 57 65
# onco_merged_chr2_unphased_20388265_21388265	keep:Japanese_Takata
# onco_merged_chr3_unphased_86610674_87967332	keep:split:Japanese_Akamatsu
# onco_merged_chr3_unphased_86610674_87967332	keep:split:Japanese_Akamatsu
# onco_merged_chr5_unphased_780028_2395829	keep:split:forPaperMainText:Japanese_Takata
# onco_merged_chr5_unphased_780028_2395829	keep:split:Japanese_Takata
# onco_merged_chr6_unphased_41036427_42043793	keep:Japanese_Takata
# onco_merged_chr6_unphased_116666036_117710052	keep:forPaperMainText:Japanese_Takata
# onco_merged_chr10_unphased_122283141_123344709	keep:Japanese_Akamatsu
# onco_merged_chr11_unphased_58415110_59610571	keep:Japanese_Akamatsu
# onco_merged_chr13_unphased_73228139_74468916	keep:Japanese_Takata
