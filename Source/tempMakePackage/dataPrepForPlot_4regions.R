# 11/08/2016
# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# ATM - chr11_107643000_108644000
# regionName = "chr11_107643000_108644000"

# CDKN1B - chr12_12371000_13372000
# regionName = "chr12_12371000_13372000"

# RFX7 - chr15_55886000_56887000
# regionName = "chr15_55886000_56887000"

# RASSF3 - chr12_64513000_65514000
# regionName = "chr12_64513000_65514000"


# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(oncofunco) # custom library
#library(cowplot)
library(rtracklayer) # bigwig file subset

# if(Sys.info()['sysname'] == "Windows") {
#   windowsFonts(Courier = windowsFont("TT Courier New"))
#   setInternet2(TRUE)
# }

setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/fix_manhattan_function")
source("plotFunctions.R")

# Data --------------------------------------------------------------------
regions <- fread("regions.csv")
folderLEData <- "C:/Users/tdadaev/Desktop/LocusExplorerPractical_v2/Data"

#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv(paste0(folderLEData,"/wgEncodeBroadHistone/wgEncodeBroadHistone.csv"),
           stringsAsFactors = FALSE)

# oncoarray finemapping annotation + TCGA
#annotOncoFinemap
#annot <- fread("annot/annotOncoFinemap.txt")
#annotOncoMeta
annot <- fread("annot/annotOncoMeta.txt")
annot$TYPEN <-
  as.numeric(
    factor(annot$TYPE,
           levels = c(
             "DNaseI","Conserved",
             #ChromHMM
             "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
             "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter")))

#tcgaEQTL <- fread("annot/LE_finemap_eQTL.tsv")
tcgaEQTL <- fread("annot/OA_NewHits_eQTLs.tsv")
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

# selected region for main paper for Zsofia
regionSelected <- c("chr11_107643000_108644000", # 107800000,108500000
                    "chr12_12371000_13372000", #12700000,13000000
                    "chr15_55886000_56887000", #56000000,
                    "chr12_64513000_65514000") #56800000

# Loop and plot -----------------------------------------------------------
for(i in regionSelected) {
  
  #define region vars
  #i = regionSelected[4]
  regionTitle = i
  regionName = i
  regionChr = unlist(strsplit(regionName,"_"))[1]
  regionChrN = as.numeric(sub("chr","",regionChr))
  regionStart = as.numeric(unlist(strsplit(regionName,"_"))[2])
  regionEnd = as.numeric(unlist(strsplit(regionName,"_"))[3])
  regionStartZoom = 64900000
  regionEndZoom =   65500000
  
  
  #Genomic ranges to subset bigwig data - wgEncodeBroadHistone
  regionGR <-
    GRanges(seqnames = regionChr,
            IRanges(start = regionStart,
                    end = regionEnd))
  
  
  fileAssoc <- paste0(folderLEData, "/ProstateData/OncoArrayMeta/LE/",regionName, "_assoc.txt")
  fileLD <- paste0(folderLEData, "/ProstateData/OncoArrayMeta/LE/",regionName, "_LD.txt")
  
  #data required for manhattan
  assoc <- fread(fileAssoc)
  assoc$PLog <- -log10(assoc$P)
  
  #region top SNP pvalue log10
  ROIPLogMax <- round(max(-log10(assoc$P)) + 5, -1)
  
  #Recombination map
  geneticMap <- seqminer::tabix.read.table(
    "GeneticMap1KG.txt.gz",
    paste0(regionChr,":", regionStart, "-", regionEnd)) %>%
    transmute(BP = V2,
              RECOMB = V3,
              RECOMB_ADJ = RECOMB * ROIPLogMax / 100)
  
  #LD data
  LD <- fread(fileLD)
  #FilterMinLD
  #RegionHitsSelected
  
  
  # Manhattan data mDat list object -----------------------------------------
  #Manhattan input data - list
  mDat <- list(
    assoc = assoc,
    geneticMap = geneticMap,
    LD = LD,
    suggestiveLine = 5,
    genomewideLine = 8,
    zoomStart = regionStartZoom, # min(assoc$BP, na.rm = TRUE),
    zoomEnd =   regionEndZoom, # max(assoc$BP, na.rm = TRUE),
    hits = unique(LD$SNP_A),
    # hits = sort(c("rs3110641","rs3110644","rs11649743","rs9901746","rs11263763","rs117836847")),
    # hits = c("rs3110641","rs3110644","rs11649743","rs9901746","rs11263763"),
    # hits = sort(c("rs3110644","rs11649743","rs9901746","rs11263763")),
    # hits = c("rs11649743","yy","rs11263763","xx"),
    hitsLabel = TRUE)
  

  
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
                 which = regionGR)) %>%
          filter(score >= 5) %>%
          transmute(BP = start,
                    ENCODE = wgEncodeBroadHistoneFileDesc[
                      which(wgEncodeBroadHistoneFileDesc$File == i),
                      "Name"],
                    SCORE = round(ifelse(score >= 100, 100, score),0))
      }))
  #merge on BP, to plot overlappling
  res <- data.frame(BP = unique(wgEncodeBroadHistone$BP))
  wgEncodeBroadHistone <-
    bind_rows(
      lapply(unique(wgEncodeBroadHistone$ENCODE), function(i){
        d <- merge(res,
                   wgEncodeBroadHistone %>%
                     filter(ENCODE == i),
                   by = "BP", all.x = TRUE)
        d$ENCODE <- i
        d$SCORE[ is.na(d$SCORE) ] <- 0
        return(d)
      }))
  rm(res)

  
  # Plot --------------------------------------------------------------------
  ggManhattan <- 
    plotManhattan(data = mDat,
                  opts = c("Recombination","LD","LDSmooth","SuggestiveLine","GenomewideLine","Hits"))
  # Switching axis messes up plot_grid alignment
  #ggManhattan <- ggdraw(switch_axis_position(ggManhattan, axis = "x"))
  
  ggLD <- plotLD(data = mDat, opts = c("LD"))
  
  ggSNPtype <- plotSNPtype(mDat, opts = c("LD"))
  
  ggGene <- plotGene(chrom = regionChr,
                     chromStart = regionStartZoom,
                     chromEnd = regionEndZoom,
                     vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  
  ggAnnot <- plotAnnot(data = annot, chrom = regionChr,
                       chromStart = regionStartZoom,
                       chromEnd = regionEndZoom,
                       #collapse = TRUE,
                       vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  
  ggAnnotCollapsed <- plotAnnot(data = annot, chrom = regionChr,
                       chromStart = regionStartZoom,
                       chromEnd = regionEndZoom,
                       collapse = TRUE,
                       vline = mDat$assoc[ SNP %in% mDat$hits, BP])
  
  ggEQTL <- plotEQTL(data = datEQTL, chromN = regionChrN,
                     chromStart = regionStartZoom, chromEnd = regionEndZoom,
                     #BP_hits = mDat$assoc[ SNP %in% mDat$hits, BP]
                     #BP_hits = mDat$LD[ R2 > 0.8, BP_B]
                     BP_hits = mDat$LD[, BP_B])
  
  ggHistone <- plotHistone(data = wgEncodeBroadHistone,
                           chromStart = regionStartZoom,
                           chromEnd = regionEndZoom)
  
  
  # Merged Plot -------------------------------------------------------------
  # using ggbio::tracks
  pdf(paste0("plot/", i, ".pdf"),width = 8.27, height = 11.69)
  tracks(ggManhattan,
         ggLD, #ggSNPtype,
         ggEQTL,
         ggHistone,
         #ggAnnot,
         ggAnnotCollapsed,
         ggGene,
         heights = c(300, #manhattan
                     25, #LD
                     50, #eqtl
                     25, #histone
                     #110, #annot
                     40, #annotcollapsed
                     120 #gene
                     ), #trackHeights(),
         track.plot.color = c("grey70", "grey80","grey70","grey80",
                              "grey70", "grey80"), #trackColours(),
         title = regionTitle,
         xlim = c(regionStartZoom, regionEndZoom)) +
   #input$downloadPlotTitle) +
    # change Y axis text to defualt font
    theme(axis.text.y = element_text(family = ''))
  dev.off()
  
  
} #END for(i in regionSelected)



# TESTING -----------------------------------------------------------------

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

