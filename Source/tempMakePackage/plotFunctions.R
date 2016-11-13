# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt
# About: function version of LocusExplorer, with no shiny.


# ggThemeLE ---------------------------------------------------------------
# custom theme
# ggThemeLE <- 
#   theme(axis.line = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(colour = "grey20"))



# plotManhattan -----------------------------------------------------------
plotManhattan <- function(data, opts){
  
  # Check input data --------------------------------------------------------
  # .... neeed adding
  # data = mDat
  
  
  # list of hits matching assoc and LD input, and sorted alpha
  hits <- unique(sort(
    intersect(
      intersect(data$hits, data$assoc$SNP),
      data$LD$SNP_A)))
  
  # Warning if input hit SNPs are not in assoc LD files
  if(length(hits) > 0) {
    if(!all(data$hits %in% hits)){ 
      warning(paste0("Some SNPs (",
                     paste(setdiff(data$hits,hits), collapse = ","),
                     ") did not match to assoc and LD data."))}}
  
  # Plot all SNPs - background ----------------------------------------------
  gg_out <-
    if(!is.null(data$assoc)){
      ggplot(data$assoc, aes(x = BP, y = PLog)) +
        # all snps grey hollow shapes
        geom_point(size = 4, colour = "#B8B8B8", shape = data$assoc$TYPED) +
        geom_hline(yintercept = seq(0,ROIPLogMax,5),
                   linetype = "dotted", col = "grey60")
      
    } else { stop("Association data is missing.") }
  
  # Recombination -----------------------------------------------------------
  if(!is.null(data$geneticMap))
    if("Recombination" %in% opts & nrow(data$geneticMap) > 2 ){
      #Recomb lines
      gg_out <- gg_out +
        geom_area(data = data$geneticMap,
                  aes(BP, RECOMB_ADJ),
                  fill = "#11d0ff", colour = "#00B4E0", alpha = 0.3)}
  
  # LD Fill -----------------------------------------------------------------
  if(length(hits) > 0 & !is.null(data$LD)){
    #colours for hits and pallete for R2 shades
    colourLD <- oncofunco::colourHue(length(hits))
    colourLDPalette <- unlist(lapply(colourLD, function(i){
      colorRampPalette(c("grey95", i))(100)}))
    
    #merge LD with assoc, to get R2 shades per point
    plotDat <- merge(
      data$LD[SNP_A %in% hits, ],
      data$assoc[, list(BP_B = BP, TYPED, PLog)],
      by = "BP_B") %>%
      # data$assoc[, list(SNP_B = SNP, TYPED, PLog)],
      # by = "SNP_B") %>%
      tbl_df %>%
      mutate(
        LDColIndex = ifelse(round(R2,2) == 0, 1, round(R2, 2) * 100),
        hitColIndex = as.numeric(factor(SNP_A, levels = hits)),
        hitCol = colourLD[hitColIndex],
        LDCol = colourLDPalette[(hitColIndex - 1) * 100 + LDColIndex],
        R2Adj = ROIPLogMax * R2 * 0.8)
    #LD fill
    if("LD" %in% opts & nrow(data$LD) > 10){
      gg_out <- gg_out +
        geom_point(data = plotDat[ plotDat$R2 > 0.1, ], aes(BP_B, PLog, colour = LDCol),
                   size = 4,
                   shape = plotDat[ plotDat$R2 > 0.1, ]$TYPED + 15,
                   col = plotDat[ plotDat$R2 > 0.1, ]$LDCol,
                   alpha = 0.8)
    } # LD fill
  } # if LD null
  
  # LD Smooth ----------------------------------------------------------------
  if(length(hits) > 0  & !is.null(data$LD))
    if("LDSmooth" %in% opts & nrow(data$LD) > 10){
      gg_out <- gg_out +
        geom_smooth(data = plotDat, aes(x = BP_B, y = R2Adj, col = hitCol),
                    method = "loess", se = FALSE)}
  
  # Suggestiveline ----------------------------------------------------------
  if(!is.null(data$suggestiveLine))
    if("SuggestiveLine" %in% opts & data$suggestiveLine > 0) {
      gg_out <- gg_out +
        geom_hline(aes(yintercept = y), data = data.frame(y = data$suggestiveLine),
                   size = 0.5,
                   colour = "#1a9641")}
  # Genomewideline ----------------------------------------------------------
  if(!is.null(data$genomewideLine))
    if("GenomewideLine" %in% opts & data$genomewideLine > 0) {
      gg_out <- gg_out +
        geom_hline(aes(yintercept = y), data = data.frame(y = data$genomewideLine),
                   size = 0.5,
                   colour = "#ca0020")}
  
  # Mark Hits: shape and vline ----------------------------------------------
  if("Hits" %in% opts & length(hits) > 0){
    gg_out <- gg_out +
      #mark hit SNPs - outline shapes
      geom_point(data = data$assoc[ SNP %in% hits, ],
                 aes(x = BP, y = PLog), size = 4, colour = "black",
                 shape = data$assoc[ SNP %in% hits, TYPED]) +
      #mark hit SNPs - vertical lines
      geom_segment(data = data$assoc[ SNP %in% hits, ],
                   aes(x = BP, y = 0, xend = BP, yend = PLog),
                   colour = "black",
                   linetype = "dashed")}
  
  
  # Mark Hits: Labels -------------------------------------------------------
  # SNP names on the plot for hits
  if("Hits" %in% opts & length(hits) > 0)
    if(!is.null(data$hitsLabel))
      if(data$hitsLabel){
        gg_out <- gg_out +
          geom_text_repel(
            aes(BP, PLog, label = SNP),
            data = data$assoc[ data$assoc$SNP %in% hits, ])}
  
  # Annotations -------------------------------------------------------------
  # if("Haploreg" %in% opts & nrow(data$haploreg) > 0){
  #   haploregPlotDat <- merge(data$assoc[, list(SNP, BP, PLog)], 
  #                            data$haploreg[ SNP %in% data$hits ], by = "SNP")
  #   typeN <- length(unique(haploregPlotDat$ANNOT))
  #   1:typeN
  #   data$zoomStart + 1:typeN * (data$zoomEnd - data$zoomStart)/typeN
  #   
  #   # haploregPlotDat <- 
  #   #   haploregPlotDat[, list(ANNOT = paste(ANNOT, collapse = ",")), by = SNP]
  # 
  #   gg_out <- gg_out +
  #     geom_text_repel(
  #       aes(BP, PLog, label = ANNOT), col = "red",
  #       #arrow = arrow(length = unit(0.02, "npc"),angle = 45),
  #       data = haploregPlotDat)
  #     
  #   
  # }
  # 
  
  # General options ---------------------------------------------------------
  # Zoom
  gg_out <- gg_out +
    xlim(c(data$zoomStart, data$zoomEnd)) +
    scale_y_continuous(
      limits = c(0, ROIPLogMax),
      breaks = seq(0, ROIPLogMax, 5),
      #labels = oncofunco::strPadLeft(seq(0, ROIPLogMax, 5)),
      labels = seq(0, ROIPLogMax, 5),
      name = expression(-log[10](italic(p)))) +
    scale_colour_identity() +
    theme_LE()
  
  # xlab("") +
  # theme(axis.line = element_blank(),
  #       axis.title.x = element_text(size = 0)) +
  #ggThemeLE() #+
  #theme(axis.line = element_blank(),
  #      axis.title.x = element_blank())
  #theme(
  #Y Axis font
  #axis.text.y = element_text(family = "Courier"),
  #panel.background = element_rect(fill = "white")) 
  
  
  # Output ------------------------------------------------------------------
  gg_out
}


# plotLD ------------------------------------------------------------------
plotLD <- function(data, opts){
  
  # Check input data --------------------------------------------------------
  # .... neeed adding
  
  
  # list of hits matching assoc and LD input, and sorted alpha
  hits <- sort(intersect(data$hits, data$LD$SNP_A))
  
  if(length(hits) > 0){
    
    # Warning if input hit SNPs are not in assoc LD files
    if(!all(data$hits %in% hits)){ 
      warning(paste0("Some SNPs (",
                     paste(setdiff(data$hits,hits), collapse = ","),
                     ") did not match to assoc and LD data."))}
    
    # LD track
    # plot LD per hit SNPs on seperate Yaxis 1,2,3,etc
    if(!is.null(data$LD)){
      #colours for hits and pallete for R2 shades
      colourLD <- oncofunco::colourHue(length(hits))
      colourLDPalette <- unlist(lapply(colourLD, function(i){
        colorRampPalette(c("grey95", i))(100)}))
      
      # LD with R2 shades per segment matching Manhattan plot colours
      plotDat <- data$LD %>% 
        filter(SNP_A %in% data$hits) %>% 
        transmute(
          x = BP_B,
          xend = BP_B,
          y = as.numeric(factor(SNP_A, levels = hits)),
          yend = y + 1,
          LDColIndex = ifelse(round(R2,2) == 0, 1, round(R2, 2) * 100),
          hitColIndex = as.numeric(factor(SNP_A, levels = hits)),
          hitCol = colourLD[hitColIndex],
          LDCol = colourLDPalette[(hitColIndex - 1) * 100 + LDColIndex],
          R2Adj = ROIPLogMax * R2 * 0.8)
      
      # plot segments
      ggplot(data = plotDat,
             aes(x = x, xend = xend,
                 y = y, yend = yend,
                 colour = LDCol)) +
        geom_hline(yintercept = 1:length(hits) + 0.5,
                   linetype = "dotted", col = "grey60") +
        geom_segment() +
        scale_colour_identity() +
        scale_y_continuous(breaks = (1:length(hits)) + 0.5,
                           labels = substr(hits,1,20),
                           limits = c(0, length(hits) + 1),
                           #labels = oncofunco::strPadLeft(hits),
                           #name="xxx"
                           name = expression(R^2)) +
        #general options
        xlim(c(data$zoomStart, data$zoomEnd)) +
        theme_LE()
    } # END if is null LD data
  } else {
    # plot blank if hits are missing
    ggplot() + 
      geom_blank() +
      annotate("text",
               x = chromStart + (chromEnd - chromStart)/2,
               y = 0.6,
               label = "No hits") +
      scale_y_continuous(breaks = c(0.6),
                         labels = "",
                         limits = c(-0.25, 1),
                         name = expression(R^2)) +
      xlim(c(chromStart, chromEnd)) +
      theme_LE() +
      theme(axis.text.x = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.ticks.x = element_blank())
    
  }
}


# plotLDarc ---------------------------------------------------------------
plotLDarc <- function(data, opts, minR2 = 0.2, hitsOnly = FALSE){

  #function(LD, stats, minR2 = 0.2, hitsOnly = FALSE){
  #minR2 = 0.1; hitsOnly = TRUE
  # LD object, plink style
  # stats, LE stats output file
  # minR2, filter LD data
  # hitsOnly, SNP_A and SNP_B in hit list
  
  # data <- mDat; minR2 = 0.2
  # data prep --------------
  hitStepwise <- unique(data$stats[ Method == "Stepwise", SNP])
  hitJAM <- unique(data$stats[ Method == "JAM", SNP])
  hits <- unique(c(hitStepwise, hitJAM))
  
  # label <- filter(stats, Method %in% c("Stepwise", "JAM")) %>% 
  #   dplyr::select(SNP, BP) %>% 
  #   unique %>% 
  #   arrange(Xstart)
  # label$Xend <- seq(min(label$Xstart), max(label$Xstart), length.out = nrow(label))
  
  arcDat <- 
    list(arc =
           rbind(
             data$LD %>% 
               filter(BP_A != BP_B &
                        R2 >= minR2 & 
                        SNP_A %in% hitStepwise) %>% 
               mutate(From = pmax(BP_A, BP_B),
                      To = pmin(BP_A, BP_B),
                      R2 = round(R2, 2)),
             data$LD %>% 
               filter(BP_A != BP_B &
                        R2 >= minR2 & 
                        SNP_A %in% hitJAM) %>% 
               mutate(From = pmin(BP_A, BP_B),
                      To = pmax(BP_A, BP_B),
                      R2 = round(R2, 2))
             
           ) %>% arrange(R2),
         label = filter(data$stats, Method %in% c("Stepwise", "JAM")) %>% 
           transmute(SNP = substr(SNP, 1, 13),
                     BP,
                     #Ylabel = if_else(Method == "Stepwise", 0.75, 0.25),
                     Y = if_else(Method == "Stepwise", 0.75, 0.25),
                     Yend = 0.5
           ) %>% 
           unique %>%  mutate(Y = if_else(row_number() %% 2 == 0, Y - 0.05, Y))
    )
  
  if(hitsOnly){
    arcDat$arc <- arcDat$arc %>% filter(SNP_A %in% hits & SNP_B %in% hits)
  }
  
  # plot ------------
  ggplot(arcDat$arc, aes(x = From, xend = To, y = 0.5, yend = 0.5, col = R2)) +
    geom_segment(aes(x = BP, xend = BP,
                     y = Y, yend = Yend),
                 linetype = "dashed", col = "grey60",
                 data = arcDat$label, inherit.aes = FALSE) +
    geom_curve(curvature = 1, ncp = 1000, lineend = 'butt') +
    geom_hline(yintercept = 0.5, col = "grey60") +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0.25, 0.75),
                       labels = c("JAM", "Stepwise")) +
    geom_label_repel(aes(x = BP, y = Y, label = SNP), data = arcDat$label,
                     inherit.aes = FALSE) +
    scale_color_gradient2(low = "grey90", mid = "yellow", high = "red",
                          limits = c(0, 1), midpoint = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank()) +
    theme(axis.line = element_blank(),
          #axis.title.x = element_blank(),
          axis.text.y = element_text(colour = "grey20"))
  
  #ggThemeLE() 
  #xlim(c(min(assoc$BP), max(assoc$BP)))
} # END plotArc

# plotSNPtype -------------------------------------------------------------
plotSNPtype <- function(data, opts){
  ggplot(data = data$assoc,
         aes(x = BP, xend = BP, 
             y = TYPED-1, yend = TYPED, 
             colour = ifelse(TYPED == 2, "#003333", "#669999"))) +
    geom_segment() +
    geom_hline(yintercept = 0:1 + 0.5,
               linetype = "dotted", col = "grey60") +
    scale_color_identity() +
    scale_y_continuous(breaks = c(0.5, 1.5),
                       labels = c("Imputed", "Typed"),
                       limits = c(-0.90, 2),
                       #labels=udf_pad(c("Imputed","Typed")),
                       name = "SNP") +  
    xlim(c(data$zoomStart, data$zoomEnd)) +
    theme_LE() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
}


# plotGene ----------------------------------------------------------------
plotGene <- function(chrom = NULL,
                     chromStart = NULL,
                     chromEnd = NULL,
                     vline = NULL){
  plotDatGeneN <- 1
  
  #get granges collapsed genes for ggplot+ggbio
  plotDatGene <- udf_GeneSymbol(chrom, chromStart, chromEnd)
  if(is.null(plotDatGene)){
    gg_out <- ggplot() + 
      geom_blank() +
      annotate("text",
               x = chromStart + (chromEnd - chromStart)/2,
               y = 0.35,
               label = "No gene") +
      #general options
      ylab("") + 
      xlim(c(chromStart, chromEnd)) +
      theme(legend.position = "none",
            panel.background = element_rect(fill="white"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.border = element_blank(),
            #Y Axis font
            axis.text = element_blank())
  } else {
    
    #number of genes in zoomed region, if no genes then 1
    plotDatGeneN <- try({
      length(unique(plotDatGene@elementMetadata$gene_id))}, silent = TRUE)
    if(class(plotDatGeneN) == "try-error"){ plotDatGeneN <- 1 }
    
    
    # return ggbio:gene plot
    gg_out <- ggplot() + 
      geom_hline(yintercept = c(1:plotDatGeneN),
                 col = "grey60", linetype = "dotted") +
      #mark hit SNPs
      # geom_vline(xintercept = plotDatStats %>% 
      #              filter(SNP %in% RegionHitsSelected) %>% .$BP,
      #            col="black",linetype="dashed") +
      #ggbio plot genes
      geom_alignment(data = plotDatGene,aes(group = gene_id,
                                            fill = strand, col = strand)) +
      #scale_y_continuous(#breaks = c(1:plotDatGeneN),
      #                  limits = c(0, plotDatGeneN)
      #labels=udf_pad(c("Imputed","Typed")),
      #name = "SNP"
      #                 ) +  
      #general options
      ylab("") + 
      xlim(c(chromStart, chromEnd)) +
      theme(legend.position="none",
            panel.background = element_rect(fill="white"),
            panel.grid.minor=element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey60", linetype = "dotted"),
            axis.title.x=element_blank(),
            #axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            axis.line = element_blank(),
            panel.border = element_blank(),
            #Y Axis font
            axis.text.y = element_text(#family = "Courier",
              colour = "grey20"))
    # mark hit SNPs - vline
    if(!is.null(vline)){
      gg_out <- gg_out +
        geom_vline(xintercept = vline,
                   col = "black", #col = "#4daf4a",
                   linetype = "dashed")
    }
  }
  #return output ggplot
  return(list(genePlot = gg_out,
              geneCnt = plotDatGeneN))
  
  
} # END plotGene


# plotAnnot ---------------------------------------------------------------
plotAnnot <- function(data,
                      chrom = NULL,
                      chromStart = NULL,
                      chromEnd = NULL,
                      vline = NULL,
                      collapse = FALSE){
  #subset data for zoomed region
  data <- data[ data$CHR == chrom &
                  data$START >= chromStart &
                  data$END <= chromEnd, ]
  
  if(collapse){ data[ FILE == "ChromHMM", TYPEN := 3] } 
  
  #plot
  gg_out <- 
    ggplot(data,
           aes(xmin = START, xmax = END, ymin = TYPEN - 1, ymax = TYPEN,
               col = NULL, fill = COLOUR_HEX )) +
    geom_rect(alpha = ifelse(collapse, 0.5, 1)) 
  
  if(collapse){
    gg_out <- gg_out +
      scale_y_continuous(
        limits = c(-1, 3),
        breaks = c(0:2) + 0.5, 
        labels = c("DNaseI", "Conserved", "ChromHMM"),
        name = "")
  } else {
    gg_out <- gg_out +
      scale_y_continuous(
        limits = c(-1, 11),
        breaks = c(0:10) + 0.5, 
        labels = c("DNaseI","Conserved",
                   #ChromHMM
                   "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
                   "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter"),
        name = "")
    
  }
  
  #prettify
  gg_out <- gg_out +
    scale_color_identity() +
    # general options
    xlim(c(chromStart, chromEnd)) +
    theme_LE() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  # mark hit SNPs - vline
  if(!is.null(vline)){
    gg_out <- gg_out +
      geom_vline(xintercept = vline,
                 col = "black", #col = "#4daf4a",
                 linetype = "dashed")
  }
  
  #return output ggplot
  return(gg_out)
} # END plotAnnot


# plotEQTL ----------------------------------------------------------------
plotEQTL <- function(data,
                     LD,
                     chromN = NULL,
                     chromStart = NULL,
                     chromEnd = NULL,
                     BP_hits = NULL,
                     BP_hits_tag = NULL
){
  
  # BP_hits = unique(mDat$LD[ SNP_A %in% mDat$hits, BP_A]),
  # BP_hits_tag = mDat$LD[ SNP_A %in% mDat$hits & R2 > 0.6, BP_B]
  # data = datEQTL
  # chromN = regionChrN
  # chromStart = regionStartZoom
  # chromEnd = regionEndZoom
  # BP_hits = mDat$LD[ SNP_A %in% mDat$hits & R2 > 0.6, BP_B]
  #subset eQTL TCGA for zoom
  plotDat <- data %>%
    filter(CHR == chromN &
             x %in% BP_hits_tag) %>% 
    mutate(col = if_else(x %in% BP_hits, "#e41a1c", "grey70"))
  
  #merge(plotDat, unique(LD[ ,list(BP_B, R2)]), by.x = "x", by.y = "BP_B")
  
  # if there is no eqtl then plot warning.
  if(nrow(plotDat) == 0){
    gg_out <- ggplot() + 
      geom_blank() +
      annotate("text",
               x = chromStart + (chromEnd - chromStart)/2,
               y = 0.35,
               label = "No eQTL")
  } else {
    # data is used to add Xaxis labels - Gene names, also to add
    #  vertical lines on gene plot to mark genes with a line.
    geneLabel <- plotDat %>%
      dplyr::transmute(x = xend, eQTL_TCGA = 0.35, label = GENE) %>%
      base::unique()
    
    gg_out <- 
      ggplot(data = plotDat,
             aes(x = x, xend = xend, y = eQTL_TCGA, yend = yend, col = col)) +
      geom_segment(linetype = "dashed") +
      #geom_segment(linetype = "dashed", col = "#e41a1c") +
      geom_text(aes(x = x, y = eQTL_TCGA, label = label),
                data = geneLabel, inherit.aes = FALSE, angle = 90,
                hjust = 1) +
      scale_color_identity()
  }
  
  #prettify
  gg_out <- gg_out +
    # scale_x_continuous(breaks = plotDatAnnotEQTL_GeneLabel()$x,
    #                    labels = plotDatAnnotEQTL_GeneLabel()$label) +
    scale_y_continuous(name = "",
                       breaks = c(0.35),
                       labels = c("eQTL TCGA"),
                       limits = c(-0.25, 1)) +
    xlim(c(chromStart, chromEnd)) +
    theme_LE() +
    theme(axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.x = element_blank())
  
  # return ggplot object
  return(gg_out)
  
} #END plotEQTL

# plotHistone -------------------------------------------------------------
plotHistone <- function(data = NULL,
                        chromStart = NULL,
                        chromEnd = NULL){
  
  ggplot(data,
         aes(BP, SCORE, fill = ENCODE)) +
    geom_area(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = wgEncodeBroadHistoneFileDesc$ColourDark) +
    scale_y_continuous(breaks = 50,
                       labels = "Histone",
                       limits = c(-10, 100),
                       name = "") +
    #general options
    xlim(c(chromStart, chromEnd)) +
    theme_LE() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
}

# GeneSymbol --------------------------------------------------------------
# https://github.com/oncogenetics/R_UDF/blob/master/GeneSymbol.R
udf_GeneSymbol <- function(chrom = NA, chromStart = NA, chromEnd = NA,
                           geneSymbolOnly = TRUE){
  # chrom = regionChr
  # chromStart = regionStartZoom
  # chromEnd = regionEndZoom
  # geneSymbolOnly = TRUE, collapse TX only if it matches with geneSymbol,
  #                  to avoid genes with no names.
  #                  FALSE, collapse TX on geneSymbol, if no match, then
  #                  use TX name as gene name.
  
  require(dplyr)
  require(ggbio)
  require(GenomicFeatures)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db) # gene symobols
  require(DBI)
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  #valid chrom names
  chr <- paste0("chr",c(1:22,"X","Y"))
  
  # Input checks ----------------------------------------------------------
  if(!is.character(chrom) |
     length(chrom)!=1 |
     !chrom %in% chr) stop("chrom: must be character class with length(chrom)==1, e.g.: chr1")
  
  if(is.na(chromStart) |
     chromStart < 0) warning("chromStart: setting to default value 0")
  
  if(is.na(chromEnd)) warning("chromEnd: default value last position for chromosome")
  
  # Validate start end position for chrom - using txdb chrominfo table
  chrominfo <- 
    dbGetQuery(txdb$conn,
               paste0("select * from chrominfo where chrom='",
                      chrom,"'"))
  chromStart <- ifelse(is.na(chromStart) | 
                         chromStart < 0,
                       0,chromStart)
  chromEnd <- ifelse(is.na(chromEnd) | 
                       chromEnd > chrominfo$length |
                       chromEnd < chromStart,
                     chrominfo$length, chromEnd)
  
  #Print summary for selection region of interest
  #print(paste0("Collapsing to gene symbols, region: ",
  #             chrom,":",chromStart,"-",chromEnd))
  
  # Get chromosome start end positions to subset TXDB for transcripts
  roi_chr <- GRanges(seqnames=chrominfo$chrom,
                     IRanges(start=chromStart,
                             end=chromEnd,
                             names=chrom))
  # Collapse to gene symbol -----------------------------------------------
  # Subset txdb over overlaps for chr-start-end
  keys_overlap_tx_id <- 
    subsetByOverlaps(transcripts(txdb), roi_chr) %>% 
    as.data.frame() %>% 
    .$tx_id %>% as.character
  
  # if there is a tx_id then collapse to geneSymbol else return NULL
  if(length(keys_overlap_tx_id) > 0){
    # Match TX ID to GENEID
    TXID_GENEID <- AnnotationDbi::select(txdb,
                                         keys = keys_overlap_tx_id,
                                         columns = c("TXNAME","GENEID"),
                                         keytype = "TXID")
    
    # Select transcipts from txdb which have GENEID
    Trans <- AnnotationDbi::select(txdb,
                                   keys = as.character(TXID_GENEID$TXID),
                                   columns = columns(txdb),
                                   keytype = "TXID")
    
    if(length(Trans[ !is.na(Trans$GENEID), "GENEID"]) > 0){
      
      # Get gene symbol
      gene_symbol <- unique(
        AnnotationDbi::select(org.Hs.eg.db,
                              #keys = Trans[ !is.na(Trans$GENEID), "GENEID"],
                              keys = Trans$GENEID,
                              columns = "SYMBOL",
                              keytype = "ENTREZID"))
      
      # Match GENEID, SYMBOL
      TXID_GENEID <- left_join(TXID_GENEID, gene_symbol,
                               by = c("GENEID" = "ENTREZID"))
      
      # # If not match on gene symbol, then TXNAME is gene symbol
      # TXID_GENEID$SYMBOL <- ifelse(is.na(TXID_GENEID$SYMBOL),
      #                              TXID_GENEID$TXNAME,TXID_GENEID$SYMBOL)
      # merge to add gene symbol
      Trans <- left_join(Trans,TXID_GENEID,
                         by=c("TXID","GENEID","TXNAME")) 
        
      # If not match on gene symbol, then TXNAME is gene symbol
      if(geneSymbolOnly){
        Trans <- Trans %>% filter(!is.na(SYMBOL))
        } else {
        Trans$SYMBOL <- ifelse(is.na(Trans$SYMBOL),
                               Trans$TXNAME,Trans$SYMBOL)
        }
      
      #Make Granges object 
      CollapsedGenes <- 
        GRanges(seqnames=Trans$EXONCHROM,
                IRanges(start=Trans$EXONSTART,
                        end=Trans$EXONEND),
                strand=Trans$EXONSTRAND)
      CollapsedGenes$gene_id <- Trans$SYMBOL
      
      # pad gene names to align
      # CollapsedGenes@elementMetadata$gene_id <-
      #   oncofunco::strPadLeft(labels = CollapsedGenes@elementMetadata$gene_id)
      
      # Output ----------------------------------------------------------------
      #return collapsed genes per CHR
      return(CollapsedGenes)
    } # END if(length(Trans[ !is.na(Trans$GENEID), "GENEID"]) > 0)
    # END if(length(keys_overlap_tx_id) > 0)
  } else {
    #no matching transcript/gene retunr null
    return(NULL)}
} # END udf_GeneSymbol



