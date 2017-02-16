#' LocusExplorer - Manhattan plot
#'
#' Manhattan plot for LocusExplorer.
#' @param assoc SNP association results, data.frame object with c("SNP","BP","P") columns. Required.
#' @param LD plink LD output format, data.frame object with c("BP_A","SNP_A","BP_B","SNP_B","R2") columns. Optional/recommended.
#' @param geneticMap Recombination map, data.frame object with c("BP", "RECOMB") columns. Subset of one of genetic_map_*_combined_b37.txt, at http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/ . Optional.
#' @param suggestiveLine Suggestive line, default is 5.
#' @param genomewideLine Genomewide link, dafault is 8.
#' @param xStart,xEnd Region range, zoom, minimum BP and maximum BP, advised to keep this less than 5Mb.
#' @param hits SNP names to label in the plot. Must be present in assoc data.frame.
#' @param hitsName alternative SNP names to label in the plot. Default same as `hits`
#' @param hitsLabel Default is TRUE, set to FALSE not to show SNP names on the plot.
#' @param pad Default is TRUE, to align plots pad strings with spaces, using oncofunco::strPadLeft().
#' @param title character string for plot title. Default is NULL, i.e.: no plot title. 
#' @param opts Default is c("Recombination","LD","LDSmooth","SuggestiveLine","GenomewideLine","Hits"), parts of plot to display.
#' @export plotManhattan
#' @author Tokhir Dadaev
#' @return a \code{ggplot} object
#' @keywords manhattan plot SNP genetics


plotManhattan <- function(
  assoc = NULL,
  LD = NULL,
  geneticMap = NULL,
  suggestiveLine = 5,
  genomewideLine = 8,
  xStart = NULL,
  xEnd = NULL,
  hits = NULL,
  hitsName = hits,
  hitsLabel = TRUE,
  pad = TRUE,
  title = NULL,
  opts = c("Recombination","LD","LDSmooth","SuggestiveLine",
           "GenomewideLine","Hits")){
  
  # Check input - assoc -------------------------------------------------------
  #check assoc
  if(is.null(assoc)) stop("assoc is missing, must provide assoc with columns: c('SNP','BP','P')")
  if(!all(c("SNP","BP","P") %in% colnames(assoc))) stop("assoc must have columns: c('SNP','BP','P')")
  # if SNP type is missing set all as typed
  if(!"TYPED" %in% colnames(assoc)){assoc$TYPED <- 2}
  assoc <- as.data.frame(assoc)
  #set plog max
  assoc$PLog <- -log10(assoc$P)
  
  # XY range ------------------------------------------------------------------
  if(is.null(xStart))xStart <- min(assoc$BP, na.rm = TRUE)
  if(is.null(xEnd))xEnd <- max(assoc$BP, na.rm = TRUE)
  yMax <- ceiling(max(c(10, assoc$PLog)))
  yRange <- c(0, max(c(10, ceiling((yMax + 1)/5) * 5)))
  xRange <- c(xStart, xEnd)
  
  #Check input - recomb -------------------------------------------------------
  if("Recombination" %in% opts){
    if(is.null(geneticMap)) stop("geneticMap data is missing for recombination, must have columns: c('BP', 'RECOMB')")
    if(!all(c("BP", "RECOMB") %in% colnames(geneticMap))) stop("geneticMap data must have columns: c('BP', 'RECOMB')")
    geneticMap <- data.frame(
      BP = geneticMap$BP,
      #adjust recomb value to pvalue
      RECOMB_ADJ = geneticMap$RECOMB * yMax / 100)}
  
  # Plot all SNPs - background ------------------------------------------------
  gg_out <-
    ggplot(assoc, aes(x = BP, y = PLog)) +
    # all snps grey hollow shapes
    geom_point(size = 4, colour = "#B8B8B8", shape = assoc$TYPED) +
    geom_hline(yintercept = seq(0, yMax, 5),
               linetype = "dotted", col = "grey60")
  
  # Plot - Recombination ------------------------------------------------------
  if("Recombination" %in% opts & nrow(geneticMap) > 2 ){
    gg_out <- gg_out +
      geom_area(data = geneticMap,
                aes(BP, RECOMB_ADJ),
                fill = "#11d0ff", colour = "#00B4E0", alpha = 0.3)}
  
  
  # Check input - LD ----------------------------------------------------------
  if("LD" %in% opts | "LDSmooth" %in% opts){
    if(is.null(LD)) stop("LD is missing, must have columns: c('BP_A','SNP_A','BP_B','SNP_B','R2')")
    if(!all(c("BP_A","SNP_A","BP_B","SNP_B","R2") %in% colnames(LD)))
      stop("LD must have columns: c('BP_A','SNP_A','BP_B','SNP_B','R2')")
    
    LD <- as.data.frame(LD)
    
    if(is.null(hits)){
      hits <- unique(LD$SNP_A)
      hits <- hits[1:min(5, length(hits))]
      warning(
        paste("hits missing, selected first <5 SNPs as hits from LD$SNP_A, n = :",
              length(unique(LD$SNP_A))))
    } else {
      hits <- sort(intersect(hits, LD$SNP_A))
      
      colourLD <- oncofunco::colourHue(length(hits))
      colourLDPalette <- unlist(lapply(colourLD, function(i){
        colorRampPalette(c("grey95", i))(100)}))
      
      #merge LD with assoc, to get R2 shades per point
      plotDat <- merge(
        LD[ LD$SNP_A %in% hits, c("BP_A","SNP_A","BP_B","SNP_B","R2")],
        assoc[, c("BP", "TYPED", "PLog")],
        by.x = "BP_B", by.y = "BP") %>% 
        mutate(
          LDColIndex = ifelse(round(R2,2) == 0, 1, round(R2, 2) * 100),
          hitColIndex = as.numeric(factor(SNP_A, levels = hits)),
          hitCol = colourLD[hitColIndex],
          LDCol = colourLDPalette[(hitColIndex - 1) * 100 + LDColIndex],
          R2Adj = yMax * R2 * 0.8)
      # Plot - LD Fill & LD Smooth --------------------------------------------
      #LD fill
      if("LD" %in% opts){
        gg_out <- gg_out +
          geom_point(data = plotDat, aes(BP_B, PLog, colour = LDCol),
                     size = 4,
                     shape = plotDat$TYPED + 15,
                     col = plotDat$LDCol,
                     alpha = 0.8)
      }
      #LDSmooth
      if("LDSmooth" %in% opts){
        gg_out <- gg_out +
          geom_smooth(data = plotDat, aes(x = BP_B, y = R2Adj, col = hitCol),
                      method = "loess", se = FALSE)
      }
    }
  } # END if("LD" %in% opts | "LDSmooth" %in% opts)
  
  
  # Suggestiveline ----------------------------------------------------------
  if("SuggestiveLine" %in% opts & 
     !is.null(suggestiveLine) &
     suggestiveLine > 0){
    gg_out <- gg_out +
      geom_hline(aes(yintercept = y), data = data.frame(y = suggestiveLine),
                 size = 0.5,
                 colour = "#1a9641")}
  # Genomewideline ----------------------------------------------------------
  if("GenomewideLine" %in% opts &
     !is.null(genomewideLine) &
     genomewideLine > 0){
    gg_out <- gg_out +
      geom_hline(aes(yintercept = y), data = data.frame(y = genomewideLine),
                 size = 0.5,
                 colour = "#ca0020")}
  
  # Mark Hits: shape and vline ----------------------------------------------
  if("Hits" %in% opts & !is.null(hits) & any(hits %in% assoc$SNP)){
    gg_out <- gg_out +
      #mark hit SNPs - outline shapes
      geom_point(data = assoc[ assoc$SNP %in% hits, ],
                 aes(x = BP, y = PLog, shape = TYPED),
                 size = 4, colour = "black") +
      scale_shape_identity() +
      #mark hit SNPs - vertical lines
      geom_segment(data = assoc[ assoc$SNP %in% hits, ],
                   aes(x = BP, y = 0, xend = BP, yend = PLog),
                   colour = "black",
                   linetype = "dashed")}
  
  
  # Mark Hits: Labels -------------------------------------------------------
  # SNP names on the plot for hits,
  # if alternative names given then use those, hitsName
  if("Hits" %in% opts & length(hits) > 0)
    if(!is.null(hitsLabel))
      if(hitsLabel){
        plotDat <- assoc[ assoc$SNP %in% hits, ]
        
        if(all(hits == hitsName)) {
          plotDat$label <- plotDat$SNP
        } else {
          plotDat$label <- hitsName[match(plotDat$SNP, hits)]
        }
        
        
        
        gg_out <- 
          gg_out +
          geom_text_repel(
            aes(BP, PLog, label = label),
            data = plotDat)
        
      }
  
  # Add title ---------------------------------------------------------------
  if(!is.null(title)) gg_out <- gg_out + ggtitle(title)
  
  # General options ---------------------------------------------------------
  gg_out <- gg_out +
    coord_cartesian(
      xlim = xRange,
      ylim = yRange) +
    scale_y_continuous(
      breaks = seq(0, yMax, 5),
      #labels = oncofunco::strPadLeft(seq(0, ROIPLogMax, 5)),
      labels = if(pad){strPadLeft(seq(0, yMax, 5))} else {
        seq(0, yMax, 5)},
      name = expression(-log[10](italic(p)))) +
    scale_colour_identity()
  
  # Output ------------------------------------------------------------------
  gg_out
} #END plotManhattan



