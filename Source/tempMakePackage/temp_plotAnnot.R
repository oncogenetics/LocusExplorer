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
  
  #table(data$TYPEN)
  #collapse=TRUE
  
  #plot
  gg_out <- 
    ggplot(data,
           aes(xmin = START, xmax = END, ymin = TYPEN - 1, ymax = TYPEN,
               col = NULL, fill = COLOUR_HEX )) +
    geom_rect(alpha = ifelse(collapse, 0.5, 1)) 
  
  if(collapse){
    gg_out <- gg_out +
      scale_y_continuous(
        limits = c(-1, 4),
        breaks = c(0:3) + 0.5, 
        labels = c("DNaseI", "Conserved", "ChromHMM", "eQTL"),
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
