# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

if(input$dataType %in% c("OncoArray", "OncoArrayMeta")){
  #if input data is OncoArray, then plot TCGA data eQTLs (Eze data)
  ggplot(data = plotDatAnnotEQTL(),
         aes(x = x, xend = xend, y = eQTL_TCGA, yend = yend)) +
    geom_segment(linetype = "dashed", col = "#e41a1c") +
    geom_text(aes(x = x, y = eQTL_TCGA, label = label),
              data = plotDatAnnotEQTL_GeneLabel(), inherit.aes = FALSE, angle = 90,
              hjust = 1) +
    # scale_x_continuous(breaks = plotDatAnnotEQTL_GeneLabel()$x,
    #                    labels = plotDatAnnotEQTL_GeneLabel()$label) +
    scale_y_continuous(breaks = c(0.2, 0.9),
                       labels = udf_pad(c("Gene", "SNP")),
                       limits = c(-0.25, 1)) +
    xlim(c(zoomStart(),zoomEnd())) +
    udf_theme() +
    theme(axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.x = element_blank())
  } else {
  #add barplot - Transcripts
  ggplot(ROIdatBedGraph(),
         aes(xmin=START,xmax=END,ymin=0,ymax=SCORE)) +
    geom_hline(aes(yintercept=0),linetype="dotted",col="grey60") +
    geom_rect(col="#A71930",fill="#A71930") +
    scale_y_continuous(limits=c(-1,1),
                       breaks=c(0),
                       labels=udf_pad(input$FileBedGraphName),
                       name="") +
    xlim(c(zoomStart(),zoomEnd())) +
    udf_theme() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
  
}