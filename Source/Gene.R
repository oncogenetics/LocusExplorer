# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

gg_out <- 
  ggplot() + 
  geom_hline(yintercept=c(1:plotDatGeneN()),col="grey60",linetype="dotted") +
  #mark hit SNPs
  geom_vline(xintercept=plotDatStats() %>% 
               filter(SNP %in% RegionHitsSelected()) %>% .$BP,
             col="#4daf4a", linetype="dashed") +
  #ggbio plot genes
  geom_alignment(data = plotDatGene(),aes(group=gene_id,
                                          fill=strand, col=strand)) +
  #facet_grid(strand~.) +
  #general options
  ylab("") + 
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() 

#annotate vertical line eQTL over genes
if("BedGraph" %in% input$ShowHideTracks){ 
  gg_out <- gg_out +
    geom_vline(xintercept = ROIdatAnnotEQTL_GeneLabel()$x,
               col = "#e41a1c", linetype = "dashed")
    }

#return
gg_out