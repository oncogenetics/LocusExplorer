ggplot() + 
  geom_hline(yintercept=c(1:plotDatGeneN()),col="grey80",linetype="dashed") +
  #ggbio
  geom_alignment(data = plotDatGene(),aes(group=gene_id,
                     fill=strand, col=strand)) +
  #facet_grid(strand~.) +
  #general options
  ylab("xxx") + 
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red")
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() 
  