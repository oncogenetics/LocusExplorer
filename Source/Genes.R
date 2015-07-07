ggplot() + 
  geom_hline(yintercept=c(1:plotDatGeneN()),col="grey80",linetype="dashed") +
  #ggbio
  geom_alignment(data = plotDatGene(),aes(group=gene_id,
                     fill=strand, col=strand)) +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  ylab("xxx") + 
  udf_theme() 

