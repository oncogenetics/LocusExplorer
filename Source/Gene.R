# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

ggplot() + 
  geom_hline(yintercept=c(1:plotDatGeneN()),col="grey80",linetype="dashed") +
  #mark hit SNPs
  geom_vline(xintercept=plotDatStats() %>% 
               filter(SNP %in% RegionHitsSelected()) %>% .$BP,
             col="black",linetype="dashed") +
  #ggbio plot genes
  geom_alignment(data = plotDatGene(),aes(group=gene_id,
                     fill=strand, col=strand)) +
  #facet_grid(strand~.) +
  #general options
  ylab("") + 
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() 
