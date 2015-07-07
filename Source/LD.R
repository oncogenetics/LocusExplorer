# LD track
# plot LD per hit SNPs on seperate Yaxis 1,2,3,etc


ggplot(data=plotDatManhattan(),
       aes(x=BP,xend=BP,
           y=as.numeric(as.factor(LDSNP)), yend=as.numeric(as.factor(LDSNP))+1,
           colour=LDSmoothCol)) +
  geom_hline(yintercept=1:length(RegionHitsSelected())+0.5,
             linetype="dashed",col="grey80") +
  geom_segment() +
  scale_colour_identity() +
  scale_y_continuous(breaks=(1:length(RegionHitsSelected()))+0.5,
                     labels=udf_pad(RegionHitsSelected()),
                    name="xxx") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()
