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
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  theme(legend.position="none",
        panel.background = element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        #Y Axis font
        axis.text.y=element_text(family="Courier")
  )
