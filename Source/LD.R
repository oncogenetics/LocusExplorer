# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# LD track
# plot LD per hit SNPs on seperate Yaxis 1,2,3,etc

ggplot(data=plotDatLD(),
       aes(x=BP,xend=BP,
           y=as.numeric(as.factor(LDSNP)), yend=as.numeric(as.factor(LDSNP))+1,
           colour=LDCol)) +
  geom_hline(yintercept=1:length(RegionHitsSelected())+0.5,
             linetype="dotted",col="grey60") +
  geom_segment() +
  scale_colour_identity() +
  scale_y_continuous(breaks=(1:length(RegionHitsSelected()))+0.5,
                     labels=udf_pad(RegionHitsSelected()),
                    #name="xxx"
                    name=expression(R^2)
                    ) +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
