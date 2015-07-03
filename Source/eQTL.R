#add barplot - Transcripts
ggplot(ROIdatEQTL(),
       aes(xmin=START,xmax=END,ymin=0,ymax=DIRECTION)) +
  geom_hline(aes(yintercept=0),linetype="dashed",col="grey80") +
  geom_rect(col="#6E273D",fill="#6E273D") +
  scale_y_continuous(limits=c(-1,1),
                     breaks=c(0),
                     labels=udf_pad("eQTL"),
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
