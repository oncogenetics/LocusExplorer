# ENCODE wgEncodeRegDnaseClustered

ggplot(ROIdatwgEncodeRegDnaseClustered(),
       aes(xmin=START,xmax=END,ymin=0,ymax=1,fill=SCORE)) +
  geom_rect(alpha=0.5) +
  scale_fill_continuous(limits=c(500,1000),low="white",high="#A71930",na.value = "white") +
  scale_y_continuous(breaks=0.5,
                     labels=udf_pad("Dnase"),
                     name="") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()
