#subset region - heatmap style track - density of region.
#make density range
temp_bw <- (RegionEnd()-RegionStart())/100
temp <- density(ROIdatLNCAP()$BP, bw=temp_bw, n=temp_bw)
temp <- data.table(BP=round(temp$x,0),LNCAP=temp$y)

#plot
ggplot(data=temp,
       aes(x=BP,xend=BP,y=0,yend=1,colour=LNCAP)) +
  geom_segment() +
  scale_color_continuous(low="white",high="#A71930") +
  scale_y_continuous(limits=c(0,1),breaks=c(0.5), 
                     labels=udf_pad("LNCAP"),
                     name="xxx") +
  #testing plot alignment
  geom_vline(xintercept=173000000,col="red") +
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()

#   theme(legend.position="none",
#         panel.background = element_rect(fill="white"),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         axis.title=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line=element_blank(),
#         panel.border=element_blank(),
#         #Y Axis font
#         axis.text.y=element_text(family="Courier")
#   )



