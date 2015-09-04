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
                     name="") +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()
