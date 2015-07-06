ggplot(data=ROIdatStats(),
       aes(x=BP, xend=BP, 
           y=TYPED-1,yend=TYPED, 
           colour=ifelse(TYPED==2,"#616365","#ADAFAF"))) +
  geom_segment() +
  scale_color_identity() +
  scale_y_continuous(breaks=c(0.5,1.5),
                     labels=udf_pad(c("Imputed","Typed")),
                     name="xxx") +  
  #testing plot alignment
  geom_vline(xintercept=173000000,col="red") +
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()

# 
#     theme(legend.position="none",
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
