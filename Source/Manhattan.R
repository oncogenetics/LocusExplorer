gg_out <-
  ggplot(plotDatManhattan(), aes(x=BP,y=PLog)) + 
  # all snps grey hollow circles
  geom_point(size=4,colour="grey80",shape=plotDatManhattan()$TYPED) +
  #Recomb lines in background
  geom_area(data=plotDatGeneticMap(),
            aes(BP,RecombAdj),
            fill="#11d0ff",colour="#11d0ff",alpha=0.2) +
  # add R2 filled shapes
  geom_point(data=plotDatManhattan(),aes(x=BP,y=PLog),
             size=4,
             shape=ifelse(plotDatManhattan()$TYPED==2,17,16),
             col=plotDatManhattan()$LDCol,
             alpha=0.8)

# LD smooth per hit SNP - optional
if("LDSmooth" %in% input$ShowHideTracks) {
  gg_out <- gg_out +
    geom_smooth(data=plotDatManhattan(),aes(x=BP,y=R2_Adj,col=LDSmoothCol),
                method="loess",se=FALSE)}

# add other plots
gg_out <- gg_out +
  #mark hit SNPs
  geom_segment(data=plotDatManhattan() %>% 
                 filter(SNP==LDSNP),
               aes(x=BP, y=0, xend=BP, yend=PLog),
               colour="black",
               linetype = "dashed") +
  #Label hits
  geom_text(data=plotDatManhattan() %>% 
              filter(SNP==LDSNP),
            aes(BP,PLog,label=LDSNP),
            vjust=1.1,hjust=-0.1,colour="black") +
  #zoom
  xlim(c(zoomStart(),zoomEnd())) +
  scale_y_continuous(
    limits=c(0,ROIPLogMax()),
    breaks=seq(0,ROIPLogMax(),5),
    labels=udf_pad(seq(0,ROIPLogMax(),5)),
    name="xxx"
    #name=expression(-log[10](italic(p)))
  ) +
  scale_colour_identity() +
  theme(
    axis.title=element_blank(),
    #Y Axis font
    axis.text.y=element_text(family="Courier"),
    panel.background = element_rect(fill="white")
  ) +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  udf_theme() 


#result
gg_out
