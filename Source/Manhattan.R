ggplot(plotDatManhattan(), aes(x=BP,y=PLog)) + 
  # all snps grey hollow circles
  geom_point(size=4,colour="grey80",shape=plotDatManhattan()$TYPED) +
  #Recomb lines in background
  geom_area(data=plotDatGeneticMap(),
            aes(BP,RecombAdj),
            #fill="deepskyblue",colour="deepskyblue1",alpha=0.2) +
            fill="#11d0ff",colour="#11d0ff",alpha=0.2) +
  # add R2 filled shapes
  geom_point(data=plotDatManhattan(),aes(x=BP,y=PLog),
             size=4,
             shape=ifelse(plotDatManhattan()$TYPED==2,17,16),
             #fill=plotDatManhattan()$LDCol,
             col=plotDatManhattan()$LDCol,#"grey80",
             alpha=0.8) +
  #LD smooth per hit SNP
  geom_smooth(data=plotDatManhattan(),aes(x=BP,y=R2_Adj,col=LDSmoothCol),
              method="loess",se=FALSE) +
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
    axis.text.y=element_text(family="mono"),
    panel.background = element_rect(fill="white")
  ) 

