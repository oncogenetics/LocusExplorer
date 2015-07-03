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
  
  xlim(c(zoomStart(),zoomEnd())) +
  scale_y_continuous(
    limits=c(0,ROIPLogMax()),
    breaks=seq(0,ROIPLogMax(),5),
    labels=udf_pad(seq(0,ROIPLogMax(),5)),
    name="xxx"
    #name=expression(-log[10](italic(p)))
  ) +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  #general options
  scale_colour_identity() +
  theme(
    axis.title=element_blank(),
    #Y Axis font
    axis.text.y=element_text(family="Courier"),
    panel.background = element_rect(fill="white")
  )



# #ggplot(df1[ df1$LD > input$LD,],aes(x=BP,y=Pvalue)) + geom_point() +
# ggplot(datRegionAssocPlot(),aes(x=BP,y=PLog)) + 
#   #all snps grey hollow circles
#   geom_point(size=4,colour="grey80",shape=datRegionAssocPlot()$SNPType) +
#   geom_rug(sides="b") +
#   #LD colours per hit SNP
#   geom_point(data=datLDColour(),aes(BP_B,PLog),
#              size=4,shape=datLDColour()$SNPType+15,
#              colour=datLDColour()$LDCol,
#              alpha=0.8) +
#   xlim(roi()) +
#   scale_y_continuous(
#     limits=c(0,regionPLog_Max()),
#     breaks=seq(0,regionPLog_Max(),5),
#     labels=udf_pad(seq(0,regionPLog_Max(),5)),
#     name="pvalue"
#     #name=expression(-log[10](italic(p)))
#     ) +
#   scale_colour_identity() +
#   theme(
#     #Y Axis font
#     axis.text.y=element_text(family="mono"),
#     panel.background = element_rect(fill="white")
#     )



# gg_main_Ymax <- regionP_value_Log_Max()
# gg_main <- 
#   #geom_smooth per hit SNP - excluding R2 for the same SNP
#   ggplot(data=LD_1KG()[ SNP_A != SNP_B,],
#          aes(BP_B,R2_Adj,colour=LDSmoothCol)) +
#   #all snps grey hollow circles
#   geom_point(data=dat_region(),aes(position,P_value_Log),
#              size=4,shape=dat_region()$SNPType+1,
#              colour="grey80") +
#   #Recomb lines in background
#   geom_area(data=regionRecomb(),
#             aes(position,recomb_adj),
#             fill="deepskyblue",colour="deepskyblue1",alpha=0.2) +
#   #LD colours per hit SNP
#   geom_point(data=LD_1KG(),aes(BP_B,P_value_Log),
#              size=4,shape=LD_1KG()$SNPType+16,
#              colour=LD_1KG()$LDCol,
#              alpha=0.8) +
#   #LD smooth per hit SNP
#   geom_smooth(method="loess",se=FALSE) +
#   #mark hit SNPs
#   geom_segment(data=LD_1KG()[ HitSNP==1 &
#                                 !is.na(HitSNP),],
#                aes(x=BP_B, y=0, xend=BP_B, yend=P_value_Log),
#                colour="black",
#                linetype = "dashed") +
#   #Label hits
#   geom_text(data=dat_region()[ HitSNP == 1,],
#             aes(position,P_value_Log,label=snp_name),
#             vjust=1.1,hjust=-0.1,colour="black") +
#   #general options
#   scale_colour_identity() +
#   scale_y_continuous(limits=c(0,gg_main_Ymax),
#                      name=expression(-log[10](italic(p)))) + 
#   theme_gg_main()











#scale_y_discrete(breaks=c(1000,2000,3000),labels=c("__1K","__2K","__3K")) +
#   theme(plot.margin=unit(c(0,10,0,30),"mm"),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_text(family="mono"))


#str_pad(round(seq(0,10,3)), 10,"left","_")


#ylim(c(0,regionPLog_Max())) +

#labs(x=NULL) +
#   scale_y_continuous(limits=c(0,regionPLog_Max()),
#                      name=expression(-log[10](italic(p)))) 
#scale_y_discrete(
