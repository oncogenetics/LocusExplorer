# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

gg_out <-
  ggplot(plotDatStats(), aes(x=BP,y=PLog)) + 
  # all snps grey hollow shapes
  geom_point(size=4,colour="grey80",shape=plotDatStats()$TYPED)


if("Recombination" %in% input$ShowHideTracks){ 
  #Recomb lines
  gg_out <- gg_out +
    geom_area(data=plotDatGeneticMap(),
              aes(BP,RecombAdj),
              fill="#11d0ff",colour="#11d0ff",alpha=0.2)}

if("LD" %in% input$ShowHideTracks){ 
    # LD colours - add R2 filled shapes
    gg_out <- gg_out +
      geom_point(data=plotDatLD(),aes(x=BP,y=PLog),
                 size=4,
                 shape=ifelse(plotDatLD()$TYPED==2,17,16),
                 col=plotDatLD()$LDCol,
                 alpha=0.8)}

# LD smooth per hit SNP - optional
if("LDSmooth" %in% input$ShowHideTracks) {
  gg_out <- gg_out +
    geom_smooth(data=plotDatLD(),aes(x=BP,y=R2_Adj,col=LDSmoothCol),
                method="loess",se=FALSE)}

# add other plots
gg_out <- gg_out +
  #mark hit SNPs - outline shapes
  geom_point(data=plotDatStats() %>% 
               filter(SNP %in% input$HitSNPs),
             aes(x=BP,y=PLog),size=4,colour="black",
             shape=plotDatStats() %>% 
               filter(SNP %in% input$HitSNPs) %>% 
               .$TYPED) +
  #mark hit SNPs - vertical lines
  geom_segment(data=plotDatLD() %>% 
                 filter(SNP==LDSNP),
               aes(x=BP, y=0, xend=BP, yend=PLog),
               colour="black",
               linetype = "dashed") +
  #mark hit SNPs - Label hits
  geom_text(data=plotDatLD() %>% 
              filter(SNP==LDSNP),
            aes(BP,PLog,label=LDSNP),
            vjust=1.1,hjust=-0.1,colour="black") +
  #zoom
  xlim(c(zoomStart(),zoomEnd())) +
  scale_y_continuous(
    limits=c(0,ROIPLogMax()),
    breaks=seq(0,ROIPLogMax(),5),
    labels=udf_pad(seq(0,ROIPLogMax(),5)),
    #name="xxx"
    name=expression(-log[10](italic(p)))
  ) +
  scale_colour_identity() +
  theme(
    #axis.title=element_blank(),
    #Y Axis font
    axis.text.y=element_text(family="Courier"),
    panel.background = element_rect(fill="white")
  ) +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  udf_theme() 


#result
gg_out
