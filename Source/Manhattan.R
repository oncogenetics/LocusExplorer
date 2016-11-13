# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt


g <- ggplot(plotDatAssoc(), aes(x=BP,y=PLog)) + 
  # all snps grey hollow shapes
  geom_point(size=4,colour="#B8B8B8",shape=plotDatAssoc()$TYPED)

g

# # Plot all SNPs - background ----------------------------------------------
# gg_out <-
#   if(!is.null(plotDatAssoc())){
#     ggplot(plotDatAssoc(), aes(x=BP,y=PLog)) + 
#       # all snps grey hollow shapes
#       geom_point(size=4,colour="#B8B8B8",shape=plotDatAssoc()$TYPED)
#     } else {NULL}
# 
# # Recombination -----------------------------------------------------------
# if("Recombination" %in% input$ShowHideTracks &
#    nrow(plotDatGeneticMap()) > 2 ){ 
#   #Recomb lines
#   gg_out <- gg_out +
#     geom_area(data=plotDatGeneticMap(),
#               aes(BP,RecombAdj),
#               fill="#11d0ff",colour="#00B4E0",alpha=0.3)}
# 
# # LD ----------------------------------------------------------------------
# if("LD" %in% input$ShowHideTracks &
#    nrow(plotDatLD()) > 0){ 
#     # LD colours - add R2 filled shapes
#     gg_out <- gg_out +
#       geom_point(data=plotDatLD(),aes(x=BP,y=PLog),
#                  size=4,
#                  shape=ifelse(plotDatLD()$TYPED==2,17,16),
#                  col=plotDatLD()$LDCol,
#                  alpha=0.8)}
# 
# # LD Smooth ---------------------------------------------------------------
# # LD smooth per hit SNP - optional
# if("LDSmooth" %in% input$ShowHideTracks &
#    nrow(plotDatLD()) > 0) {
#   gg_out <- gg_out +
#     geom_smooth(data=plotDatLD(),aes(x=BP,y=R2_Adj,col=LDSmoothCol),
#                 method="loess",se=FALSE)}
# 
# # Suggestiveline and Genomewideline ---------------------------------------
# # Do not plot if set to zero
# 
# # suggestiveLine
# if(input$suggestiveLine != 0) {
#   gg_out <- gg_out +
#     geom_hline(aes(yintercept = y), data = data.frame(y=input$suggestiveLine),
#                #linetype = "dashed",
#                size = 0.5,
#                colour = "#1a9641")}
# # genomewideLine
# if(input$genomewideLine != 0) {
#   gg_out <- gg_out +
#     geom_hline(aes(yintercept = y), data = data.frame(y=input$genomewideLine),
#                #linetype = "dashed",
#                size = 0.5,
#                colour = "#ca0020")}
# 
# # Mark Hits, shape and vline ----------------------------------------------
# gg_out <- gg_out +
#   #mark hit SNPs - outline shapes
#   geom_point(data=plotDatAssoc() %>% 
#                filter(SNP %in% RegionHitsSelected()),
#              aes(x=BP,y=PLog),size=4,colour="black",
#              shape=plotDatAssoc() %>% 
#                filter(SNP %in% RegionHitsSelected()) %>% 
#                .$TYPED) +
#   #mark hit SNPs - vertical lines
#   geom_segment(data=plotDatAssoc() %>% 
#                  filter(SNP %in% RegionHitsSelected()),
#                aes(x=BP, y=0, xend=BP, yend=PLog),
#                colour="#4daf4a",
#                linetype = "dashed")
# 
# # Mark Hits, Labels -------------------------------------------------------
# #mark hit SNPs - Label hits - use repulsion force to push labels away
# # this needs more work, repulsion force should be calcualted automatically based
# # on size of the region... currently input is a slider "input$adjustLabels"
# hitLabels <- plotDatAssoc() %>% 
#   filter(SNP %in% RegionHitsSelected())
# 
# if(input$adjustLabels & nrow(hitLabels) > 1) {
#   x.fact <- 1000/max(hitLabels$BP)
#   y.fact <- 100/max(hitLabels$PLog)
#   coords <- FFieldPtRep(coords = cbind(hitLabels$BP * x.fact, 
#                                        hitLabels$PLog * y.fact),
#                         rep.fact = input$repFact, 
#                         rep.dist.lmt = 5, 
#                         attr.fact = 0.2, 
#                         adj.max = 0.1, 
#                         adj.lmt = 0.5,
#                         iter.max = 10000)
#   
#   hitLabels$x.t <- coords$x/x.fact
#   hitLabels$y.t <- coords$y/y.fact
#   
#   gg_out <- gg_out +
#     geom_text(data=hitLabels,
#               aes(x.t,y.t,label=SNP),
#               vjust=-1, colour="black") +
#     geom_segment(data=hitLabels,aes(xend = x.t, yend = y.t)) 
#   } else {
#   gg_out <- gg_out +
#        geom_text(data=hitLabels,
#                  aes(BP,PLog,label=SNP),
#                  vjust=1.1,hjust=-0.1) 
#     }
# 
# 
# # Mark TCGA eQTL SNPs -----------------------------------------------------
# # Data is only available for OncoArray
# if(input$dataType == "OncoArray" &
#    "BedGraph" %in% input$ShowHideTracks){
#   
#   dd <- ROIdatStats() %>% filter(BP  %in% plotDatAnnotEQTL()$x)
#   
#   gg_out <- gg_out +
#     #mark hit SNPs - outline shapes
#     geom_point(data = dd,
#                aes(x = BP, y = PLog), size = 4, colour = "#e41a1c",
#                shape = 23 #ifelse(dd$TYPED == 2, 17, 16)
#                ) +
#     #mark hit SNPs - vertical lines
#     geom_segment(data = dd,
#                  aes(x = BP, y = 0, xend = BP, yend = PLog),
#                  colour = "#e41a1c",
#                  linetype = "dashed")
#   
#   #add SNP names on manhattan for eQTL SNPs in red
#   
#   gg_out <- gg_out +
#     geom_text(data = dd,
#               aes(x = BP, y = PLog, label = SNP),
#               vjust =-1, colour = "#e41a1c") 
#   
#   
#   # hitLabels <- plotDatAssoc() %>% 
#   #   filter(SNP %in% RegionHitsSelected())
#   # 
#   # if(input$adjustLabels & nrow(hitLabels) > 1) {
#   #   x.fact <- 1000/max(hitLabels$BP)
#   #   y.fact <- 100/max(hitLabels$PLog)
#   #   coords <- FFieldPtRep(coords = cbind(hitLabels$BP * x.fact, 
#   #                                        hitLabels$PLog * y.fact),
#   #                         rep.fact = input$repFact, 
#   #                         rep.dist.lmt = 5, 
#   #                         attr.fact = 0.2, 
#   #                         adj.max = 0.1, 
#   #                         adj.lmt = 0.5,
#   #                         iter.max = 10000)
#   #   
#   #   hitLabels$x.t <- coords$x/x.fact
#   #   hitLabels$y.t <- coords$y/y.fact
#   
#   } # END if(input$dataType == "OncoArray" ...
# 
# # General options ---------------------------------------------------------
# # Zoom
# gg_out <- gg_out +
#   xlim(c(zoomStart(),zoomEnd())) +
#   scale_y_continuous(
#     limits=c(0,ROIPLogMax()),
#     breaks=seq(0,ROIPLogMax(),5),
#     labels=udf_pad(seq(0,ROIPLogMax(),5)),
#     name=expression(-log[10](italic(p)))) +
#   scale_colour_identity() +
#   theme(
#     #Y Axis font
#     axis.text.y=element_text(family="Courier"),
#     panel.background = element_rect(fill="white")) +
#   udf_theme()
# 
# # Output ------------------------------------------------------------------
# gg_out
