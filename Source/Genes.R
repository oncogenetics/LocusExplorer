ggplot() + 
  geom_hline(yintercept=c(1:plotDatGeneN()),col="grey80",linetype="dashed") +
  #ggbio
  geom_alignment(data = plotDatGene(),aes(group=gene_id,
                     fill=strand, col=strand)) +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  ylab("xxx") + 
  udf_theme() +
  #testing plot alignment
  geom_vline(xintercept=173000000,col="red")
  
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

# #subset region from list Grange object - "plotGenes"
# 
# #subset genes for roi
# regionGenes <- subsetByOverlaps(plotGenes[[as.character(seqnames(roi_gene()))]],
#                                 roi_gene())
# 
# #check if there are any genes in roi
# if(mcols(regionGenes) %>% nrow == 0){
#   #no genes to plot - plot text - "No Genes"
#   ggplot(data=data.frame(x=c(input$Xrange[1],input$Xrange[2]),y=1),
#          aes(x=x, y=y)) +
#     scale_y_continuous(breaks=NULL) +
#     geom_blank() +
#     annotate("text",
#              #put text in the middle
#              x=input$Xrange[1]+(input$Xrange[2]-input$Xrange[1])/2,
#              y=1,label = "No Genes")
# }else{
#   #pad gene names
#   regionGenes$gene_id <- udf_pad(regionGenes$gene_id)
#   ggplot(regionGenes) + 
#     geom_alignment(aes(group=gene_id,colour=strand,fill=strand), 
#                    ylab="ddd") +
#     xlim(roi()) + 
#     theme(legend.position="none",
#           #Y Axis font
#           axis.text.y=element_text(family="mono"),
#           panel.grid.major.x=element_blank(),
#           panel.grid.minor.x=element_blank(),
#           panel.border=element_blank())
# } #if
# 
