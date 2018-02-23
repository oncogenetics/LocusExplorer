# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt
plotList <-
  list(
    if("Manhattan" %in% input$ShowHideManhattanPvalues) plotObjManhattanPvalues() else NA,
    if("Manhattan" %in% input$ShowHideManhattanPostProbs) plotObjManhattanPostProbs() else NA,
    if("SNPType" %in% input$ShowHideTracks) plotObjSNPType() else NA,
    if("LD" %in% input$ShowHideTracks) plotObjSNPLD() else NA,
    if("wgEncodeBroadHistone" %in% input$ShowHideTracks) plotObjwgEncodeBroadHistone() else NA,
    if("annotOncoFinemap" %in% input$ShowHideTracks) plotObjAnnotOncoFinemap() else NA,
    if("Gene" %in% input$ShowHideTracks) plotObjGenePlot() else NA
  )

names(plotList) <- 
  c("ManhattanPvalues", "ManhattanPostProbs", "SNPType", "LD",
    "wgEncodeBroadHistone", "annotOncoFinemap", "Gene")

plotHeights <- c(400,
                 200,
                 10 * )
  1

plotList <- plotList[!is.na(plotList)]

tracks(
  plotList  
  #heights = c(7,3)
  )



# Old version -------------------------------------------------------------

 
# #make dynamic merge plot string
# 
# # exclude LDSmooth - it is part of Manhattan.R
# plotList <- 
#   input$ShowHideTracks[ !input$ShowHideTracks %in% 
#                           c("LDSmooth","Recombination")]
# 
# 
# 
# # character string to parse, plot only selected tracks
# trackString <-
#   paste(
#     "tracks(",
#     paste(
#       paste(sapply(plotList,function(i)
#         paste0("source('Source/",i,".R',local=TRUE)$value")
#       ),collapse=","),
#       "heights = trackHeights()",
#       "track.plot.color = trackColours()",
#       sep=","),
#     ", title = input$downloadPlotTitle)",
#     
#     # change Y axis text to defualt font
#     "+ theme(axis.text.y = element_text(family=''))",
#     
#     # Clear theme
#     ifelse(input$PlotTheme == "2",
#            "+ theme(panel.border = element_rect(fill = NA, color = 'grey80'))",
#            "")
#     )
# 
# #run plot string
# eval(parse(text=trackString))
# 
# 
# 
# # background grey fill optional
# #ifelse(input$PlotTrackShade == "1",
# #       "track.plot.color = trackColours()",
# #      "track.plot.color = NULL"),
