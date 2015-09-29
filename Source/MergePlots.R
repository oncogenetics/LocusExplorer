# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

#make dynamic merge plot string

# exclude LDSmooth - it is part of Manhattan.R
plotList <- 
  input$ShowHideTracks[ !input$ShowHideTracks %in% 
                          c("LDSmooth","Recombination")]



# character string to parse, plot only selected tracks
trackString <-
  paste(
    "tracks(",
    paste(
      paste(sapply(plotList,function(i)
        paste0("source('Source/",i,".R',local=TRUE)$value")
      ),collapse=","),
      "heights = trackHeights()",
      "track.plot.color = trackColours()",
      sep=","),
    ", title = input$downloadPlotTitle)",
    
    # Clear theme
    ifelse(input$PlotTheme == "4",
           "+theme(panel.border = element_rect(fill = NA, color = 'grey80'))",
           "")
    )

#run plot string
eval(parse(text=trackString))



# background grey fill optional
#ifelse(input$PlotTrackShade == "1",
#       "track.plot.color = trackColours()",
#      "track.plot.color = NULL"),
