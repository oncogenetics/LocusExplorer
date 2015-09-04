#make dynamic merge plot string

# exclude LDSmooth - it is part of Manhattan.R
plotList <- 
  input$ShowHideTracks[ !input$ShowHideTracks %in% 
                          c("LDSmooth","Recombination")]

# character string to parse, plot only selected tracks
trackString <-
  paste("tracks(",
        paste(
          paste(sapply(plotList,function(i)
            paste0("source('Source/",i,".R',local=TRUE)$value")
          ),collapse=","),
          "heights = trackHeights()",
          "track.plot.color = trackColours()",
          sep=","),",title = input$downloadPlotTitle)")

#run plot string
eval(parse(text=trackString))
