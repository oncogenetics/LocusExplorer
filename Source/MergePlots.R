#make dynamic merge plot string
trackString <-
  paste("tracks(",
        paste(
          paste(sapply(input$ShowHideTracks,function(i)
            paste0("source('Source/",i,".R',local=TRUE)$value")
          ),collapse=","),
          "heights = trackHeights()",
          "track.plot.color = trackColours()",
          sep=","),",title = input$downloadPlotTitle)")

#run plot string
eval(parse(text=trackString))
