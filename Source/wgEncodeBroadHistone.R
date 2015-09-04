# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# ENCODE bigwig data 
# 7 bigwig filtered data, score > 5, similar to UCSC track, overlaying area plots

# if bigWig files are missing plot warning.
if(nrow(plotDatwgEncodeBroadHistone())==0){
  gg_out <- ggplot() + 
    geom_blank() +
    annotate("text",
             x=zoomStart()+(zoomEnd()-zoomStart())/2,
             y=50,label="Please download Histone bigWig files!")
} else {
  gg_out <- ggplot(data=plotDatwgEncodeBroadHistone(),
                   aes(BP,SCORE,fill=ENCODE)) +
    geom_area(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = wgEncodeBroadHistoneFileDesc$ColourDark) }

#output
gg_out +
  scale_y_continuous(breaks=50,
                     labels=udf_pad("Histone"),
                     name="") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() 

