# ENCODE wgEncodeRegDnaseClustered

ggplot(ROIdatwgEncodeRegDnaseClustered(),
       aes(xmin=START,xmax=END,ymin=0,ymax=1,fill=SCORE)) +
  geom_rect(alpha=0.5) +
  scale_fill_continuous(limits=c(500,1000),low="white",high="#A71930",na.value = "white") +
  scale_y_continuous(breaks=0.5,
                     labels=udf_pad("Dnase"),
                     name="xxx") +
  #general options
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme()



# TESTING -----------------------------------------------------------------
# ggplot(x,aes(BP,SCORE,fill=ENCODE)) +
#   geom_area(alpha = 0.5, position = "identity") +
#   scale_fill_manual(values = EncodeFileDesc$ColourDark) +
#   theme_classic()
# 
# ggplot(x,aes(BP,SCORE,fill=ENCODE)) +
#   geom_area(position="stack") +
#   theme_classic()
# 
# ggplot(x,aes(BP,SCORE,fill=ENCODE)) +
#   geom_area(col="grey") +
#   theme_classic()
# 
# ggplot(x,aes(BP,SCORE,fill=ENCODE)) +
#   geom_area(col="grey") +
#   facet_grid(ENCODE~.) +
#   theme_classic()
# 
# # roi <- GRanges(seqnames = "chrX",
# #                IRanges(start = 66820000,
# #                        end = 66840000))

