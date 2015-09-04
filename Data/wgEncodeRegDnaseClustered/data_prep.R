# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# 14/08/2015
# Process DnaseCluster bed file for tabixing
# input bed: 
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz


# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)
#library(ggplot2)


setwd("D:/Tokhir/TEMP/LocusExplorer")


# Data prep ---------------------------------------------------------------
dat <- fread("Data/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed")

dat <- dat %>% 
  transmute(CHR=V1,
            START=V2,
            END=V3,
            SCORE=V5)

write.table(dat,"Data/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# head(dat)
# summary(dat$SCORE)

# TESTING -----------------------------------------------------------------

# plotting example
# #TERT chr5:1,170,000-1,410,000
# plotDat <- dat %>% filter(CHR=="chr5",
#                           START>=1170000,
#                           END<=1410000)

# ggplot(plotDat,aes(xmin=START,xmax=END,ymin=0,ymax=1,fill=SCORE)) +
#   geom_rect(alpha=0.5) +
#   scale_fill_continuous(limits=c(500,1000),low="grey50",high="red",na.value = "grey50") +
#   theme_classic()


# ggplot(plotDat,aes(START,SCORE)) +
#   geom_area()
# 
#   
# ggplot(plotDat,aes(START+(END-START)/2,y = SCORE)) +
#   geom_line()
# 
# 
# geom_line(aes(x = START+(END-START)/2,y = SCORE/1000,alpha=0.5),plotDat) +