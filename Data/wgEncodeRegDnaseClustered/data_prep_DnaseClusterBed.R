
library(data.table)
library(ggplot2)
library(dplyr)


setwd("D:/Tokhir/TEMP/LocusExplorer")


dat <- fread("Data/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed")


dat <- dat %>% 
  transmute(CHR=V1,
            START=V2,
            END=V3,
            SCORE=V5)


plot(density(dat$SCORE))
summary(dat$SCORE)
chr5:1,170,000-1,410,000

plotDat <- dat %>% filter(CHR=="chr5",
                        START>=1170000,
                        END<=1410000)

ggplot(plotDat,aes(START,SCORE)) +
  geom_area()


summary(plotDat$SCORE)

ggplot(plotDat,aes(xmin=START,xmax=END,ymin=0,ymax=1,fill=SCORE)) +
  geom_rect(alpha=0.5) +
  scale_fill_continuous(limits=c(500,1000),low="grey50",high="red",na.value = "grey50") +
  theme_classic()
  
ggplot(plotDat,aes(START+(END-START)/2,y = SCORE)) +
  geom_line()


geom_line(aes(x = START+(END-START)/2,y = SCORE/1000,alpha=0.5),plotDat) +