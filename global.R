# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "shinyjs", "data.table", "dplyr", "lazyeval","tidyr",
#                    "ggplot2", "knitr", "markdown", "stringr","DT","seqminer",
#                    "lattice", "cluster", "DBI", "networkD3", "scales",
#                    "googleVis"))

library(shiny)
library(shinyjs) #colourInput
library(dtplyr)
library(data.table)
library(dplyr) 
library(lazyeval) # dplyr
library(tidyr)
library(ggplot2)
library(ggrepel) # geom_text_repel
library(knitr)
library(markdown)
library(stringr)
library(DT)
library(seqminer) # query tabix
library(lattice) # ggbio
library(acepack) # ggbio
library(cluster)
library(DBI) # gene plot

#Bioconductor packages
#source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "org.Hs.eg.db","rtracklayer"))

library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols
library(rtracklayer) # bigwig

#GitHub packages
# install.packages("devtools")
# devtools::install_github("oncogenetics/oncofunco")
library(oncofunco)

# increase upload limit to 20Mb
options(shiny.maxRequestSize = 20 * 1024 ^ 2)

#Custom functions - replaced by "oncofunco" package
#source("Source/UDF.R", local = TRUE)

# Data --------------------------------------------------------------------
# prostate data region names
regions <- fread("Data/ProstateData/regions.csv",
                 data.table = FALSE, header = TRUE) %>% 
  mutate(REGIONBED = paste0(CHR,"_",START,"_",END))


annot <- fread("Data/Annotation/Supp_Table2.csv", na.strings = c("", "NA"),
               data.table = FALSE)
annotCols <- setNames(
  c("#D3D3D3", "#FFA500", "#7070FF", "#00DD00", "#0000AA", "#009900",
    "#E0E0E0", "#ED3333", "#4FFFA4", "#C9DD03", "#F9A100"),
  c("Heterochromatin",	"CTCF", "CTCF+Enhancer", "Promoter",
    "Enhancer", "Poised_Promoter", "Transcribed", "Repressed",
    "CTCF+Promoter", "DNaseI", "Conserved"))

#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
           stringsAsFactors = FALSE)

