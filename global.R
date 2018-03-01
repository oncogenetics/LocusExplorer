# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "shinyjs", "dtplyr", "data.table", "dplyr", "lazyeval","tidyr",
#                    "ggplot2", "ggrepel", "knitr", "markdown", "stringr","DT","seqminer",
#                    "lattice", "cluster", "DBI", "networkD3", "scales",
#                    "googleVis"))

library(shiny)
#library(shinyjs) #colourInput
#library(dtplyr)
library(dplyr) 
library(tidyr)
library(lazyeval) # dplyr
library(data.table)
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
library(colourpicker)

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


#Map fonts to Windows
if(Sys.info()['sysname'] == "Windows") {
  windowsFonts(Courier=windowsFont("TT Courier New")) 
}


#Custom functions - replaced by "oncofunco" package
#source("Source/UDF.R", local = TRUE)

# Data --------------------------------------------------------------------
# 01. prostate data region names ------------------------------------------
regions <- fread("Data/ProstateData/regions.csv",
                 data.table = FALSE, header = TRUE) %>% 
  mutate(REGIONBED = paste0(CHR,"_",START,"_",END))

# 02. OncoFinemap: annot ---------------------------------------------------
annotOncoFinemap <- fread("Data/Annotation/Supp_Table2_clean_annot.csv")
annotOncoFinemapEQTL <- fread("Data/Annotation/Supp_Table2_clean_eqtl.csv")

# 03. wgEncodeBroadHistone --------------------------------------------------
#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
           stringsAsFactors = FALSE)

# UDF ---------------------------------------------------------------------

