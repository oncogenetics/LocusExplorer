# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "data.table", "dplyr", "lazyeval",
#                    "ggplot2", "ggrepel", "knitr", "markdown","DT",
#                    "lattice", "acepack", "cluster", "DBI", "scales",
#                    "colourpicker"))

library(shiny)
library(dplyr) 
library(tidyr)
library(lazyeval) # dplyr
library(data.table)
library(ggplot2)
library(ggrepel) # geom_text_repel
library(knitr)
library(markdown)
library(DT)
library(lattice) # ggbio
library(acepack) # ggbio
library(cluster)
library(DBI) # gene plot
library(colourpicker)

library(igraph)
library(visNetwork)

#Bioconductor packages
# source("http://bioconductor.org/biocLite.R")
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




# Data --------------------------------------------------------------------
# 01. prostate data region names ------------------------------------------
regions <- fread("Data/ProstateData/regions.csv",
                 data.table = FALSE, header = TRUE) %>% 
  mutate(REGIONBED = paste0(CHR,"_",START,"_",END))

# 02. OncoFinemap: annot ---------------------------------------------------
annotOncoFinemap <- fread("Data/Annotation/Supp_Data1_clean_annot.csv")
annotOncoFinemapEQTL <- fread("Data/Annotation/Supp_Data1_clean_eqtl.csv")

# 03. wgEncodeBroadHistone ------------------------------------------------
#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
           stringsAsFactors = FALSE)
# 04. GeneticMap1KG -------------------------------------------------------
GeneticMap1KG <- readRDS("Data/GeneticMap1KG/GeneticMap1KG.RData")


# TESTING -----------------------------------------------------------------
#removed dependency on
#library(tidyr)
#library(shinyjs) #colourInput
#library(dtplyr)
#library(stringr)
#library(seqminer) # query tabix
#Custom functions - replaced by "oncofunco" package
#source("Source/UDF.R", local = TRUE)

#Custom functions - replaced by "oncofunco" package
#source("Source/UDF.R", local = TRUE)
