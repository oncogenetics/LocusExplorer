# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "shinyjs", "data.table", "dplyr", "lazyeval","tidyr",
#                    "ggplot2", "knitr", "markdown", "stringr","DT","seqminer",
#                    "lattice", "cluster", "DBI", "networkD3", "scales",
#                    "googleVis"))

library(shiny)
#library(shinyjs)
library(dtplyr)
library(data.table)
library(dplyr) 
library(lazyeval) # dplyr
library(tidyr)
library(ggplot2)
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

# increase upload limit to 20Mb
options(shiny.maxRequestSize = 20 * 1024 ^ 2)

#Custom functions
source("Source/UDF.R", local = TRUE)

# Data --------------------------------------------------------------------
# prostate data region names
regions <- fread("Data/ProstateData/regions.csv",
                 data.table = FALSE, header = TRUE) %>% 
  mutate(REGIONBED = paste0(CHR,"_",START,"_",END))

# oncoarray finemapping annotation + TCGA
annotOncoFinemap <- fread("Data/Annotation/OA_meta_All_SNPs_Annotated_150716.tsv")
# annotOncoFinemap$Gene <-
#   sapply(annotOncoFinemap$eqtl.Eze.gene,
#          function(i)unlist(strsplit(i, split = "|", fixed = TRUE))[1])

annotOncoFinemapEQTL <- fread("Data/Annotation/LE_finemap_eQTL.tsv")
#annotOncoFinemapEQTL[!is.na(annotOncoFinemapEQTL$geneStart),]

annotOncoNewHits <- fread("Data/Annotation/LE_OA_65_new_hits_Annotation_2016-07-22.tsv")

annotOncoNewHitsEQTL <- fread("Data/Annotation/OA_NewHits_eQTLs.tsv")


#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
           stringsAsFactors = FALSE)

