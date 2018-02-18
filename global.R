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

# 02 OncoFinemap: annot ---------------------------------------------------
annotOncoFinemap <- fread("Data/Annotation/Supp_Table2_clean_annot.csv")
annotOncoFinemapEQTL <- fread("Data/Annotation/Supp_Table2_clean_eqtl.csv")

# 04. wgEncodeBroadHistone --------------------------------------------------
#wgEncodeBroadHistone bigwig data description
wgEncodeBroadHistoneFileDesc <-
  read.csv("Data/wgEncodeBroadHistone/wgEncodeBroadHistone.csv",
           stringsAsFactors = FALSE)


# UDF ---------------------------------------------------------------------
plotAnnot <- function(data,
                      chrom = NULL,
                      chromStart = NULL,
                      chromEnd = NULL,
                      vline = NULL,
                      collapse = FALSE){
  #subset data for zoomed region
  data <- data[ CHR == chrom &
                  BP >= chromStart &
                  BP <= chromEnd, ]
  
  if(collapse){ data[ TYPE1 == "ChromHMM", TYPE2N := 3] } 
  
  #plot
  gg_out <- 
    ggplot(data,
           aes(xmin = BP, xmax = BP, ymin = TYPE2N - 1, ymax = TYPE2N,
               col = COLOUR_HEX, fill = COLOUR_HEX )) +
    geom_rect(alpha = ifelse(collapse, 0.5, 1)) 
  
  if(collapse){
    gg_out <- gg_out +
      scale_y_continuous(
        limits = c(-1, 4),
        breaks = c(0:3) + 0.5, 
        labels = c("DNaseI", "Conserved", "ChromHMM", "eQTL"),
        name = "")
  } else {
    gg_out <- gg_out +
      scale_y_continuous(
        limits = c(-1, 11),
        breaks = c(0:10) + 0.5, 
        labels = c("DNaseI","Conserved",
                   #ChromHMM
                   "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
                   "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter"),
        name = "")
    
  }
  
  #prettify
  gg_out <- gg_out +
    scale_color_identity() +
    # general options
    coord_cartesian(xlim = c(chromStart, chromEnd))
    # xlim(c(chromStart, chromEnd)) +
    # theme_LE() +
    # theme(
    #   axis.text.x = element_blank(),
    #   axis.ticks.x = element_blank())
  
  # mark hit SNPs - vline
  if(!is.null(vline)){
    gg_out <- gg_out +
      geom_vline(xintercept = vline,
                 col = "black", #col = "#4daf4a",
                 linetype = "dashed")
  }
  
  #return output ggplot
  return(gg_out)
} # END plotAnnot



NCBIrsidHyperlink <- function(SNP,
                              link = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
                              html = TRUE){
  
  if(html){
    ifelse(substr(SNP, 1, 2) == "rs",
           paste0('<a href="', link, gsub("rs", "", SNP, fixed = TRUE),
                  '" target="_blank">', SNP, '</a>'),
           SNP)
  } else {
    ifelse(substr(SNP, 1, 2) == "rs",
           paste0(link, gsub("rs", "", SNP, fixed = TRUE)),
           SNP)
  }
} #END NCBIrsidHyperlink


theme_LE <- function(){
  # LocusExplorer ggplot custom theme
  # General options for all plots -------------------------------------------
  # Usage: ggplot() + theme_LE()
  
  theme(legend.position = "none",
        #panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill="darkseagreen"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",
                                          linetype = "dotted"),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(colour = "grey20",
                                   margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "cm"),
                                   face = "italic"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 2, unit = "cm")),
        plot.margin = unit(c(0, 0, 0, 1), "cm")
        #plot.margin = unit(c(1, 1, 0.5, 0.5), "lines")
        #plot.margin = rep(unit(0,"null"),4)
        )
  
  }
