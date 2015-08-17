# About -------------------------------------------------------------------
# Date: 17/08/2015
# prepare LNCAP data

# Workspace ---------------------------------------------------------------
setwd("N:/Translational Cancer Genetics Team/04 Zsofia/20 Fine Map/LocusExplorer/LocusExplorer/")

library(data.table)
library(dplyr)

# Data prep ---------------------------------------------------------------
# LNCAP file
LNCAP <- fread("Data/ProstateLNCAP/tokhir.matrix.csv")
#head(LNCAP)
res <- LNCAP %>% 
  transmute(
    CHR=chromosome,
    BP=end,
    ANNOTATE=ifelse(Func.refGene %in% c("exonic","intergenic",
                                        "intronic","ncRNA_exonic",
                                        "ncRNA_intronic","splicing",
                                        "UTR3","UTR5"),
                    Func.refGene,NA)) %>% 
  #sorty by chr and bpo
  mutate(CHRN=gsub("chr","",CHR),
         CHRN=as.integer(ifelse(CHRN=="X",23,CHRN))) %>% 
  arrange(CHRN,BP) %>% 
  select(CHR,BP,ANNOTATE)

# Output ------------------------------------------------------------------
#output file for tabix
write.table(res,"Data/ProstateLNCAP/ProstateLNCAP.txt",
            quote=FALSE,row.names=FALSE,sep="\t")
