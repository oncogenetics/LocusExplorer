# About -------------------------------------------------------------------
# 09/03/2018
# - Supp_Data1, prepare data for LE annotation tracks

# Workspace ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)

setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/LocusExplorer")
#setwd("N:/Translational Cancer Genetics Team/Bioinformatics/Development/R_ShinyApps/LocusExplorer")

# Data --------------------------------------------------------------------
annot <- fread("Data/Annotation/Supp_Data1.csv", na.strings = c("", "NA"),
               data.table = FALSE, header = TRUE)

# Annotation --------------------------------------------------------------
annotCols <- setNames(
  c("#D3D3D3", "#FFA500", "#7070FF", "#00DD00", "#0000AA", "#009900",
    "#E0E0E0", "#ED3333", "#4FFFA4", "#C9DD03", "#F9A100"),
  c("Heterochromatin",	"CTCF", "CTCF+Enhancer", "Promoter",
    "Enhancer", "Poised_Promoter", "Transcribed", "Repressed",
    "CTCF+Promoter", "DNaseI", "Conserved"))

annotOncoFinemap <- annot %>% 
  transmute(CHR = Variant_Chromosome,
            BP = Variant_Position,
            DNaseI = (GSM1008595_RWPE1_DNaseI + GSM1024742_PrEC_DNaseI +
                        GSM1024743_PrEC_DNaseI + GSM736565_LNCaP_DNaseI +
                        GSM736603_LNCaP_DNaseI + GSM816634_LNCaP_androgen_DNaseI +
                        GSM816637_LNCaP_DNaseI) > 0,
            Conserved = (GERP__ + SiPhy_Omega + SiPhy_Pi + PhastCons) > 0,
            #PrEC_ChromHMM
            ChromHMM = PrEC_ChromHMM) 
# wide to long
annotOncoFinemap <- gather(annotOncoFinemap, key = TYPE1, value = TYPE2, -c(CHR, BP)) %>% 
  filter(TYPE2 != "FALSE") %>%
  mutate(TYPE2 = if_else(TYPE2 == "TRUE", TYPE1, TYPE2)) %>% 
  arrange(CHR, BP)
# add colors names
annotOncoFinemap$COLOUR_HEX <- annotCols[ annotOncoFinemap$TYPE2 ]
# TypeN, used for plotting in order on yaxis
annotOncoFinemap$TYPE2N <-
  as.numeric(
    factor(annotOncoFinemap$TYPE2,
           levels = c(
             "DNaseI","Conserved",
             #ChromHMM
             "Heterochromatin","CTCF","CTCF+Enhancer","Promoter","Enhancer",
             "Poised_Promoter","Transcribed","Repressed","CTCF+Promoter")))

#head(annotOncoFinemap, 2)
#   CHR        BP    TYPE1           TYPE2 COLOUR_HEX TYPE2N
# 1   1 150267316 ChromHMM        Promoter    #00DD00      6
# 2   1 150270088 ChromHMM Heterochromatin    #D3D3D3      3

# eQTL --------------------------------------------------------------------
annotOncoFinemapEQTL <-
  annot %>%
  filter( !is.na(TCGA_eQTL_coloc0.9_genes)) %>%
  dplyr::transmute(CHR = Variant_Chromosome,
                   SNPhit = PP_best_tag_dbSNP,
                   SNPhit_BP = PP_best_tag_position,
                   SNP = dbSNP,
                   SNP_BP = Variant_Position,
                   R2 = PP_tag_r.2,
                   GENE = TCGA_eQTL_coloc0.9_genes
                   ) %>%
  base::unique()

# split stack genes
annotOncoFinemapEQTL <- cSplit(annotOncoFinemapEQTL, splitCols = "GENE",
                               sep = ";", direction = "long")
annotOncoFinemapEQTL <- annotOncoFinemapEQTL %>% 
  mutate(R2 = ifelse(SNPhit == SNP, 1, R2)) %>% 
  arrange(CHR, SNP_BP)

# head(annotOncoFinemapEQTL, 2)
#   CHR      SNPhit SNPhit_BP         SNP    SNP_BP R2     GENE
# 1   1  rs61305632 150267316  rs61305632 150267316  1 ADAMTSL4
# 2   1 rs570851450 150334248 rs570851450 150334248  1 ADAMTSL4

# Output ------------------------------------------------------------------
write.csv(annotOncoFinemap, "Data/Annotation/Supp_Data1_clean_annot.csv", row.names = FALSE)
write.csv(annotOncoFinemapEQTL, "Data/Annotation/Supp_Data1_clean_eqtl.csv", row.names = FALSE)

