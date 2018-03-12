# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Genetic Map, recombination data
# - download 1KG gemetic map files chr1:23
#   http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1interim_jun2011_impute/
# - read into R, clean up, output as one object RData file


# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)

# setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/LocusExplorer")

res <- lapply(1:23, function(i){

  if(i == 23) {CHR <- "X_nonPAR"} else {CHR <-  i}
  genLink <- paste0("http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr", CHR, "_combined_b37.txt")
  genMap <- fread(genLink)
  
  genMap %>% 
    transmute(CHR = i,
              BP = position,
              RECOMB = round(`COMBINED_rate(cM/Mb)`, 2)) %>% 
    filter(RECOMB > 0) %>% 
    arrange(BP)
})

GeneticMap1KG <- rbindlist(res)
setkeyv(GeneticMap1KG, c("CHR", "BP"))

saveRDS(GeneticMap1KG, "Data/GeneticMap1KG/GeneticMap1KG.RData")
