# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Genetic Map, recombination data
# - download 1KG gemetic map files chr1:23
#   http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1interim_jun2011_impute/
# - read into R, clean up
# - use bgzip and tabix
#      - bgzip GeneticMap1KG.txt
#      - tabix -s 1 -b 2 -e 2 GeneticMap1KG.txt.gz
#      - issue "tabix.read.table" function throwing error... 
#              - manual conversion to data.frame in server.R


# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)

setwd("PATH TO GENETIC MAP FILES")

for(i in paste0("chr",1:23)) {
  #i="chr20"
  genMapFile <- paste0("Data/GeneticMap1KG/genetic_map_",i,"_combined_b37.txt")
  genMap <- fread(genMapFile,sep = " ", header = TRUE)
  write.table(genMap %>% 
                transmute(CHR=i,
                          BP=position,
                          RECOMB=round(`COMBINED_rate(cM/Mb)`,2)) %>% 
                filter(RECOMB>0) %>% 
                arrange(BP),
              "Data/GeneticMap1KG/GeneticMap1KG.txt",
              quote=FALSE,row.names=FALSE,col.names = FALSE,sep="\t",
              append= TRUE)
  rm(genMap)
  gc()
}





