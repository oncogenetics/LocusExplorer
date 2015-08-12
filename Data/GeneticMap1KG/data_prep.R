# About -------------------------------------------------------------------
# Genetic Map, recombination data
# - download 1KG gemetic map files chr1:23
#   http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1interim_jun2011_impute/
# - read into R, clean up
# - use bgzip and tabix
#      - bgzip GeneticMap1KG.txt
#      - tabix -s 1 -b 2 -e 2 GeneticMap1KG.txt.gz
#      - issue "tabix.read.table" function throwing error... 
#              - manual conversion to data.frame


# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)

setwd("N:/Translational Cancer Genetics Team/04 Zsofia/20 Fine Map/LocusExplorer/LocusExplorer")

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


# TESTING -----------------------------------------------------------------
# library(seqminer)
# library(tidyr)

# x <- tabix.read("Data/GeneticMap1KG/GeneticMap1KG.txt.gz","chr1:721290-729948")
# x <- data.frame(temp=gsub("\r","",x,fixed=TRUE))
# x <- separate(x,col=temp,into=c("CHR","BP","RECOMB"),sep="\t",convert=TRUE)
# str(x)

 
# works, returns char vector...
# tabix.read("Data/GeneticMap1KG/GeneticMap1KG.txt.gz","chr1:721290-729948")

# FAILs
# tabix.read.table("Data/GeneticMap1KG/GeneticMap1KG.txt.gz","chr1:721290-729948", col.names = FALSE)

# example from manuals
# fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
# snp <- tabix.read.table(fileName, "1:196623337-196632470")




