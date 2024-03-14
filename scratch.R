

install.packages('BiocManager')
BiocManager::install(c('tidyverse','stringr','biovizBase','GenomicRanges','ggbio','knitr'))
install.packages(c("shiny","dplyr","tidyr","lazyeval","data.table","ggplot2","ggrepel","knitr","markdown","DT","lattice","acepack","cluster","DBI","colourpicker","igraph","visNetwork", "devtools"))
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","rtracklayer"))


#devtools::install_github("oncogenetics/oncofunco")

detach("package:oncofunco", unload = TRUE)
devtools::install_github("rachelicr/oncofunco@dev",force=TRUE)