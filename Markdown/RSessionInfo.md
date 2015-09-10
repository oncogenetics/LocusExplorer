### Required libraries
```R
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "tidyr", "ggplot2",
#                    "knitr", "markdown", "stringr","DT","seqminer",
#                    "lattice","cluster"),
#                  dependencies = TRUE)
library(shiny)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(markdown)
library(stringr)
library(DT)
library(seqminer)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "org.Hs.eg.db","rtracklayer"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols
library(rtracklayer) # bigwig
```
### Session info
Application tested on 10/09/2015 16:24

```R
> sessionInfo()
R version 3.2.2 (2015-08-14)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rtracklayer_1.28.10                     org.Hs.eg.db_3.1.2                      RSQLite_1.0.0                          
 [4] DBI_0.3.1                               TxDb.Hsapiens.UCSC.hg19.knownGene_3.1.2 GenomicFeatures_1.21.16                
 [7] AnnotationDbi_1.31.17                   Biobase_2.28.0                          GenomicRanges_1.20.6                   
[10] GenomeInfoDb_1.5.12                     IRanges_2.3.18                          S4Vectors_0.7.13                       
[13] ggbio_1.16.1                            BiocGenerics_0.15.6                     seqminer_4.7                           
[16] DT_0.1                                  stringr_1.0.0                           markdown_0.7.7                         
[19] knitr_1.11                              ggplot2_1.0.1                           tidyr_0.3.1                            
[22] dplyr_0.4.3                             data.table_1.9.4                        shiny_0.12.2                           

loaded via a namespace (and not attached):
 [1] splines_3.2.2              Formula_1.2-1              assertthat_0.1             latticeExtra_0.6-26       
 [5] RBGL_1.45.1                BSgenome_1.37.4            Rsamtools_1.21.15          lattice_0.20-33           
 [9] biovizBase_1.16.0          chron_2.3-47               digest_0.6.8               RColorBrewer_1.1-2        
[13] XVector_0.9.1              colorspace_1.2-6           htmltools_0.2.6            httpuv_1.3.3              
[17] plyr_1.8.3                 OrganismDbi_1.11.42        XML_3.98-1.3               biomaRt_2.25.1            
[21] zlibbioc_1.15.0            xtable_1.7-4               scales_0.3.0               BiocParallel_1.3.48       
[25] SummarizedExperiment_0.3.3 nnet_7.3-11                proto_0.3-10               survival_2.38-3           
[29] magrittr_1.5               mime_0.4                   GGally_0.5.0               MASS_7.3-43               
[33] foreign_0.8-66             graph_1.47.2               BiocInstaller_1.18.4       tools_3.2.2               
[37] munsell_0.4.2              cluster_2.0.3              lambda.r_1.1.7             Biostrings_2.36.4         
[41] futile.logger_1.4.1        grid_3.2.2                 RCurl_1.95-4.7             dichromat_2.0-0           
[45] VariantAnnotation_1.15.26  htmlwidgets_0.5            bitops_1.0-6               gtable_0.1.2              
[49] reshape_0.8.5              reshape2_1.4.1             R6_2.1.1                   GenomicAlignments_1.4.1   
[53] gridExtra_2.0.0            Hmisc_3.16-0               futile.options_1.0.0       stringi_0.5-5             
[57] Rcpp_0.12.0                rpart_4.1-10               acepack_1.3-3.3  
```
