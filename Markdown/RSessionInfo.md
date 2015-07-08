### Required libraries
```R
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "ggplot2", "ggbio", "knitr", "markdown"),dependencies = TRUE)
library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(markdown)

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols
```
### Session info
```R
R version 3.2.1 (2015-06-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stringr_1.0.0                           org.Hs.eg.db_3.1.2                      RSQLite_1.0.0                          
 [4] DBI_0.3.1                               TxDb.Hsapiens.UCSC.hg19.knownGene_3.1.2 GenomicFeatures_1.20.1                 
 [7] AnnotationDbi_1.30.1                    Biobase_2.28.0                          GenomicRanges_1.20.5                   
[10] GenomeInfoDb_1.4.1                      IRanges_2.2.5                           S4Vectors_0.6.1                        
[13] ggbio_1.16.0                            BiocGenerics_0.14.0                     markdown_0.7.7                         
[16] knitr_1.10.5                            ggplot2_1.0.1                           dplyr_0.4.2                            
[19] data.table_1.9.4                        shiny_0.12.1                           

loaded via a namespace (and not attached):
 [1] jsonlite_0.9.16          splines_3.2.1            Formula_1.2-1            assertthat_0.1           latticeExtra_0.6-26      RBGL_1.44.0             
 [7] BSgenome_1.36.1          Rsamtools_1.20.4         lattice_0.20-31          biovizBase_1.16.0        chron_2.3-47             digest_0.6.8            
[13] RColorBrewer_1.1-2       XVector_0.8.0            colorspace_1.2-6         htmltools_0.2.6          httpuv_1.3.2             plyr_1.8.3              
[19] OrganismDbi_1.10.0       XML_3.98-1.3             biomaRt_2.24.0           zlibbioc_1.14.0          xtable_1.7-4             scales_0.2.5            
[25] BiocParallel_1.2.6       lazyeval_0.1.10          nnet_7.3-9               proto_0.3-10             survival_2.38-1          magrittr_1.5            
[31] mime_0.3                 GGally_0.5.0             MASS_7.3-40              foreign_0.8-63           Cairo_1.5-6              graph_1.46.0            
[37] tools_3.2.1              munsell_0.4.2            cluster_2.0.1            lambda.r_1.1.7           Biostrings_2.36.1        futile.logger_1.4.1     
[43] grid_3.2.1               RCurl_1.95-4.7           rstudioapi_0.3.1         dichromat_2.0-0          VariantAnnotation_1.14.4 rmarkdown_0.7           
[49] labeling_0.3             bitops_1.0-6             gtable_0.1.2             reshape_0.8.5            reshape2_1.4.1           R6_2.0.1                
[55] GenomicAlignments_1.4.1  gridExtra_0.9.1          rtracklayer_1.28.5       Hmisc_3.16-0             futile.options_1.0.0     stringi_0.5-5           
[61] Rcpp_0.11.6              rpart_4.1-9              acepack_1.3-3.3    
```

