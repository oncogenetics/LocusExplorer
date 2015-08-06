### Required libraries
```R
#CRAN
# install.packages(c("shiny", "data.table", "dplyr", "tidyr", ggplot2",
#                    "knitr", "markdown", "stringr","DT"),
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

#Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "org.Hs.eg.db"))
library(ggbio)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # gene symobols
```
### Session info
Application tested on 06/08/2015 14:42

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
 [1] org.Hs.eg.db_3.1.2                      RSQLite_1.0.0                           DBI_0.3.1                              
 [4] TxDb.Hsapiens.UCSC.hg19.knownGene_3.1.2 GenomicFeatures_1.20.1                  AnnotationDbi_1.30.1                   
 [7] Biobase_2.28.0                          GenomicRanges_1.20.5                    GenomeInfoDb_1.4.1                     
[10] IRanges_2.2.5                           S4Vectors_0.6.1                         ggbio_1.16.0                           
[13] BiocGenerics_0.14.0                     DT_0.1                                  stringr_1.0.0                          
[16] markdown_0.7.7                          knitr_1.10.5                            ggplot2_1.0.1                          
[19] tidyr_0.2.0                             dplyr_0.4.2                             data.table_1.9.4                       
[22] shiny_0.12.1                           

loaded via a namespace (and not attached):
 [1] splines_3.2.1            Formula_1.2-1            assertthat_0.1           latticeExtra_0.6-26     
 [5] RBGL_1.44.0              BSgenome_1.36.2          Rsamtools_1.20.4         lattice_0.20-31         
 [9] biovizBase_1.16.0        chron_2.3-47             digest_0.6.8             RColorBrewer_1.1-2      
[13] XVector_0.8.0            colorspace_1.2-6         htmltools_0.2.6          httpuv_1.3.2            
[17] plyr_1.8.3               OrganismDbi_1.10.0       XML_3.98-1.3             biomaRt_2.24.0          
[21] zlibbioc_1.14.0          xtable_1.7-4             scales_0.2.5             BiocParallel_1.2.7      
[25] nnet_7.3-9               proto_0.3-10             survival_2.38-3          magrittr_1.5            
[29] mime_0.3                 GGally_0.5.0             MASS_7.3-42              foreign_0.8-63          
[33] graph_1.46.0             tools_3.2.1              munsell_0.4.2            cluster_2.0.1           
[37] lambda.r_1.1.7           Biostrings_2.36.1        futile.logger_1.4.1      grid_3.2.1              
[41] RCurl_1.95-4.7           dichromat_2.0-0          VariantAnnotation_1.14.6 htmlwidgets_0.5         
[45] bitops_1.0-6             gtable_0.1.2             reshape_0.8.5            reshape2_1.4.1          
[49] R6_2.1.0                 GenomicAlignments_1.4.1  gridExtra_0.9.1          rtracklayer_1.28.6      
[53] Hmisc_3.16-0             futile.options_1.0.0     stringi_0.5-5            Rcpp_0.11.6             
[57] rpart_4.1-9              acepack_1.3-3.3     
```

