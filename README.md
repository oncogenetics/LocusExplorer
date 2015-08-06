Locus Explorer v0.2
=============
## An interactive graphical illustration of genetic associations and their biological context.

### Disclaimer
Locus Explorer should be used for illustrative purposes only. Any results provided by Locus Explorer should be used with caution. 

### Availability  
The source code and installation instructions for Locus Explorer are available at https://github.com/oncogenetics/LocusExplorer.

Locus Explorer is made available under the xxxx license.

### Required Software
Locus Explorer runs in the R environment but is designed to be an easy to use interface that does not require familiarity with R as a prerequisite. Locus explorer is platform agnostic and able to run on any operating system for which R is available.

The latest version of R can be downloaded by following the instructions at [https://www.r-project.org/](https://www.r-project.org/])

After installation of  the R software, R packages used by Locus Explorer must be installed prior to use. This may take a few minutes, but is only required on the first occasion. To install packages, open the R program, copy the following code into the R console and hit Return:
```R
install.packages(c("shiny", "data.table", "dplyr", "tidyr", "ggplot2", "knitr", "markdown", "stringr","DT"), dependencies = TRUE)

source("http://bioconductor.org/biocLite.R")
biocLite(c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db"))
```

### Launching the Locus Explorer Application
To start Locus Explorer, open R, copy the following code into the console and hit Return:
```R
library(shiny)  
runGitHub("LocusExplorer", "oncogenetics")
```

Locus Explorer runs through a web browser and uses an intuitive interface that does not require high level computational skills to operate

#### Troubleshooting
If you see the following error:
> Error in download.file(URL, destfile = ...) : 
>   unsupported URL scheme

Try running:
```R
setInternet2(TRUE)
```

### To Cite Locus Explorer
[Locus Explorer: a user-friendly tool for integrated visualisation of genetic association data and biological annotations](url to add) Tokhir Dadaev<sup>1</sup>, Daniel A Leongamornlert<sup>1</sup>, Edward J Saunders<sup>1</sup>, Zsofia Kote-Jarai<sup>1</sup>
<sup>1</sup>The Institute of Cancer Research, London
#### Abstract
**Summary:** In this article we present Locus Explorer, a data visualisation and exploration tool for genetic association data. Locus Explorer is written in R using the Shiny library, providing access to powerful R-based functions and libraries through a simple user interface. Locus Explorer allows users to simultaneously display genetic, statistical and biological data in a single image and allows dynamic zooming and customisation of the plot display. Publication quality plots may be downloaded in a variety of file formats.  
**Availability:** Locus Explorer is open source, runs through R and a web browser, and available at https://github.com/oncogenetics/LocusExplorer, where additional README information and example data are also provided.

### Additional Publications That Use Locus Explorer
[Multiple novel prostate cancer susceptibility signals identified by fine-mapping of known risk loci among Europeans.](http://www.ncbi.nlm.nih.gov/pubmed/26025378) Al Olama AA *et al.*

### Contact  
Questions, suggestions, and bug reports are welcome and appreciated. 
- Submit [suggestions and bug-reports](https://github.com/oncogenetics/LocusExplorer/issues)
- Send [pull request](https://github.com/oncogenetics/LocusExplorer/pulls)
- Contact email [Dadaev T](mailto: tokhir.dadaev@icr.ac.uk)

### To-do list
-	LNCAP arrows for prostate data.
-	Plot legend
-	Make LD file optional
-	Bed track custom, multiple files
-	Input file validation
-	Encode - wiggle data as bed track
-	LD tutorial add tabix/plink option
-	Multiple region upload
-	Hide red error messages



