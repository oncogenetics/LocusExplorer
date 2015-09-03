Locus Explorer v0.2
=============
## An interactive graphical illustration of genetic associations and their biological context.

### Disclaimer
Locus Explorer should be used for illustrative purposes only. Any results provided by Locus Explorer should be used with caution. 

### Availability  
The source code and installation instructions for Locus Explorer are available at https://github.com/oncogenetics/LocusExplorer.

Locus Explorer is made available under the MIT license.

### Required Software
Locus Explorer runs in the R environment but is designed to be an easy to use interface that does not require familiarity with R as a prerequisite. Locus explorer is platform agnostic and able to run on any operating system for which R is available.

Locus Explorer requires R version 3.2.2 to run and can be downloaded by following the instructions at [https://www.r-project.org/](https://www.r-project.org/). Some required packages are not available for earlier versions of R.

After installation of  the R software, R packages used by Locus Explorer must be installed prior to use. This may take a few minutes, but is only required on the first occasion. To install packages, open the R program, copy the following code into the R console and hit Return:
```R
#install CRAN packages, if missing
packages <- c("shiny", "data.table", "dplyr", "tidyr", "ggplot2", "knitr", "markdown", "stringr","DT","seqminer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)  
} else { print("All required CRAN packages installed")}

#install Bioconductor packages if missing
source("https://bioconductor.org/biocLite.R")
bioc <- c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db")
if (length(setdiff(bioc, rownames(installed.packages()))) > 0) {
  biocLite(setdiff(bioc, rownames(installed.packages())))  
} else { print("All required Bioconductor packages installed")}
```
- In cases when user do not have admin rights, pop up window will prompt to set a personal library location for installation of packages, please click yes.
- If using R GUI then user might get prompted to choose CRAN mirror to use for package downloads, please choose the city nearer to your location.
- If prompted to "Update packages all/some/none [a/s/n]", type "n" and hit Return.


### Launching Locus Explorer Application
To start Locus Explorer, open R (RStudio recommended), copy the following code into the console and hit Return:
```R
library(shiny)  
runGitHub("LocusExplorer", "oncogenetics")
```

Locus Explorer runs through a web browser and uses an intuitive interface that does not require high level computational skills to operate.

#### Troubleshooting
If you see the following error:
> Error in download.file(URL, destfile = ...) : 
>   unsupported URL scheme

Try running:
```R
setInternet2(TRUE)
```

### To Cite Locus Explorer
[Locus Explorer: a user-friendly tool for integrated visualisation of genetic association data and biological annotations](url to add) Tokhir Dadaev<sup>1</sup>, Daniel A Leongamornlert<sup>1</sup>, Edward J Saunders<sup>1</sup>, Rosalind Eeles<sup>1,2</sup> , Zsofia Kote-Jarai<sup>1</sup>   

<sup>1</sup>Department of Genetics and Epidemiology, The Institute of Cancer Research, London, UK   
<sup>2</sup>Royal Marsden NHS Foundation Trust, London, UK

#### Abstract
**Summary:** In this article we present Locus Explorer, a data visualisation and exploration tool for genetic association data. Locus Explorer is written in R using the Shiny library, providing access to powerful R-based functions through a simple user interface. Locus Explorer allows users to simultaneously display genetic, statistical and biological data in a single image and allows dynamic zooming and customisation of the plot features. Publication quality plots may be downloaded in a variety of file formats.   
**Availability and implementation:** Locus Explorer is open source and runs through R and a web browser. It is available at https://github.com/oncogenetics/LocusExplorer, where user guides and example data are also provided.

### Example plot output
<img src="www/chr17_36020000_36140000.jpeg" height="800px" width="600px" />


### Publications That Use Locus Explorer
[Multiple novel prostate cancer susceptibility signals identified by fine-mapping of known risk loci among Europeans.](http://www.ncbi.nlm.nih.gov/pubmed/26025378) Al Olama AA *et al.*

### Contact  
Questions, suggestions, and bug reports are welcome and appreciated. 
- Submit [suggestions and bug-reports](https://github.com/oncogenetics/LocusExplorer/issues)
- Send [pull request](https://github.com/oncogenetics/LocusExplorer/pulls)
- Contact email [T Dadaev](mailto: tokhir.dadaev@icr.ac.uk)

### To-do list
https://github.com/oncogenetics/LocusExplorer/issues   

[![Issue Stats](http://issuestats.com/github/oncogenetics/LocusExplorer/badge/pr)](http://issuestats.com/github/oncogenetics/LocusExplorer)
[![Issue Stats](http://issuestats.com/github/oncogenetics/LocusExplorer/badge/issue)](http://issuestats.com/github/oncogenetics/LocusExplorer)

