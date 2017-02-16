# Writing an R package from scratch
# https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/


# Step 0: Packages you will need
library(devtools)
library(roxygen2)

# Step 1: Create your package directory
#setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects")
#create("oncofunco")

# Step 2: Add functions
# save them in R folder
# Step 3: Add documentation

# Step 4: Process your documentation
library(devtools)
library(roxygen2)
#setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/oncofunco/")
setwd("N:/Translational Cancer Genetics Team/Bioinformatics/Development/R_Packages/oncofunco")
document()

# Step 5: Install!
setwd("..")
install("oncofunco")

# check how Rd files look
library(oncofunco)
glmSummary()
getMAF()
oncofunco::colourHue(3)
oncofunco::strPadLeft("dd")
oncofunco::theme_LE()
oncofunco::allMatch()
oncofunco::cleanPackages()
oncofunco::gen2dose()
oncofunco::gen2vcf()
oncofunco::ICRcolours(plot = TRUE)
oncofunco::plotLD()
barplot(rep(1, 11), border = NA,
        col = oncofunco::shadeTintColour(col = "orange",
                                         change = seq(0, 100, 10)))

# install on hpc ----------------------------------------------------------
library("devtools",lib.loc="/scratch/cancgene/tdadaev/Rpackages/")
library("httr",lib.loc="/scratch/cancgene/tdadaev/Rpackages/")
library("curl",lib.loc="/scratch/cancgene/tdadaev/Rpackages/")
library("withr",lib.loc="/scratch/cancgene/tdadaev/Rpackages/")

with_libpaths("/scratch/cancgene/tdadaev/Rpackages/",install_github("oncogenetics/oncofunco"))
