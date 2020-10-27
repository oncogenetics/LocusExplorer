#### Q1. Can you tell me the naming convention for the genomic features?   
I don't know what uc031tcg.1 is (but I do understand PCAT1, I think that naming convention is RefSeq, are the others also?).

**A:** I am trying to re-build gene symbols by collapsing transcripts into genes, when transcripts do not overlap with gene symbols, they get named as transcript names - in this case something like uc031tcg.1 - this is UCSC ID.
See, [udf_GeneSymbol](https://github.com/oncogenetics/LocusExplorer/blob/master/Source/UDF.R) for details. This part of the script is quite heavy and we are [working](https://github.com/oncogenetics/LocusExplorer/issues/23) on it.

#### Q2: Can you tell me how to interpret the wavy lines that cross the plot (for text in the figure legend)?
So far I am using this sentence: "The colored lines spanning the plotting region indicate the extent of LD for the lead SNPs with the same color designation, where the height of the line represents __ and the length of the line represents ___." Can you send me a better sentence to describe how these lines should be interpreted if I am not on the right track here? 

**A:** It is a loess smoothing for matching hit SNP. If there are 2 SNPs marked with red and green shape and fill, then we will have 2 matching loess lines red and green. Smoothing is using LD values from 1000G phase EUR subset. We can use different cut-offs of LD: LD=0, LD > 0.1, LD >= 0.2, etc., usually LD = 0, i.e.: include all SNP LDs works best.

Y axis is 0 to 1, as in minimum and maximum value for LD - R^2. When *wavy lines* have similar shape, we can safely assume, those SNPs are the same signal. There is also an "R^2" track, the darker the lines the higher the LD, this track also helps visually see how hit SNPs overlap.

See [Manhattan.R](https://github.com/oncogenetics/LocusExplorer/blob/master/Source/Manhattan.R).

```R
# LD smooth per hit SNP - optional
if("LDSmooth" %in% input$ShowHideTracks) {
  gg_out <- gg_out +
    geom_smooth(data=plotDatLD(),aes(x=BP,y=R2_Adj,col=LDSmoothCol),
                method="loess",se=FALSE)}
```

#### Q3. Can you tell me what is included (maybe even just the data source) for the histone and DNase panels?
What the colour interpretation is for the histone panel?

**A:** Data from ENCODE project, see links for more info.

 - Histone - [LocusExplorer README.md](https://github.com/oncogenetics/LocusExplorer/tree/master/Data/wgEncodeBroadHistone), [UCSC tables](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?org=human&db=hg19&g=wgEncodeRegMarkH3k27ac)


 - DNaseI - [LocusExplorer README.md](https://github.com/oncogenetics/LocusExplorer/tree/master/Data/wgEncodeRegDnaseClustered), [UCSC tables](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?org=human&db=hg19&g=wgEncodeRegDnaseClustered)



#### Q4: Which operating system it is tested on?
**A:** See [RSessionInfo.md](https://github.com/oncogenetics/LocusExplorer/blob/master/Markdown/RSessionInfo.md)


#### Q5: "Please download Histone bigWig files" message on the plot?
**A:** As GitHub has limitations on size of the [repositories and files](https://help.github.com/articles/what-is-my-disk-quota/), Histone BigWig files are not included in `LocusExplorer/Data/EncodeBigWig/`. These files are public and can be downloaded from <a href="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k27ac/" target="_blank">UCSC golden path</a> - total ~2.5GB. Downloaded bigWig files must be saved in `LocusExplorer/Data/EncodeBigWig/` folder.

We are working on server version of LocusExplorer, expected to be **live by December 2015**. Keep an eye on https://github.com/oncogenetics/LocusExplorer page. This will resolve the issues of different R versions. R packages, and no limits on anntation data - such as bigWig files.


#### Q6: RStudio crashes when clicked on *Plot Settings*
After succesfully uploading data using *Input Data* tab, when clicked on *Plot Settings* RStudio crashes. This happens even when we run RGUI.

**A:** It is hard to guess the cause of the crash, it could be package dependencies with different versions. We can try to re-install packages.

#### !!! WARNING: Below steps will remove all of your installed packages !!!

Try below steps:

1. Run `.libPaths()` and get library path, choose the one where the user has a write access. e.g.: On my machine I will choose the first of those two listed folders.   

```R
> .libPaths()
[1] "C:/Users/tdadaev/Documents/R/win-library/3.2" "C:/Program Files/R/R-3.2.2/library"  

myLibraryLocation <- .libPaths()[1]
```

2. Make list of installed and required packages.

```R
# all packages excluding base
allPackages <- installed.packages()
allPackages <- allPackages[ is.na(allPackages[,4]), 1]

# CRAN packages required by this Application.
cranPackages <- c("shiny", "shinyjs", "data.table", "dplyr", "tidyr", 
                  "ggplot2", "knitr", "markdown", "stringr","DT","seqminer",
                  "lattice","cluster")
# Bioconductor packages required by this Application.
bioPackages <-  c("ggbio","GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene",
                  "org.Hs.eg.db","rtracklayer")
```

3. Remove all packages, excluding base.

```R
# remove all packages
sapply(allPackages, remove.packages)
```

4. Reinstall **only** required packages in `myLibraryLocation` folder.

```R
#reinstall CRAN packages
sapply(cranPackages, install.packages, lib = myLibraryLocation)

#reinstall Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("BiocInstaller")
biocLite(bioPackages)
```


#### Q7: Plink Error: Duplicate ID '.'
When creating LD file using "LD Tutorial", getting "Duplicate ID '.'" error.

**A:** Here is a guide on how to get LD files when we have duplicated ID issue in the 1000 genomes vcf files. I am using your example `snplist.txt` file

We need to get 1000 genomes vcf for the region, usually I give 250Kb flank based on my SNPs in `snplist.txt`.

```
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 17:1-71324000 > genotypes.vcf
```

This as you described give us an error.
```
plink --vcf  genotypes.vcf \
--r2 \
--ld-snp-list snplist.txt \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--out myLD

# ...
# Error: Duplicate ID '.'.
```

As some SNP names are actually dots `"."`, they don't have an RS ID attached to them. So we need to give them unique names.

Convert to plink format   

```
plink --vcf  genotypes.vcf \
--make-bed \
--out genotypes
```

Here is the cause of the problem, dots as SNP IDs:

```
head -n5 genotypes.bim
# 17      .       0       56      T       C
# 17      .       0       78      C       G
# 17      .       0       355     A       G
# 17      .       0       684     T       C
# 17      rs62053745      0       828     T       C
```
Run [this R script](https://github.com/oncogenetics/Miscellaneous/blob/master/convert/makeSNPnamesUnique.R) to make unique names:
```
Rscript makeSNPnamesUnique.R genotypes.bim 
head -n5 genotypes.bim 
# 17      17_56_T_C       0       56      T       C
# 17      17_78_C_G       0       78      C       G
# 17      17_355_A_G      0       355     A       G
# 17      17_684_T_C      0       684     T       C
# 17      rs62053745      0       828     T       C
```
Now, the SNP  IDs are fixed, we run plink LD as usual:
```
plink --bfile genotypes \
--r2 \
--ld-snp-list snplist.txt \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \                   # to have smaller output file: --ld-window-r2 0.2
--out myLD
```

**Note:** If the output file is too big, then consider setting higher R2 filter value, for example `--ld-window-r2 0.2`, meaning only output R2 values that are more than `0.2`.

Check the output
```
wc -l myLD.ld
# 49172 myLD.ld

head -n5 myLD.ld
#  CHR_A         BP_A             SNP_A  CHR_B         BP_B             SNP_B           R2
#     17          834         rs9747082     17           56         17_56_T_C   0.00137027
#     17          834         rs9747082     17          355        17_355_A_G   0.00151223
#     17          834         rs9747082     17          684        17_684_T_C   0.00127126
#     17          834         rs9747082     17          828        rs62053745     0.678207
```

Finally, we can use this file as input for LocusExplorer.

