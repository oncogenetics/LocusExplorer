---
title: "Response to Reviewers"
author: "T.Dadaev"
date: "13 October 2015"
output: word_document
---

## Note

We have created a GitHub issue for each of the reviewer's comment, see the status of them using below links:   
- [Reviewer: 1 | status: 100% fixed](https://github.com/oncogenetics/LocusExplorer/issues?utf8=%E2%9C%93&q=milestone%3A%22Reviewer+1%22+)   
- [Reviewer: 2 | status: 100% fixed](https://github.com/oncogenetics/LocusExplorer/issues?utf8=%E2%9C%93&q=milestone%3A%22Reviewer+2%22+)   
- [Reviewer: 3 | status: 88% fixed](https://github.com/oncogenetics/LocusExplorer/issues?utf8=%E2%9C%93&q=milestone%3A%22Reviewer+3%22+) - Two issues are known and under development, will be fixed in future releases.   

If the provided solutions are not satisfactory, please suggest re-open using GitHub interface, [Re-open](https://github.com/oncogenetics/LocusExplorer/issues) existing issues or submit as a [New issue](https://github.com/oncogenetics/LocusExplorer/issues/new).

## Reviewer issues in detail

### Reviewer: 1 

1. I would like to suggest adding lines in the manhattan plot to show the key significant value, such as 0.01 or 0.05. This could be useful for understanding the significant associated markers.   

- [#42](https://github.com/oncogenetics/LocusExplorer/issues/42) - Added `sliderInput()` for suggestive line and genomewide line, red and blue lines, with default values of 5 and 8. Set to zero, to hide suggestive lines.

2. The color of the dots in the manhattan plot was not very bright. Also, the recombination ratio curve was not very bright.   

- [#43](https://github.com/oncogenetics/LocusExplorer/issues/43) - Changed recombination and SNP outline colours to 3 shades darker.

3. I would suggest the authors figuring out what species the tool was suitable.   

- [#44](https://github.com/oncogenetics/LocusExplorer/issues/44) - We are willing to provide help to anyone who wishes to adapt this application for non-human use, but we can't develop that aspect ourselves at the moment.

4. The Letters in the figure was not very clear. Also, it may overlap when the detail information of some significant markers was displayed together.   

- [#45](https://github.com/oncogenetics/LocusExplorer/issues/45) - Added SNP adjust option `sliderInput()`, using repulsion factor input value.

5. There were some mistakes when the format of JPEG changed to TIFF.   

- [#46](https://github.com/oncogenetics/LocusExplorer/issues/46) - We could not replicate this error. Tested on:

**Windows**
```R
R version 3.2.2 (2015-08-14)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1
```

**Mac OS**
```R
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.3 (Yosemite)
```

Please, suggest re-open this issue with your `sessionInfo()`.

### Reviewer: 2 

1. Minor revisions: "supplied BED file": "BED" --> "bedGraph", since the 4th column is numeric.      

- [#40](https://github.com/oncogenetics/LocusExplorer/issues/40) - BED file format definition is changed to bedGraph

2. Replicate Figure 1 in manuscript   
2.1. I tried to reproduce Figure 1 using Locus Explorer and noted a few differences. I suggest regenerating the image for Figure 1 using the latest version.   
In Figure 1, while the genomic region is labeled 36020000 - 36140000 along the bottom, it seems that the Manhattan plot covers that range while other plots such as SNPType see to extend beyond that range. In Locus Explorer, the SNPType track covers a smaller range (contains fewer elements) than in Figure 1, but now matches the range of the Manhattan plot.   
2.2. The LDSmooth curve shapes for rs718961 and rs11649743 are different -- in Figure 1 they have one peak and the pink curve is split at the top, while in Locus Explorer, they have a large peak plus a smaller peak to the right (depending on LD). I tried adjusting LD but none of the values (0, 0.05, 0.10, 0.15, and the default 0.20) seemed to duplicate the shapes in Figure 1. 

- [#47](https://github.com/oncogenetics/LocusExplorer/issues/47) - Updated manuscript has the latest, reproducible Figure 

### Reviewer: 3 

#### 1. Comments to the Author   

We suggest to submit the package to the Bioconductor platform where it will be reviewed for consistency and tested on different platforms. This will not only improve user-friendliness but also make sure that the tool keeps running after the next update of R and make it more visible and accessible to potential users.   

- We do not have required resources to upgrade this appliaction to a bioconductor package.

#### 2. Specifically, we suggest the following   

2.1. Major comments:   
2.1.1. When trying the software in RStudio, we initially receiced the following error after running   `
```R
runGitHub("LocusExplorer", "oncogenetics")
```
to start the tool: 
  
> Error in unloadNamespace(package) :   
>  namespace 'AnnotationDbi' is imported by 'geneplotter', 'VariantAnnotation', 'GenomicFeatures', 'annotate', 'genefilter', 'biomaRt', 'OrganismDbi', 'biovizBase' so cannot be unloaded   
> Error in library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc, :   
>  Package 'AnnotationDbi' version 1.31.18 cannot be unloaded 

- [#48](https://github.com/oncogenetics/LocusExplorer/issues/48) - Start the application with latest packages and in a new clean R session. `README.md` is updated.
                 
2.1.2. With a clean session, I receive the following error:   

> Note: the specification for S3 class "AsIs" in package 'jsonlite' seems equivalent to one from package 'BiocGenerics': not turning on duplicate class definitions for this class.    
> Warning in file(filename, "r", encoding = encoding) :cannot open file 'source/UDF.R': No such file or directory 
Error in file(filename, "r", encoding = encoding) : cannot open the connection 

- [#49](https://github.com/oncogenetics/LocusExplorer/issues/49) - Typo in the code fixed - *source* to *Source*.

2.1.3. When I click on plotSettings directly after loading the app, I receive the following message:    

> Error in nrow(plotDatGeneticMap()) : 
>  error in evaluating the argument 'x' in selecting a method for function 'nrow': Error in nrow(plotDatGeneticMap()) : 
  incorrect length (0), expecting: 1309    
> Error in nrow(plotDatLD()) :    
  error in evaluating the argument 'x' in selecting a method for function 'nrow': Error in fix.by(by.y, y) : 'by' must specify a uniquely valid column 

- [#50](https://github.com/oncogenetics/LocusExplorer/issues/50) - We get this error messages, as the `ggplot` is trying to plot, while the data is still being crunched. Updated the code to hide error messages.

2.1.4. When playing around with the region that should be displayed, one may receive error messages as: 
  
> chr2:172800000-NA missing value where TRUE/FALSE needed 

- [#51](https://github.com/oncogenetics/LocusExplorer/issues/51) - Solution covered by **2.1.3**.

2.1.5. Message "Please download Histone bigWig files" appears    

- [#52](https://github.com/oncogenetics/LocusExplorer/issues/52) - Added FAQ, see FAQ-Q5.   

2.1.6. Final plots warning   

> Cannot open the connection 

- [#53](https://github.com/oncogenetics/LocusExplorer/issues/53) -  Solution covered by **2.1.3**, and **2.1.2**.  


2.1.7. The suggested filenames when saving a plot contains " characters at the beginning and end 
The " should be removed because not all file systems may handle this appropriately.     

- [#54](https://github.com/oncogenetics/LocusExplorer/issues/54) - We could not replicate this error. Tested on:

**Windows**
```R
R version 3.2.2 (2015-08-14)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1
```

**Mac OS**
```R
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.3 (Yosemite)
```

Please, suggest re-open this issue with your `sessionInfo()`.
   

2.1.8. No file is produced after saving: On my machine, after saving the file to a specific folder, there is no output file being produced. I tried with in multiple folders with multiple different file names.   

- [#55](https://github.com/oncogenetics/LocusExplorer/issues/55) - Application was probably running within RStudio Viewer window. Application needs to run on an external web browser, e.g.: Chrome, IE, etc. If error persists, suggest to re-open this issue on [GitHub issue #55](https://github.com/oncogenetics/LocusExplorer/issues/55).

2.1.9. After closing and restarting the app, jumping directly to the Final plot, I receive the following message:    

>  Error in evaluating the argument 'x' in selecting a method for function 'nrow': Error in nrow(plotDatGeneticMap()) : 
  incorrect length (0), expecting: 1309 

- [#56](https://github.com/oncogenetics/LocusExplorer/issues/56) - Solution covered by **2.1.3**.

2.1.10. When selecting *Custom* data without uploading any custom data, multiple error messages appear in the different plot sections.    

> Error in evaluating the argument 'data' in selecting a method for function 'ggplot': Error: Please upload Association file 

- [#57](https://github.com/oncogenetics/LocusExplorer/issues/57) - This is a known issue, will be addressed in future releases.   

#### 2.2.Minor comments regarding the software:  
  
2.2.1. Plot settings, SNP filter, -log10 p-value. For the -log10 p value scale, it might make sense to display which p-value the transformed values correspond to. 

- [#58](https://github.com/oncogenetics/LocusExplorer/issues/58) - For *Manhattan plots* it is a common practice to show p-values on Y axis as `-log10(pvalue)` scale.

[Manhattan plots (wikipedia)](https://en.wikipedia.org/wiki/Manhattan_plot)   

> In GWAS Manhattan plots, genomic coordinates are displayed along the X-axis, with the negative logarithm of the association P-value for each single nucleotide polymorphism (SNP) displayed on the Y-axis, meaning that each dot on the Manhattan plot signifies a SNP. Because the strongest associations have the smallest P-values (e.g., 10^-15^), their negative logarithms will be the greatest (e.g., 15).

2.2.2. Custom tracks: Would it be possible and user-friendly to add the option to display custom tracks instead of the preselected ones? Maybe the functionality of the rtracklayer package can help?   

- [#59](https://github.com/oncogenetics/LocusExplorer/issues/59) - This is a known feature request, currently in development.     

2.2.3. *Make LD file* section: Initially I thought that the "Make LD file" section is needed after the final plot in addition somehow (following the menu from left to right: input data -> Plot settings -> final plot -> make LD file?). Now, after navigating through the app, it seems that this section is more suited for the help section to explain how one can obtain or generate a LD file if not readily available. I therefore suggest to move this to the help section because it seems to belong there and not next to the final plot. 

- [#60](https://github.com/oncogenetics/LocusExplorer/issues/60) - Moved to Help dropdown menu.

2.2.4. Custom input data: Is it really necessary to require a .txt file ending? This may be typical for Windows, but not so much for Linux or Mac. 

- [#61](https://github.com/oncogenetics/LocusExplorer/issues/61) - Any extension, any delimiter should work now - `data.table::fread()`. Updated InputFileFormat.md.

2.2.5. Width of menu in the help section: The width of the menu could be increased so the menu fits on one line. I have lots of free space in the right side, but the menu still spans two lines. 

- [#62](https://github.com/oncogenetics/LocusExplorer/issues/62) - Help tabs changed to dropdown menu.

2.2.6. Input data: BED file vs. bedgraph file: This is a bit confusing to the user because on the left side, the notation is "bedGraph", but on the right side, it says "BED file" 

- [#63](https://github.com/oncogenetics/LocusExplorer/issues/63) - Typo corrected, *BED* to *bedGraph*   

2.2.7. The color panel in the final plot section could be changed so that users can enter their own favourite color as RGB code or in the traditional HTML style (hexadecimal triplets).   

- [#64](https://github.com/oncogenetics/LocusExplorer/issues/64) - `shinyjs::colourInput` solution added. See, [shinyjs author's examples](http://deanattali.com/2015/06/28/introducing-shinyjs-colourinput/).

#### 2.3. Minor comments regarding the paper: 
  
2.3.1. "however, the statistically most significant GWAS variants only tag loci within the genome that contain the underlying functional variant(s) and are rarely causal themselves." You may want to highlight more prominently that often, these loci are non-coding 

- [#65](https://github.com/oncogenetics/LocusExplorer/issues/65) - Manuscript updated.

##### In the manuscript below lines: 
> however, the statistically most significant GWAS variants only tag loci within the genome that contain the underlying functional variant(s) and are rarely causal themselves.

##### are updated to:
> however, the statistically most significant GWAS variants only tag loci within the genome that contain the underlying functional variant(s), are rarely causal themselves and the association signals are frequently situated in non-coding regions.

