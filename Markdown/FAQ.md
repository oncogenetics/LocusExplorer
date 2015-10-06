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

 - Histone - [LocusExplorer README.md](https://github.com/oncogenetics/LocusExplorer/tree/master/Data/wgEncodeBroadHistone), [UCSC tables](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=209730420_hvKa8YLni4ekR6rt9w1cSuGABF1Y&c=chr7&g=wgEncodeRegMarkH3k27ac)

 - DNaseI - [LocusExplorer README.md](https://github.com/oncogenetics/LocusExplorer/tree/master/Data/wgEncodeRegDnaseClustered), [UCSC tables](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=209730420_hvKa8YLni4ekR6rt9w1cSuGABF1Y&c=chr7&g=wgEncodeRegDnaseClustered)

#### Q4: Which operating system it is tested on?
**A:** See [RSessionInfo.md](https://github.com/oncogenetics/LocusExplorer/blob/master/Markdown/RSessionInfo.md)


#### Q5: "Please download Histone bigWig files" message on the plot?
As GitHub has limitations on size of the [repositories and files](https://help.github.com/articles/what-is-my-disk-quota/), Histone BigWig files are not included in `LocusExplorer/Data/EncodeBigWig/`. These files are public and can be downloaded from <a href="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k27ac/" target="_blank">UCSC golden path</a> - total ~2.5GB. Downloaded bigWig files must be saved in `LocusExplorer/Data/EncodeBigWig/` folder.

We are working on server version of LocusExplorer, expected to be **live by December 2015**. Keep an eye on https://github.com/oncogenetics/LocusExplorer page. This will resolve the issues of different R versions. R packages, and no limits on anntation data - such as bigWig files.


