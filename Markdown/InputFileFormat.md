#### Input File Specifications and Format:  
- Input files must be delimited flat text files  
- Headers are required and should use the exact column names described below and used in the example input files   
- Example files are available at: `LocusExplorer/Data/CustomDataExample`  


**1. Association File**  
Association File is mandatory for plot generation. All other files are optional but enhance plot aesthetics and interpretation  
 - **CHR -** Chromosome on which variant is located preceded by "chr". e.g `chr2`, `chrX`
 - **SNP -** Variant ID. e.g. `rs12345`, `chr10:104329988:D` 
 - **BP	-** Start coordinate of variant (does not include chromosome or end coordinate for in/del variants). e.g. `104356185`
 - **P -** *P*-value for specified variant
 - **TYPED -** Use code `2` for typed and `1` for imputed variants

**2. LD File**  
LD File is not mandatory but is recommended for more informative plots. If user supplied LD data is not available, see **Make LD file** tab for instructions of how LD data relative to the index SNP(s) can be obtained from the 1000 Genomes Project Phase 3 Dataset.  
 - **CHR_A -** Chromosome on which Index SNP is located (n.b. do not include "chr"). e.g. `2`, `23`
 - **BP_A	-** Index SNP start coordinate (Hg19, do not include chromosome or end coordinate for in/del variants). e.g. `104356185`
 - **SNP_A -** Index SNP ID
 - **CHR_B -** Chromosome for SNP in LD with Index SNP (SNP_A)
 - **BP_B	-** Start coordinate (Hg19, do not include chromosome or end coordinate for in/del variants) of SNP in LD with Index SNP (SNP_A). e.g. `104315667`
 - **SNP_B -** ID of SNP in LD with Index SNP (SNP_A). e.g. `rs10786679`, `chr10:104329988:D` 
 - **R2 -** LD score between SNP_A and SNP_B (0 to 1). e.g. `1`, `0.740917`

*Note:* Lead SNP must be defined relative to itself for plotting purposes, e.g.:
```
CHR_A	BP_A	SNP_A	CHR_B	BP_B	SNP_B	R2
2	173309618	rs13410475	2	173309618	rs13410475	1
2	173309618	rs13410475	2	172827293	rs148800555	0.0906124
```
When using plink or LDlink method this does not need to be manually added.

**3. Custom bedGraph Track**  
The first four required bedGraph fields are:

`chrom` - The name of the chromosome (e.g. chr3, chrY).  
`chromStart` - The starting position of the feature in the chromosome.    
`chromEnd` - The ending position of the feature in the chromosome.  
`score` - A score, any number.    

See [BedGraph Track Format](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) for more details.

File is tab separated and has no header. This file will be used to create a bar chart. Score is the height, e.g.:

```
chr2	173292313	173371181	-100
chr2	173500000	173520000	1000
```

