#### Input File Specifications and Format:  
- Input files must be **tab-delimited** and saved as **.txt**
- Headers are required and should use the exact column names described below and used in the example input files
- Example files are available at: *oncogenetics//LocusExplorer//Data//CustomDataExample*  


**1. Association File - _Association.txt_**  
Association File is mandatory for plot generation. All other files are optional but enhance plot aesthetics and interpretation  
 - **CHR -** Chromosome on which variant is located preceded by "chr". e.g `chr2`, `chrX`
 - **SNP -** Variant ID. e.g. `rs12345`, `chr10:104329988:D` 
 - **BP	-** Start coordinate of variant (does not include chromosome or end coordinate for in/del variants). e.g. `104356185`
 - **P -** *P*-value for specified variant
 - **TYPED -** Use code `2` for typed and `1` for imputed variants

**2. LD File - _LD.txt_**  
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

**3. Gene Expression BED Track - _eQTL.txt_**  
Gene expression File is not mandatory but can be used to display the positions of genes differentially regulated by variants within the plot region. Genes within the plot region for which no eQTL is odserved should be excluded from the input file. 
 - **CHR -** Chromosome on which differentially expressed gene is located preceded by "chr". e.g `chr2` `chrX`
 - **START -** Transcript start position
 - **END -** Transcript end position
 - **DIRECTION -** Use `1` for Upregulated genes and `-1` for Downregulated genes
