**Note -** The maximum region of interest that Locus Explorer can display is 2Mb  

#### Input File Specifications and Format:  
- Input files must be **tab-delimited** and saved as **.txt**
- Headers are required and should use the exact column names described below and used in the example input files
- Example files are available at: *oncogenetics//LocusExplorer//Data//CustomDataExample*  


**1. Association File - _Association.txt_**  
*Note* - Association File is mandatory for plot generation. All other files are optional but enhance plot aesthetics and interpretation  
 - **CHR -** Chromosome on which variant is located preceded by "chr". e.g `chr2`, `chrX`
 - **SNP -** Variant ID. e.g. `rs12345`, `chr10:104329988:D` 
 - **BP	-** Start coordinate of variant (does not include chromosome or end coordinate for in/del variants). e.g. `104356185`
 - **P -** *P*-value for specified variant
 - **TYPED -** Use code `2` for typed and `1` for imputed variants

**2. LD File - _LD.txt_**  
*Note* - LD File is not mandatory but is recommended for more informative plots. If user supplied LD data is not available, see **Make LD file** tab for instructions of how LD data relative to the index SNP(s) can be obtained from the 1000 Genomes Project Phase 3 Dataset.  
 - **CHR_A -** Chromosome on which Index SNP is located (n.b. do not include "chr"). e.g. `2`, `23`
 - **BP_A	-** Index SNP start coordinate (Hg19, do not include chromosome or end coordinate for in/del variants). e.g. `104356185`
 - **SNP_A -** Index SNP ID
 - **CHR_B -** Chromosome for SNP in LD with Index SNP (SNP_A)
 - **BP_B	-** Start coordinate (Hg19, do not include chromosome or end coordinate for in/del variants) of SNP in LD with Index SNP (SNP_A). e.g. `104315667`
 - **SNP_B -** ID of SNP in LD with Index SNP (SNP_A). e.g. `rs10786679`, `chr10:104329988:D` 
 - **R2 -** LD score between SNP_A and SNP_B (0 to 1). e.g. `1`, `0.740917`

**3. Biofeature BED Track - _LNCAP.txt_**  
*Note* - Biofeature annotation file is not mandatory but can be used to display the positions of regulatory regions of interest from user supplied or publically available datasets. Visit xxxx for instructions of how to obtain data of interest from the ENCODE Project if custom data is unavailable  
 - **CHR -** Chromosome on which Biofeature is located preceded by "chr". e.g `chr2`, `chrX`
 - **BP	-** Position of Biofeature
 - **ANNOTATE -** Description of Biofeature. e.g. `exonic`, `intronic`, `UTR3`, etc.
 
**4. Gene Expression BED Track - _eQTL.txt_** 
*Note* - Gene expression File is not mandatory but can be used to display the positions of genes differentially regulated by variants within the plot region. Genes within the plot region for which no eQTL is odserved should be excluded from the input file. Visit xxxx for instructions of how to obtain eQTL data from the TCGA or GTEx datasets if custom data is unavailable
 - **CHR -** Chromosome on which differentially expressed gene is located preceded by "chr". e.g `chr2` `chrX`
 - **START -** Transcript start position
 - **END -** Transcript end position
 - **DIRECTION -** Use `1` for Upregulated genes and `-1` for Downregulated genes


