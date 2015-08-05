 - Region of interest must be no more than 1-2Mb.  
 - Input files must be tab delimited.

#### Example input files:  
Files are at: *LocusExplorer\\Data\\CustomDataExample*    


1. Association File - *Association.txt*   
 - CHR - chromosome name, chr2
 - SNP - variant name, string, rs12345
 - BP	- base pair position, numeric
 - P - pvalue, numeric
 - TYPED - if the SNP was typed (2) or imputed (1) 

2. LD - *LD.txt*  
 - CHR_A - chromosome name, chr2
 - BP_A	- hit SNP position
 - SNP_A - hit SNP name
 - CHR_B - chromosome of SNP in LD with hit SNP
 - BP_B	- position of SNP in LD with hit SNPSNP_B
 - R2 - LD score, numeric, 0 to 1
- *Note:* Online resources to extract LD: 
[Biostars: 1000 Genomes Ld Calculation](https://www.biostars.org/p/2909/), [LDlink](http://analysistools.nci.nih.gov/LDlink/), [Haploreg v3](http://www.broadinstitute.org/mammals/haploreg/haploreg_v3.php), [plink](http://pngu.mgh.harvard.edu/~purcell/plink/ld.shtml), [Ensemble](http://www.ensembl.info/blog/2015/06/18/1000-genomes-phase-3-frequencies-genotypes-and-ld-data/)

3. LNCAP - *LNCAP.txt*  
 - CHR - chromosome name, chr2
 - BP	- position, numeric
 - ANNOTATE - text string: exonic, intronic, UTR3, etc.
 
4. eQTL - *eQTL.txt*  
Gene expression data.
 - CHR - chromosome name, chr2
 - START - position, transcript start position
 - END - position, trabscript end position
 - DIRECTION - overexpressed 1, underexpressed -1


