LD Tutorial
=============

### Make LD file from 1000 Genomes Data

#### 1. Using LDlink website

1. Go to [LDlink](http://analysistools.nci.nih.gov/LDlink/) website.
2. Select **LDproxy** tab.
3. Enter SNP rs number (must be in 1000 Genomes phase 3 dataset / dbSNP 142).
4. Select population(s) for LD calculation and hit **Calculate**.
5. Scroll down and right click **Download all proxy SNPs** and click **Save link as...** to save as a text file.
6. Return to the **Make LD file** tab of Locus Explorer and upload unprocessed LDlink file as input.
7. Download processed LD File for use at **Input Data** tab in Locus Explorer.

**Note:**

* This procedure will calculate LD relative to one top index SNP only. For regions with multiple index SNPs, the process can be performed separately for each individual index SNP and the individual processed LD files combined to make a single input LD file.

* LD relative to the index SNP cannot be calculated for variants that are not present in the 1000 Genomes phase 3 dataset (dbSNP 142) from publically available data and these will therefore always appear to be uncorrelated on the plot when generating LD information in this way.

#### 2. Using tabix and plink
We will need to install tabix and plink 1.9.
Download and install [HTSlib package](http://www.htslib.org/download/) which will include tabix. Then download and install [PLINK 1.9](https://www.cog-genomics.org/plink2).

Then we download 1000 Genomes VCF file using tabix and calculate LD using PLINK.

Example:

##### 2.1 Download 1000 Genomes VCF
Download vcf for region of interest `16:56995835-57017756` from 1000 genomes ftp site using tabix.  
`tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 16:56995835-57017756 > genotypes.vcf`

##### 2.1 Use plink to calculate LD
###### 2.1.1 Calcualte LD for 2 SNPs  
```
plink --vcf genotypes.vcf --ld rs9935228 rs1864163

#output
--ld rs9935228 rs1864163:

   R-sq = 7.20831e-05    D' = 0.0355584

   Haplotype     Frequency    Expectation under LE
   ---------     ---------    --------------------
          GA      0.005359                0.004853
          AA      0.249013                0.249519
          GG      0.013719                0.014225
          AG      0.731909                0.731403

   In phase alleles are GA/AG
```
###### 2.1.2 Calculate LD for list of SNP against all SNPs within 1000kb region.

We will need a list of SNPs file, one SNP per row.
```
# Example file snplist.txt
> cat snplist.txt
# rs9935228
# rs1864163

# Now we pass snplist to plink. To learn more about plink options selected below see XYZ.
plink --vcf  genotypes.vcf \
--r2 \
--ld-snp-list snplist.txt \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--out LD_rs9935228_rs1864163
```

