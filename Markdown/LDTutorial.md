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

Coming soon...



