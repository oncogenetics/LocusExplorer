library(ggplot2)
library(Rcpp)
library(ggbio)


plotDatGene <- 
  udf_GeneSymbol(chrom="chr17",
                 chromStart=41101591,
                 chromEnd=41372220) 


plotDatGeneN <- length(unique(plotDatGene@elementMetadata$gene_id))


length(plotDatGene@elementMetadata$gene_id)
length(as.character(strand(plotDatGene)))


plotDatGene@elementMetadata$strand1 <- as.character(strand(plotDatGene))
table(plotDatGene@elementMetadata$strand1)
table(plotDatGene@seqnames)


ggplot() + 
  geom_alignment(data=plotDatGene,aes(group=gene_id,fill=strand, col=strand),
                 facets=strand~.)
ggplot() + 
  geom_alignment(data=plotDatGene,aes(group=gene_id,fill=strand, col=strand)) +
  facet_grid(strand~.)

ggplot() + 
  geom_alignment(data=plotDatGene,aes(group=gene_id,fill=strand, col=strand))


  #general options
  ylab("xxx")  +
  facet_grid(strand~.,scales="free",space="free") +
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red")
#  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() 

ggplot() + 
  #geom_hline(yintercept=c(1:plotDatGeneN),col="grey80",linetype="dashed") +
  #ggbio
  geom_alignment(data = plotDatGene[plotDatGene@strand=="-",],
                 aes(group=gene_id,fill=strand, col=strand))

ggplot() + 
  #geom_hline(yintercept=c(1:plotDatGeneN),col="grey80",linetype="dashed") +
  #ggbio
  geom_alignment(data = c(plotDatGene[plotDatGene@strand=="+",],
                          plotDatGene[plotDatGene@strand=="-",]),
                 aes(group=gene_id,fill=strand, col=strand))


plotDatGene@elementMetadata$gene_id

strand(plotDatGene)
head(as.data.frame(plotDatGene))


#example
set.seed(1)
N=10
gr <- GRanges(seqnames = 
                sample(c("chr1", "chr1", "chr1"),
                       size = N, replace = TRUE),
              IRanges(
                start = sample(1:300, size = N, replace = TRUE),
                width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"), 
                              size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                            replace = TRUE))
ggplot(gr) + geom_alignment(facets = sample ~ strand)

ggplot(gr) + geom_alignment(aes(group = pair),facets = sample ~ seqnames)
