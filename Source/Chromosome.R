Ideogram(genome="hg19",
         subchr = RegionChr(), color="#F00034", fill="#F00034") +
  #highlight ROI - region
  xlim(IRanges(RegionStart(),RegionEnd())) +
  ylab(NULL)
