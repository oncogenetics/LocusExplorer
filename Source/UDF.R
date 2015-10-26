# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# General options for all plots -------------------------------------------
# ggplot() + udf_theme()
udf_theme <- function(){
    theme(legend.position="none",
          panel.background = element_rect(fill="white"),
          panel.grid.minor=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey60", linetype = "dotted"),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          #Y Axis font
          axis.text.y=element_text(family="Courier",colour = "grey20")
    )
}

# Padding function, used to add custom labels for Yaxis, to have them aligned.
udf_pad <- function(labels){
  require(stringr)
  return(str_pad(labels,15,"left",pad = " "))
}


# GeneSymbol --------------------------------------------------------------
# https://github.com/oncogenetics/R_UDF/blob/master/GeneSymbol.R
udf_GeneSymbol <- function(chrom=NA,chromStart=NA,chromEnd=NA){
  require(dplyr)
  require(ggplot2)
  require(ggbio)
  require(GenomicFeatures)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db) # gene symobols
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  #valid chrom names
  chr <- paste0("chr",c(1:22,"X","Y"))
  
  # Input checks ----------------------------------------------------------
  if(!is.character(chrom) |
     length(chrom)!=1 |
     !chrom %in% chr) stop("chrom: must be character class with length(chrom)==1, e.g.: chr1")
  
  if(is.na(chromStart) |
     chromStart < 0) warning("chromStart: setting to default value 0")
  
  if(is.na(chromEnd)) warning("chromEnd: default value last position for chromosome")
  
  # Validate start end position for chrom - using txdb chrominfo table
  chrominfo <- 
    dbGetQuery(txdb$conn,
               paste0("select * from chrominfo where chrom='",
                      chrom,"'"))
  chromStart <- ifelse(is.na(chromStart) | 
                         chromStart < 0,
                       0,chromStart)
  chromEnd <- ifelse(is.na(chromEnd) | 
                       chromEnd > chrominfo$length |
                       chromEnd < chromStart,
                     chrominfo$length, chromEnd)
  
  #Print summary for selection region of interest
  print(paste0("Collapsing to gene symbols, region: ",
               chrom,":",chromStart,"-",chromEnd))
  
  # Get chromosome start end positions to subset TXDB for transcripts
  roi_chr <- GRanges(seqnames=chrominfo$chrom,
                     IRanges(start=chromStart,
                             end=chromEnd,
                             names=chrom))
  # Collapse to gene symbol -----------------------------------------------
  # Subset txdb over overlaps for chr-start-end
  keys_overlap_tx_id <- 
    subsetByOverlaps(transcripts(txdb),roi_chr) %>% 
    as.data.frame() %>% 
    .$tx_id %>% as.character
  
  # Match TX ID to GENEID
  TXID_GENEID <- AnnotationDbi::select(txdb,
                                       keys=keys_overlap_tx_id,
                                       columns=c("TXNAME","GENEID"),
                                       keytype="TXID")
  
  # Select transcipts from txdb which have GENEID
  Trans <- AnnotationDbi::select(txdb,
                                 keys=as.character(TXID_GENEID$TXID),
                                 columns=columns(txdb),
                                 keytype = "TXID")
  
  # Get gene symbol
  gene_symbol <- unique(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys=Trans[ !is.na(Trans$GENEID),"GENEID"],
                          columns="SYMBOL",
                          keytype="ENTREZID"))
  
  # Match GENEID, SYMBOL
  TXID_GENEID <- left_join(TXID_GENEID,gene_symbol,
                           by=c("GENEID" = "ENTREZID"))
  
  # If not match on gene symbol, then TXNAME is gene symbol
  TXID_GENEID$SYMBOL <- ifelse(is.na(TXID_GENEID$SYMBOL),
                               TXID_GENEID$TXNAME,TXID_GENEID$SYMBOL)
  
  # merge to add gene symbol
  Trans <- left_join(Trans,TXID_GENEID,
                     by=c("TXID","GENEID","TXNAME"))
  
  # If not match on gene symbol, then TXNAME is gene symbol
  Trans$SYMBOL <- ifelse(is.na(Trans$SYMBOL),
                         Trans$TXNAME,Trans$SYMBOL)
  
  #Make Granges object 
  CollapsedGenes <- 
    GRanges(seqnames=Trans$EXONCHROM,
            IRanges(start=Trans$EXONSTART,
                    end=Trans$EXONEND),
            strand=Trans$EXONSTRAND)
  CollapsedGenes$gene_id <- Trans$SYMBOL
  
  # pad gene names to align
  CollapsedGenes@elementMetadata$gene_id <- 
    udf_pad(CollapsedGenes@elementMetadata$gene_id)
  
  # Output ----------------------------------------------------------------
  #return collapsed genes per CHR
  return(CollapsedGenes)
}


# SNP label reposition ----------------------------------------------------
# Adapted form FField package by Grigori Kapoustin
# https://cran.r-project.org/web/packages/FField/
FFieldPtRep <- function (coords, rep.fact = 20, rep.dist.lmt = 10, attr.fact = 0.2, 
                         adj.max = 0.1, adj.lmt = 0.5, iter.max = 10000) 
{
  if (length(dim(coords)) != 2) {
    stop("FFieldPtRep: dim(coords) must be 2\n")
  }
  if (ncol(coords) < 2) {
    stop("FFieldPtRep: ncol(coords) must be >= 2\n")
  }
  if (nrow(coords) < 2) {
    stop("FFieldPtRep: nrow(coords) must be >= 2\n")
  }
  coords <- as.data.frame(coords)
  colnames(coords)[(1:2)] <- c("x", "y")
  coords.orig <- coords
  FVCalc <- function(vects.x, vects.y, f.fact, f.type = "invsq") {
    d.sq <- (vects.x^2 + vects.y^2)
    d <- sqrt(d.sq)
    vects.x <- vects.x/d
    vects.y <- vects.y/d
    if (f.type == "invsq") {
      d.sq[d >= rep.dist.lmt] <- Inf
      vect.f.x <- vects.x/d.sq * f.fact
      vect.f.y <- vects.y/d.sq * f.fact
    }
    else if (f.type == "lin") {
      vect.f.x <- vects.x * d * f.fact
      vect.f.y <- vects.y * d * f.fact
    }
    else {
      stop("FFieldPtRep: Unexpected f.type\n")
    }
    vect.f.x[is.na(vect.f.x)] <- 0
    vect.f.y[is.na(vect.f.y)] <- 0
    f.vect <- cbind(colSums(vect.f.x), colSums(vect.f.y))
    return(f.vect)
  }
  iter <- 0
  repeat {
    vects.x <- apply(coords, 1, function(c) (c[1] - coords$x))
    vects.y <- apply(coords, 1, function(c) (c[2] - coords$y))
    f.rep.v <- FVCalc(vects.x = vects.x, vects.y = vects.y, 
                      f.fact = rep.fact, f.type = "invsq")
    vects.orig <- coords.orig - coords
    f.attr.v <- FVCalc(vects.x = t(as.matrix(vects.orig$x)), 
                       vects.y = t(as.matrix(vects.orig$y)), f.fact = attr.fact, 
                       f.type = "lin")
    f.v <- f.rep.v + f.attr.v
    if (all(abs(f.v) <= adj.lmt)) {
      (break)()
    }
    mv.vect <- apply(f.v, c(1, 2), function(x) sign(x) * 
                       min(abs(x), adj.max))
    coords <- coords + mv.vect
    if ((iter <- iter + 1) > iter.max) {
      warning("FFieldPtRep: Maximum iterations exceeded ", 
              "without convergence.\n")
      (break)()
    }
  }
  return(coords)
}


# Shade tint HEX input ----------------------------------------------------
udf_shadeTintColor <- function(color, change = 25) {  
  # Example: shadeColor("#5CFF5C", 25)
  # positive shade
  # negative tints
  RGB <- col2rgb(color)
  RGB <- RGB - RGB/100*change
  RGB <- ifelse(RGB < 0, 0, RGB)
  RGB <- ifelse(RGB > 255, 255, RGB)
  
  return(rgb(RGB[1],RGB[2],RGB[3],maxColorValue = 255))
}
