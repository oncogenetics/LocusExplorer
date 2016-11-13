# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt



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





# Data for forceNetwork() -------------------------------------------------
#Hits[method] + LD plots
udf_NodesAndLinks <- function(method=NA){
  
  #method <- "Stepwise Forward"
  # datLD <- fread("Data/ProstateData/OncoArray/LE/chr17_35547276_36603565_LD.txt",
  #                header=TRUE, data.table=FALSE)
  # 
  # methodStats <- fread("Data/ProstateData/OncoArray/LE/chr17_35547276_36603565_stats.txt",
  #                      header=TRUE, data.table=FALSE) %>% 
  #   filter(Method == method)
  
  datLD <- datLD()
  
  methodStats <- datHitSNPStats() %>% 
    filter(Method == method)
  
  
  methodHits <- methodStats$SNP
  
  #make all combos for hits, in case there is no LD, then we set it to 0
  # so they show up on graph plot even if they are not connected to anything.
  methodHitsCombn <-
    as.data.frame(t(combn(methodHits,2)))
  colnames(methodHitsCombn) <- c("SNP_A","SNP_B")
  methodHitsCombn$SNP_A <- as.character(methodHitsCombn$SNP_A)
  methodHitsCombn$SNP_B <- as.character(methodHitsCombn$SNP_B)
  
  datGraph <- merge(
    datLD[,c("SNP_A","SNP_B","R2")] %>% 
      filter(SNP_A!=SNP_B & 
               R2 >=  input$FilterMinLD &
               (SNP_A %in% methodHits |
                  SNP_B %in% methodHits)) %>% 
      group_by(SNP_A) %>% 
      filter(n()>1),
    methodHitsCombn, 
    all=TRUE)
  #set missing R2 to zero - 0
  datGraph$R2 <- ifelse(is.na(datGraph$R2),0,datGraph$R2)
  
  #prepare graph plot data
  # Nodes and Links
  MisNodes <- data.table(name=unique(c(methodHits,
                                       datGraph$SNP_A,datGraph$SNP_B)),
                         #id must be zero based
                         id=0:(length(unique(c(methodHits,
                                               datGraph$SNP_A,datGraph$SNP_B)))-1))
  MisLinks <- merge(datGraph,MisNodes,by.x="SNP_A",by.y="name",all.x=TRUE)
  MisLinks <- merge(MisLinks,MisNodes,by.x="SNP_B",by.y="name",all.x=TRUE)
  MisLinks <- MisLinks %>% 
    transmute(source=ifelse(id.x<=id.y,id.x,id.y),
              target=ifelse(id.y>id.x,id.y,id.x),
              value=R2) %>% unique
  #group nodes
  MisNodes <- MisNodes %>% 
    mutate(group=ifelse(name %in% methodHits,"Hits","Other")) 
  
  #get Pvalues, use to set the size of nodes
  MisNodes <- merge(MisNodes,
                    methodStats[,c("SNP","Stats")],
                    by.x="name",by.y="SNP",all.x=TRUE) %>% unique
  MisNodes$size <- abs(ifelse(is.na(MisNodes$Stats),1,
                          round(log10(abs(MisNodes$Stats))))
                       ) + 5
  MisNodes$size <- ifelse(MisNodes$size > 30, 30, MisNodes$size)
  MisNodes <- MisNodes %>% arrange(id)
  
  #return list
  return(list(MisNodes=MisNodes,
              MisLinks=MisLinks))
  
}

# Data for arcplot() -------------------------------------------------
#Hits[method] + LD plots
udf_arcPlotData <- function(method=NA){
  #testing
  #method <- "Stepwise Forward"
  # datLD <- fread("Data/ProstateData/OncoArray/LE/chr17_35547276_36603565_LD.txt",
  #                header=TRUE, data.table=FALSE)
  # datHitSNPStats <- fread("Data/ProstateData/OncoArray/LE/chr17_35547276_36603565_stats.txt",
  #                         header=TRUE, data.table=FALSE) 
  # RegionChrN <- 17
  
  
  #input data
  datLD <- datLD()
  datHitSNPStats <- datHitSNPStats() 
  RegionChrN <- RegionChrN()
  
  hits <- unique(datHitSNPStats$SNP)
  methodHits <- datHitSNPStats %>% filter(Method %in% method) %>% .$SNP
  
  # Data 
  datArc <- datLD %>% 
    filter(
      SNP_A == SNP_B | (
        SNP_B %in% hits &
          R2 >= input$FilterMinLD
      )) 
  #if hits are not in LD file then add
  datMissingHits <-
    datHitSNPStats %>% 
    transmute(CHR_A=RegionChrN,
              BP_A=BP,
              SNP_A=SNP,
              CHR_B=RegionChrN,
              BP_B=BP,
              SNP_B=SNP,
              R2=1)
  datArc <- rbind(datArc,datMissingHits) %>% 
    unique %>% 
    arrange(BP_B,BP_A)
  #keep R2 for methodhits only, the rest 0
  datArc$R2 <- ifelse(datArc$SNP_B %in% methodHits,datArc$R2,0)
  
  #Data for arcplot
  # edgelist
  lab <- cbind(datArc$SNP_A,datArc$SNP_B)
  # make graph from edgelist
  glab <- graph.edgelist(lab, directed = F)
  
  E(glab)$weight <- datArc$R2
  
  #reorder based on BP - position of SNPs
  myV <- data.frame(SNP=as.character(names(V(glab))),
                    stringsAsFactors = FALSE)
  myV$rn <- 1:nrow(myV)
  
  orderV <- left_join(myV,
                      datHitSNPStats[,c("SNP","BP")],by="SNP") %>%
    unique %>% arrange(BP) %>% 
    mutate(orderV=row_number()) %>% 
    arrange(rn) %>% .$orderV
  
  # get clusters and rename to hits and other SNPs
  gclus <- clusters(glab)
  gclus$membership <- 
    ifelse(names(gclus$membership) %in% methodHits,1,2)
  # colour clusters
  cols <- c("green","grey50")[gclus$membership]
  
  return(list(lab=lab,
              glab=glab,
              cols=cols,
              orderV=orderV))
  }
  