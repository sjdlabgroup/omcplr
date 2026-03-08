#' Bin the mapping mapID contig with a resolution of 10kb.
#' This function takes genomic bin and xmap data and performs binning with a specified window size.
#' @param genome_bin A file path or data frame containing genomic bin information with columns "chr", "start", and "end".
#' @param xmap A file path or data frame containing xmap information with columns "mapID", "chr", "start_q", "end_q", "start", "end", and "strand".
#' @param window_size The size of the binning window in base pairs.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom data.table setDT
#' @export
overlap_bin_mapID <- function(genome_bin, xmap, window_size = 10000, ...) {
    library(GenomicRanges)
    library(data.table)
    
# read input data    
  if (exists("genome_bin", mode = "list")) {
    genome_bin <- genome_bin
  } else if (file.exists(genome_bin)) {
    genome_bin <-
      read.table(
        genome_bin,
        sep = "\t",
        header = F,
        stringsAsFactors = F
      )
  }
    colnames(genome_bin) <- c("chr", "start", "end")
    
  # read input data   
    if (exists("xmap", mode = "list")) {
      contig_xmap <- xmap
    } else if (file.exists(xmap)) {
      contig_xmap <-
        utils::read.table(
          xmap,
          sep = "\t",
          header = F,
          stringsAsFactors = F
        )
    }
    colnames(contig_xmap) <- c("mapID", "chr", "start_q", "end_q", "start", "end", "strand")
    contig_xmap <- contig_xmap[, c("chr", "start", "end", "mapID")]
    contig_xmap$chr <- paste0("chr", contig_xmap$chr)
    contig_xmap[which(contig_xmap$chr == "chr23"), "chr"] = "chrX"
    
    # identify overlaps between xmap position and reference bins 
    refGR <-
      makeGRangesFromDataFrame(genome_bin,
                               keep.extra.columns = TRUE,
                               ignore.strand = FALSE)
    contig_xmapGR <-
      makeGRangesFromDataFrame(contig_xmap,
                               keep.extra.columns = TRUE,
                               ignore.strand = FALSE)
    hits <- findOverlaps(refGR, contig_xmapGR)
    genome_bin_tmp <- contig_xmap[subjectHits(hits),]
    contig_xmap_tmp <- genome_bin[queryHits(hits),]
    contig_xmap_tmp$chr <- gsub("chr", "C", contig_xmap_tmp$chr)
    contig_xmap_tmp$combine <-
      paste(contig_xmap_tmp$chr,
            ceiling(contig_xmap_tmp$start / window_size),
            sep = "_")
    
    combine_tmp <- as.data.frame(cbind(contig_xmap_tmp$combine, genome_bin_tmp$mapID))
    combine_tmp <- unique(combine_tmp)
    combine_tmp <- data.table(combine_tmp)
#   combine_tmp[, `:=`(type_merge, paste0(as.character(V2), collapse = ",")),by = .(V1)]
    combine_tmp[, type_merge := paste0(as.character(V2), collapse = ","), by = .(V1)]
    combine_tmp <- unique(as.data.frame(combine_tmp)[, c("V1", "type_merge")])
    combine_tmp$V1 <- as.character(combine_tmp$V1)
    return(combine_tmp)
  }

