# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.19
# Function: Retrieve signal records from BwgTrack objects
# Input: BwgTrack object and genomic region
# Output: data.table, BwgTrack, or GRanges

#' Retrieve signal records from a BwgTrack object
#'
#' @description
#' `retrieve_bwg()` is the unified signal retrieval API in GeneTrackR. It
#' retrieves signal intervals from bedGraph, wig, or bigWig tracks by genomic
#' region. It replaces the older `query_bwg()` interface.
#'
#' @param object A BwgTrack object returned by `read_bwg()`.
#' @param chrom Chromosome name.
#' @param start Region start in 1-based closed coordinates.
#' @param end Region end in 1-based closed coordinates.
#' @param samples Optional character vector of sample IDs to retrieve. If NULL,
#' all samples are queried.
#' @param strand Strand selector. Use `"ignore"` for unstranded retrieval,
#' `"+"` or `"-"` for strand-specific tracks, `"both"` for both strands, or
#' `"auto"` when called by higher-level gene-aware plotting functions.
#' @param strand_policy How to handle strand filtering for unstranded signal
#' files. `"ignore_unstranded"` returns unstranded tracks for any strand request;
#' `"strict"` only returns explicitly matching strand records.
#' @param as Output type. `"data.table"` returns a signal table, `"BwgTrack"`
#' returns an in-memory BwgTrack sub-object, and `"GRanges"` returns a GRanges
#' object.
#' @param verbose Logical. Whether to print progress messages.
#' @param progress Logical. Whether to show a text progress bar for multi-sample
#' queries.
#' @param keep_empty_samples Logical. If TRUE, samples without signal in the
#' requested region are returned as zero-valued placeholder intervals. This is
#' useful for preserving empty facets in plots.
#' @param tabix_empty_fallback Optional logical. If TRUE, empty tabix results are
#' verified by full-file fread. Keep FALSE for speed unless debugging tabix
#' coordinate/index issues.
#'
#' @return A data.table, BwgTrack object, or GRanges object.
#'
#' @examples
#' \dontrun{
#' bg <- read_bwg("sample.bedgraph.gz", format = "bedgraph", mode = "lazy")
#'
#' dt <- retrieve_bwg(
#'   bg,
#'   chrom = "chr1",
#'   start = 1,
#'   end = 10000
#' )
#'
#' bg_sub <- retrieve_bwg(
#'   bg,
#'   chrom = "chr1",
#'   start = 1,
#'   end = 10000,
#'   as = "BwgTrack"
#' )
#' }
#'
#' @export
retrieve_bwg <- function(object,
                         chrom,
                         start,
                         end,
                         samples = NULL,
                         strand = c("ignore", "+", "-", "both", "auto"),
                         strand_policy = c("ignore_unstranded", "strict"),
                         as = c("data.table", "BwgTrack", "GRanges"),
                         verbose = FALSE,
                         progress = interactive() && isTRUE(verbose),
                         keep_empty_samples = FALSE,
                         tabix_empty_fallback = NULL) {
  as <- match.arg(as)
  strand <- match.arg(strand)
  strand_policy <- match.arg(strand_policy)
  dt <- .query_bwg_internal(
    object = object,
    chrom = chrom,
    start = start,
    end = end,
    samples = samples,
    strand = strand,
    strand_policy = strand_policy,
    verbose = verbose,
    progress = progress,
    keep_empty_samples = keep_empty_samples,
    tabix_empty_fallback = tabix_empty_fallback
  )
  if (as == "data.table") return(dt)
  if (as == "GRanges") {
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(dt$start, dt$end),
      sample_id = dt$sample_id,
      score = dt$value,
      strand = dt$strand
    ))
  }
  sample_tbl <- object$samples
  if (!is.null(samples)) {
    sample_tbl <- sample_tbl[sample_tbl[["sample_id"]] %in% as.character(samples)]
  }
  BwgTrack(samples = sample_tbl, data = dt, meta = modifyList(object$meta, list(mode = "memory")), validation = make_empty_validation())
}

