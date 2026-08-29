# Author: Rensc
# Date: 2026-08-29
# Version: dev003
# Function: Retrieve signal records from BwgTrack objects
# Input: BwgTrack object and genomic region
# Output: data.table, BwgTrack, or GRanges

#' Retrieve signal records from a BwgTrack object
#'
#' @description
#' `retrieve_bwg()` is the unified signal retrieval API in GeneTrackR. It
#' retrieves signal intervals from bedGraph, wig, or bigWig tracks by genomic
#' region.
#'
#' @param object A BwgTrack object returned by `read_bwg()`.
#' @param chrom Chromosome name. Required for direct region queries.
#' @param start Region start in 1-based closed coordinates. Required for direct region queries.
#' @param end Region end in 1-based closed coordinates. Required for direct region queries.
#' @param annotation Optional GenePred/Feature annotation object used for gene/transcript-aware retrieval.
#' @param gene_id Optional gene ID. When supplied, `annotation` is used to resolve the gene range.
#' @param transcript_id Optional transcript ID. When supplied, `annotation` is used to resolve the transcript range.
#' @param upstream Upstream flanking length in bp for gene/transcript queries.
#' @param downstream Downstream flanking length in bp for gene/transcript queries.
#' @param strand_aware Logical. Whether upstream/downstream should follow gene/transcript strand direction.
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
#' @return A data.table, BwgTrack object, or GRanges object. BwgTrack subsets
#' retain sequence metadata for the selected samples and chromosome when
#' available.
#'
#' @examples
#' \dontrun{
#' rnaseq_files <- system.file(
#'   "extdata",
#'   c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
#'   package = "GeneTrackR"
#' )
#' rnaseq <- read_bwg(
#'   rnaseq_files,
#'   format = "bedgraph",
#'   sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#'
#' dt <- retrieve_bwg(
#'   rnaseq,
#'   chrom = "chr1",
#'   start = 12339001,
#'   end = 12352000,
#'   strand = "+"
#' )
#'
#' anno <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt",
#'   verbose = FALSE
#' )
#' gene_signal <- retrieve_bwg(
#'   rnaseq,
#'   annotation = anno,
#'   gene_id = "GeneA",
#'   upstream = 500,
#'   downstream = 500,
#'   strand = "auto"
#' )
#' }
#'
#' @export
retrieve_bwg <- function(object,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         annotation = NULL,
                         gene_id = NULL,
                         transcript_id = NULL,
                         upstream = 0L,
                         downstream = 0L,
                         strand_aware = TRUE,
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

  region <- resolve_retrieve_gene_region(
    annotation = annotation,
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    start = start,
    end = end,
    upstream = upstream,
    downstream = downstream,
    strand_aware = strand_aware
  )

  chrom <- region$chrom
  start <- region$start
  end <- region$end
  if (identical(strand, "auto") && !is.null(region$strand) && !is.na(region$strand) && region$strand %in% c("+", "-")) {
    strand <- region$strand
  }

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
  seqinfo <- subset_bwg_seqinfo(
    object,
    sample_ids = sample_tbl[["sample_id"]],
    chrom = chrom
  )
  BwgTrack(
    samples = sample_tbl,
    data = dt,
    meta = modifyList(object$meta, list(mode = "memory")),
    validation = make_empty_validation(),
    seqinfo = seqinfo
  )
}


resolve_retrieve_gene_region <- function(annotation = NULL,
                                         gene_id = NULL,
                                         transcript_id = NULL,
                                         chrom = NULL,
                                         start = NULL,
                                         end = NULL,
                                         upstream = 0L,
                                         downstream = 0L,
                                         strand_aware = TRUE) {
  has_gene <- !is.null(gene_id)
  has_tx <- !is.null(transcript_id)
  has_direct_region <- !is.null(chrom) || !is.null(start) || !is.null(end)

  if (has_gene || has_tx) {
    stop_if_not(!has_direct_region, "Do not mix `gene_id`/`transcript_id` with direct `chrom`/`start`/`end` region arguments.")
    return(resolve_haplotype_gene_region(
      annotation = annotation,
      gene_id = gene_id,
      transcript_id = transcript_id,
      upstream = upstream,
      downstream = downstream,
      strand_aware = strand_aware
    ))
  }

  stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end), "Specify either `gene_id`/`transcript_id` with `annotation`, or direct `chrom`, `start`, and `end`.")
  check_region(chrom, start, end)

  list(
    locator = "region",
    id = paste0(as.character(chrom)[1L], ":", as.integer(start)[1L], "-", as.integer(end)[1L]),
    chrom = as.character(chrom)[1L],
    start = as.integer(start)[1L],
    end = as.integer(end)[1L],
    core_start = as.integer(start)[1L],
    core_end = as.integer(end)[1L],
    upstream = 0L,
    downstream = 0L,
    strand = NA_character_,
    strand_aware = FALSE
  )
}
