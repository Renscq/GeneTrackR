# Author: Rensc
# Date: 2026-08-31
# Version: dev002
# Function: Unified summary API for Feature-compatible annotation objects
# Input: Feature, FeatureTrack, or GenePred-compatible annotation objects
# Output: Summary tables

#' Summarize a Feature-compatible annotation object
#'
#' @description
#' `summary_feature()` is the unified summary API for GenePred, GTF, GFF, and
#' BED-derived annotation objects.
#'
#' @param object A Feature-compatible annotation object.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start.
#' @param end Optional region end.
#' @param level Summary level. Use `feature`, `gene`, `transcript`, or `exon`.
#' @param by Grouping columns used when `level = "feature"`.
#'
#' @return A data.table summary.
#' @examples
#' gp_file <- system.file(
#'   "extdata", "gtr_demo.genePredExt", package = "GeneTrackR"
#' )
#' gp <- read_genepred(
#'   gp_file, format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' summary_feature(gp, level = "gene")
#' @export
summary_feature <- function(object,
                            chrom = NULL,
                            start = NULL,
                            end = NULL,
                            level = c("feature", "gene", "transcript", "exon"),
                            by = c("chrom", "type")) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"),
              "`object` must be a Feature-compatible annotation object.")
  level <- match.arg(level)

  if (level == "feature") {
    dt <- retrieve_feature(object, chrom = chrom, start = start, end = end, level = "feature", as = "data.table")
    if (nrow(dt) == 0L) return(data.table::data.table())
    by <- intersect(as.character(by), names(dt))
    if (length(by) == 0L) by <- "chrom"
    out <- dt[, .(
      n_features = as.integer(.N),
      median_feature_length = as.numeric(stats::median(as.numeric(end - start + 1L), na.rm = TRUE))
    ), by = by]
    data.table::setorderv(out, by)
    return(out[])
  }

  gp <- as_genepred(object)

  if (level == "gene") {
    dt <- retrieve_feature(gp, chrom = chrom, start = start, end = end, level = "gene", as = "data.table")
    if (nrow(dt) == 0L) {
      return(data.table::data.table(
        chrom = character(),
        strand = character(),
        n_genes = integer(),
        n_coding_genes = integer(),
        n_noncoding_genes = integer(),
        median_gene_length = numeric()
      ))
    }
    return(dt[, .(
      n_genes = as.integer(.N),
      n_coding_genes = as.integer(sum(gene_type == "coding", na.rm = TRUE)),
      n_noncoding_genes = as.integer(sum(gene_type != "coding" | is.na(gene_type), na.rm = TRUE)),
      median_gene_length = as.numeric(stats::median(as.numeric(gene_end - gene_start + 1L), na.rm = TRUE))
    ), by = .(chrom, strand)][order(chrom, strand)])
  }

  if (level == "transcript") {
    dt <- retrieve_feature(gp, chrom = chrom, start = start, end = end, level = "transcript", as = "data.table")
    if (nrow(dt) == 0L) {
      return(data.table::data.table(
        chrom = character(),
        strand = character(),
        n_transcripts = integer(),
        n_coding_transcripts = integer(),
        n_noncoding_transcripts = integer(),
        median_transcript_length = numeric(),
        median_exon_count = numeric()
      ))
    }
    return(dt[, .(
      n_transcripts = as.integer(.N),
      n_coding_transcripts = as.integer(sum(gene_type == "coding", na.rm = TRUE)),
      n_noncoding_transcripts = as.integer(sum(gene_type != "coding" | is.na(gene_type), na.rm = TRUE)),
      median_transcript_length = as.numeric(stats::median(as.numeric(tx_end - tx_start + 1L), na.rm = TRUE)),
      median_exon_count = as.numeric(stats::median(as.numeric(exon_count), na.rm = TRUE))
    ), by = .(chrom, strand)][order(chrom, strand)])
  }

  dt <- retrieve_feature(gp, chrom = chrom, start = start, end = end, level = "exon", as = "data.table")
  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      chrom = character(),
      strand = character(),
      n_exons = integer(),
      median_exon_length = numeric()
    ))
  }
  dt[, .(
    n_exons = as.integer(.N),
    median_exon_length = as.numeric(stats::median(as.numeric(exon_end - exon_start + 1L), na.rm = TRUE))
  ), by = .(chrom, strand)][order(chrom, strand)]
}
