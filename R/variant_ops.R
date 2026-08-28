# Author: Rensc
# Date: 2026-05-28
# Version: dev002
# Function: Unified summary API for genome-level VariantTrack objects
# Input: VariantTrack objects
# Output: Variant summary tables

#' Summarize a VariantTrack object
#'
#' @description
#' `summary_vcf()` is the unified summary API for VCF-derived VariantTrack
#' objects. When no genomic range is supplied, the full in-memory object is
#' summarized. For lazy tracks, an omitted range causes the source VCF to be
#' read before summarization.
#'
#' @param object A VariantTrack object.
#' @param chrom Optional chromosome filter. May be used without `start`/`end`.
#' @param start Optional 1-based start coordinate. Must be paired with `end`.
#' @param end Optional 1-based end coordinate. Must be paired with `start`.
#' @param by Grouping columns. Default `c("chrom", "variant_type")`.
#'
#' @return A data.table summary.
#' @examples
#' vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
#' summary_vcf(vcf)
#' summary_vcf(vcf, chrom = "chr1")
#' @export
summary_vcf <- function(object, chrom = NULL, start = NULL, end = NULL, by = c("chrom", "variant_type")) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  dt <- retrieve_vcf(object, chrom = chrom, start = start, end = end, verbose = FALSE)
  if (nrow(dt) == 0L) return(data.table::data.table())
  by <- intersect(as.character(by), names(dt))
  if (length(by) == 0L) by <- "chrom"
  out <- dt[, .(n_variants = as.integer(.N)), by = by]
  data.table::setorderv(out, by)
  out[]
}
