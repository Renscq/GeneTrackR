# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.14
# Function: Unified summary API for genome-level VariantTrack objects
# Input: VariantTrack objects
# Output: Variant summary tables

#' Summarize a VariantTrack object
#'
#' @description
#' `summary_vcf()` is the unified summary API for VCF-derived VariantTrack
#' objects.
#'
#' @param object A VariantTrack object.
#' @param chrom Optional chromosome filter.
#' @param start Optional start coordinate.
#' @param end Optional end coordinate.
#' @param by Grouping columns. Default `c("chrom", "variant_type")`.
#'
#' @return A data.table summary.
#' @export
summary_vcf <- function(object, chrom = NULL, start = NULL, end = NULL, by = c("chrom", "variant_type")) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  dt <- retrieve_vcf(object, chrom = chrom, start = start, end = end)
  if (nrow(dt) == 0L) return(data.table::data.table())
  by <- intersect(as.character(by), names(dt))
  if (length(by) == 0L) by <- "chrom"
  out <- dt[, .(n_variants = as.integer(.N)), by = by]
  data.table::setorderv(out, by)
  out[]
}
