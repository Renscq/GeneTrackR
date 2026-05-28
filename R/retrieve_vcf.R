# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.19
# Function: Retrieve variants from VariantTrack objects
# Input: VariantTrack object and filters
# Output: data.table or VariantTrack object

retrieve_vcf <- function(object,
                         pattern = NULL,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         variant_id = NULL,
                         variant_type = NULL,
                         ignore_case = TRUE,
                         fixed = FALSE,
                         as = c("data.table", "VariantTrack")) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  as <- match.arg(as)
  dt <- data.table::copy(object$data)
  if (!is.null(chrom)) dt <- dt[dt[["chrom"]] %in% as.character(chrom)]
  if (!is.null(start) && !is.null(end)) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    dt <- dt[dt[["pos"]] >= s & dt[["pos"]] <= e]
  }
  if (!is.null(variant_id)) dt <- dt[dt[["variant_id"]] %in% as.character(variant_id)]
  if (!is.null(variant_type)) dt <- dt[dt[["variant_type"]] %in% as.character(variant_type)]
  dt <- match_pattern_internal(dt, pattern = pattern,
                               fields = c("variant_id", "ref", "alt", "info", "variant_type"),
                               ignore_case = ignore_case, fixed = fixed)
  if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "pos", "variant_id"))
  if (as == "data.table") return(dt[])
  VariantTrack(dt, meta = object$meta)
}

