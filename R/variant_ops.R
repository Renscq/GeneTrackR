# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.2
# Function: Operations for genome-level VariantTrack objects
# Input: VariantTrack objects
# Output: Summaries, sliced objects, merged objects, and VCF-like files

#' Summarize a VariantTrack object
#'
#' @param object A VariantTrack object.
#' @param chrom Optional chromosome.
#' @param start Optional start coordinate.
#' @param end Optional end coordinate.
#' @return A data.table summary by chromosome and variant type.
#' @export
summary_variant_track <- function(object, chrom = NULL, start = NULL, end = NULL) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  dt <- slice_variant_track(object, chrom = chrom, start = start, end = end)$data
  if (nrow(dt) == 0L) return(data.table::data.table())
  dt[, .(n_variants = as.integer(.N)), by = .(chrom, variant_type)]
}

#' Slice a VariantTrack object
#'
#' @param object A VariantTrack object.
#' @param chrom Optional chromosome.
#' @param start Optional start coordinate.
#' @param end Optional end coordinate.
#' @return A VariantTrack object.
#' @export
slice_variant_track <- function(object, chrom = NULL, start = NULL, end = NULL) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  dt <- data.table::copy(object$data)
  if (!is.null(chrom)) dt <- dt[dt[["chrom"]] == as.character(chrom)[1L]]
  if (!is.null(start) && !is.null(end)) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    dt <- dt[dt[["pos"]] >= s & dt[["pos"]] <= e]
  }
  VariantTrack(dt, meta = object$meta)
}

#' Merge VariantTrack objects
#'
#' @param ... VariantTrack objects.
#' @param source_names Optional source names to store in the `track_source` column.
#' @return A VariantTrack object.
#' @export
merge_variant_track <- function(..., source_names = NULL) {
  tracks <- list(...)
  stop_if_not(length(tracks) > 0L, "At least one VariantTrack object is required.")
  stop_if_not(all(vapply(tracks, inherits, logical(1L), "VariantTrack")), "All inputs must be VariantTrack objects.")
  if (is.null(source_names)) source_names <- paste0("track", seq_along(tracks))
  out <- lapply(seq_along(tracks), function(i) {
    dt <- data.table::copy(tracks[[i]]$data)
    dt[, "track_source" := as.character(source_names[i])]
    dt
  })
  VariantTrack(data.table::rbindlist(out, fill = TRUE), meta = list(format = "merged", coordinate_internal = "1-based position"))
}

