# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.19
# Function: Merge VariantTrack objects
# Input: VariantTrack objects
# Output: Merged VariantTrack object

merge_vcf <- function(...,
                      source_names = NULL,
                      conflict = c("keep_all", "error", "rename", "keep_first"),
                      rename_prefix = NULL,
                      sort = TRUE) {
  tracks <- list(...)
  stop_if_not(length(tracks) > 0L, "At least one VariantTrack object is required.")
  stop_if_not(all(vapply(tracks, inherits, logical(1L), "VariantTrack")), "All inputs must be VariantTrack objects.")
  conflict <- match.arg(conflict)
  if (is.null(source_names)) source_names <- paste0("track", seq_along(tracks))
  stop_if_not(length(source_names) == length(tracks), "`source_names` must match the number of input objects.")
  if (conflict == "rename" && is.null(rename_prefix)) rename_prefix <- paste0(source_names, "_")
  if (!is.null(rename_prefix)) stop_if_not(length(rename_prefix) == length(tracks), "`rename_prefix` must match the number of input objects.")

  out <- lapply(seq_along(tracks), function(i) {
    dt <- data.table::copy(tracks[[i]]$data)
    if (conflict == "rename" && "variant_id" %in% names(dt)) {
      idx <- !is.na(dt[["variant_id"]]) & nzchar(dt[["variant_id"]])
      dt[idx, "variant_id" := paste0(rename_prefix[i], dt[["variant_id"]])]
    }
    dt[, "track_source" := as.character(source_names[i])]
    dt
  })
  dt <- data.table::rbindlist(out, fill = TRUE)

  if (conflict == "error" && "variant_id" %in% names(dt)) {
    dup <- dt[duplicated(variant_id) | duplicated(variant_id, fromLast = TRUE), unique(variant_id)]
    stop_if_not(length(dup) == 0L, "Duplicated variant IDs found. Use `conflict = 'rename'`, `keep_all`, or `keep_first`.")
  }
  if (conflict == "keep_first" && "variant_id" %in% names(dt)) {
    dt <- dt[!duplicated(variant_id)]
  }
  if (isTRUE(sort) && nrow(dt) > 0L) {
    data.table::setorderv(dt, intersect(c("chrom", "pos", "variant_id"), names(dt)))
  }
  VariantTrack(dt, meta = list(format = "merged", coordinate_internal = "1-based position"))
}

