# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.13
# Function: Validate VCF-like VariantTrack objects
# Input: VariantTrack objects
# Output: Validation reports

#' Validate a VariantTrack object
#'
#' @description
#' `validate_vcf()` validates genome-level variant records stored in a
#' `VariantTrack` object. It checks chromosome and position fields, allele
#' fields, duplicated variant IDs, and consistency between `pos`, `start`, and
#' `end` when those columns are available.
#'
#' @param object A VariantTrack object.
#'
#' @return A validation list with `invalid_records`, `invalid_summary`, and
#' `warnings`.
#'
#' @export
validate_vcf <- function(object) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  dt <- data.table::copy(data.table::as.data.table(object$data))
  if (!"row_id" %in% names(dt)) dt[, "row_id" := .I]
  invalid <- list()
  warnings <- character()

  required <- c("chrom", "pos", "variant_id", "ref", "alt")
  missing_required <- setdiff(required, names(dt))
  if (length(missing_required) > 0L) {
    invalid[[length(invalid) + 1L]] <- data.table::data.table(
      row_id = integer(),
      reason = paste0("missing required columns: ", paste(missing_required, collapse = ", "))
    )
  } else {
    invalid <- append_invalid(invalid, dt[is.na(chrom) | chrom == "", row_id], "empty chromosome")
    invalid <- append_invalid(invalid, dt[is.na(pos) | pos < 1L, row_id], "invalid variant position")
    invalid <- append_invalid(invalid, dt[is.na(variant_id) | variant_id == "", row_id], "empty variant_id")
    invalid <- append_invalid(invalid, dt[is.na(ref) | ref == "", row_id], "empty REF allele")
    invalid <- append_invalid(invalid, dt[is.na(alt) | alt == "", row_id], "empty ALT allele")

    dup <- dt[!is.na(variant_id) & variant_id != "" & duplicated(variant_id), row_id]
    invalid <- append_invalid(invalid, dup, "duplicated variant_id")

    if (all(c("start", "end") %in% names(dt))) {
      invalid <- append_invalid(invalid, dt[is.na(start) | is.na(end) | start > end, row_id], "invalid start/end coordinates")
      invalid <- append_invalid(invalid, dt[!is.na(pos) & (!is.na(start) & pos != start | !is.na(end) & pos != end), row_id], "pos does not match start/end")
    }
    if ("variant_type" %in% names(dt)) {
      unknown_type <- dt[is.na(variant_type) | variant_type == "", row_id]
      invalid <- append_invalid(invalid, unknown_type, "empty variant_type")
    }
  }

  invalid_records <- data.table::rbindlist(invalid, fill = TRUE)
  make_validation_report(invalid_records, warnings)
}
