# Author: Rensc
# Date: 2026-05-28
# Version: dev001
# Function: Validate unified Feature annotation objects
# Input: Feature-compatible annotation objects
# Output: Validation reports

make_validation_report <- function(invalid_records, warnings = character()) {
  invalid_records <- data.table::as.data.table(invalid_records)
  if (nrow(invalid_records) == 0L) {
    if (!"row_id" %in% names(invalid_records)) invalid_records[, "row_id" := integer()]
    if (!"reason" %in% names(invalid_records)) invalid_records[, "reason" := character()]
  }
  invalid_summary <- invalid_records[, .N, by = "reason"]
  if (nrow(invalid_summary) > 0L) data.table::setorder(invalid_summary, -N, reason)
  list(
    invalid_records = invalid_records[],
    invalid_summary = invalid_summary[],
    warnings = as.character(warnings)
  )
}

append_invalid <- function(container, rows, reason) {
  rows <- unique(rows[!is.na(rows)])
  if (length(rows) == 0L) return(container)
  container[[length(container) + 1L]] <- data.table::data.table(row_id = rows, reason = reason)
  container
}

validate_feature_table <- function(dt) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  if (!"row_id" %in% names(dt)) dt[, "row_id" := .I]
  invalid <- list()
  warnings <- character()

  required <- c("chrom", "start", "end")
  missing_required <- setdiff(required, names(dt))
  if (length(missing_required) > 0L) {
    invalid[[1L]] <- data.table::data.table(
      row_id = integer(),
      reason = paste0("missing required columns: ", paste(missing_required, collapse = ", "))
    )
    return(make_validation_report(invalid[[1L]], warnings))
  }

  invalid <- append_invalid(invalid, dt[is.na(chrom) | chrom == "", row_id], "empty chromosome")
  invalid <- append_invalid(invalid, dt[is.na(start) | is.na(end) | start > end, row_id], "invalid coordinates")

  if ("strand" %in% names(dt)) {
    invalid <- append_invalid(invalid, dt[is.na(strand) | !strand %in% c("+", "-", "*", "."), row_id], "invalid strand")
  }
  if ("feature_id" %in% names(dt)) {
    duplicated_id <- dt[!is.na(feature_id) & feature_id != "" & duplicated(feature_id), row_id]
    invalid <- append_invalid(invalid, duplicated_id, "duplicated feature_id")
  }
  if ("type" %in% names(dt) && any(is.na(dt$type) | dt$type == "")) {
    warnings <- c(warnings, "Some feature records have missing feature type.")
  }

  invalid_records <- data.table::rbindlist(invalid, fill = TRUE)
  make_validation_report(invalid_records, warnings)
}

validate_gene_model_tables <- function(transcripts, exons) {
  tx <- data.table::copy(data.table::as.data.table(transcripts))
  ex <- data.table::copy(data.table::as.data.table(exons))
  invalid <- list()
  warnings <- character()

  if (!"row_id" %in% names(tx)) {
    tx[, "row_id" := .I]
  }
  if (!"row_id" %in% names(ex)) {
    if ("transcript_id" %in% names(ex) && "transcript_id" %in% names(tx)) {
      ex <- merge(ex, tx[, .(row_id, transcript_id)], by = "transcript_id", all.x = TRUE)
    } else {
      ex[, "row_id" := NA_integer_]
    }
  }

  required_tx <- c("chrom", "transcript_id", "gene_id", "strand", "tx_start", "tx_end")
  missing_tx <- setdiff(required_tx, names(tx))
  if (length(missing_tx) > 0L) {
    invalid[[length(invalid) + 1L]] <- data.table::data.table(
      row_id = integer(),
      reason = paste0("missing transcript columns: ", paste(missing_tx, collapse = ", "))
    )
    return(make_validation_report(data.table::rbindlist(invalid, fill = TRUE), warnings))
  }

  invalid <- append_invalid(invalid, tx[is.na(chrom) | chrom == "", row_id], "empty chromosome")
  invalid <- append_invalid(invalid, tx[is.na(transcript_id) | transcript_id == "", row_id], "empty transcript_id")
  invalid <- append_invalid(invalid, tx[is.na(gene_id) | gene_id == "", row_id], "empty gene_id")
  invalid <- append_invalid(invalid, tx[is.na(strand) | !strand %in% c("+", "-"), row_id], "invalid strand")
  invalid <- append_invalid(invalid, tx[is.na(tx_start) | is.na(tx_end) | tx_start > tx_end, row_id], "invalid transcript coordinates")

  if (all(c("gene_type", "cds_start", "cds_end") %in% names(tx))) {
    coding_bad <- tx[gene_type == "coding" & (is.na(cds_start) | is.na(cds_end) | cds_start > cds_end | cds_start < tx_start | cds_end > tx_end), row_id]
    invalid <- append_invalid(invalid, coding_bad, "invalid CDS coordinates")
  }

  if (nrow(ex) > 0L) {
    required_ex <- c("exon_start", "exon_end")
    missing_ex <- setdiff(required_ex, names(ex))
    if (length(missing_ex) > 0L) {
      invalid[[length(invalid) + 1L]] <- data.table::data.table(
        row_id = integer(),
        reason = paste0("missing exon columns: ", paste(missing_ex, collapse = ", "))
      )
    } else {
      if ("exon_count" %in% names(tx)) {
        exon_n <- ex[, .N, by = "row_id"]
        tx_n <- merge(tx[, .(row_id, exon_count)], exon_n, by = "row_id", all.x = TRUE)
        tx_n[is.na(N), "N" := 0L]
        invalid <- append_invalid(invalid, tx_n[exon_count != N, row_id], "exonCount does not match exon records")
      }
      invalid <- append_invalid(invalid, ex[is.na(exon_start) | is.na(exon_end) | exon_start > exon_end, unique(row_id)], "invalid exon coordinates")

      ex_tx <- merge(ex, tx[, .(row_id, tx_start, tx_end)], by = "row_id", all.x = TRUE)
      invalid <- append_invalid(invalid, ex_tx[exon_start < tx_start | exon_end > tx_end, unique(row_id)], "exon outside transcript region")

      overlap_bad <- ex[order(row_id, exon_start, exon_end), {
        bad <- if (.N <= 1L) FALSE else any(exon_start[-1L] <= exon_end[-.N])
        .(bad = bad)
      }, by = "row_id"][bad == TRUE, row_id]
      invalid <- append_invalid(invalid, overlap_bad, "overlapping exons")
    }
  }

  multi_gene <- tx[, data.table::uniqueN(gene_id), by = "transcript_id"][V1 > 1L]
  if (nrow(multi_gene) > 0L) warnings <- c(warnings, "Some transcript IDs are assigned to multiple gene IDs.")

  gene_multi_chrom <- tx[, data.table::uniqueN(chrom), by = "gene_id"][V1 > 1L]
  if (nrow(gene_multi_chrom) > 0L) warnings <- c(warnings, "Some gene IDs are assigned to multiple chromosomes.")

  gene_multi_strand <- tx[, data.table::uniqueN(strand), by = "gene_id"][V1 > 1L]
  if (nrow(gene_multi_strand) > 0L) warnings <- c(warnings, "Some gene IDs are assigned to multiple strands.")

  invalid_records <- data.table::rbindlist(invalid, fill = TRUE)
  make_validation_report(invalid_records, warnings)
}

merge_validation_reports <- function(...) {
  reports <- list(...)
  invalid <- data.table::rbindlist(lapply(reports, `[[`, "invalid_records"), fill = TRUE)
  warnings <- unique(unlist(lapply(reports, `[[`, "warnings"), use.names = FALSE))
  make_validation_report(invalid, warnings)
}

#' Validate a Feature annotation object
#'
#' @description
#' `validate_feature()` is the unified validation entry for annotation data. It
#' validates BED/GFF/GTF/GenePred-derived `Feature` objects. If gene-model tables
#' are present, transcript/exon consistency is checked in addition to the flat
#' standardized feature table.
#'
#' @param object A Feature-compatible annotation object.
#' @param check_gene_model Logical. Whether to check derived gene-model tables
#' when transcripts and exons are available.
#'
#' @return A validation list with `invalid_records`, `invalid_summary`, and
#' `warnings`.
#'
#' @export
validate_feature <- function(object, check_gene_model = TRUE) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"),
              "`object` must be a Feature-compatible annotation object.")
  feature_report <- validate_feature_table(as_feature_table(object))
  if (isTRUE(check_gene_model) && is_gene_model_feature(object)) {
    model_report <- validate_gene_model_tables(object$transcripts, object$exons)
    return(merge_validation_reports(feature_report, model_report))
  }
  feature_report
}

# Compatibility wrapper. Not exported in the unified API.
validate_genepred <- function(object) {
  validate_feature(object, check_gene_model = TRUE)
}

# Internal compatibility function used by read_genepred().
validate_genepred_tables <- function(transcripts, exons) {
  validate_gene_model_tables(transcripts, exons)
}
