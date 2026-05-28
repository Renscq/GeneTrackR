#' Validate a standardized Feature annotation object
#'
#' @param object A Feature-compatible annotation object.
#' @return A validation list with invalid records, summaries, and warnings.
#' @export
validate_feature <- function(object) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"), "`object` must be a Feature-compatible annotation object.")
  if (is_gene_model_feature(object)) {
    return(validate_genepred(object))
  }
  dt <- data.table::copy(as_feature_table(object))
  if (!"row_id" %in% names(dt)) dt[, row_id := .I]
  invalid <- list()
  add_invalid <- function(rows, reason) {
    rows <- unique(rows[!is.na(rows)])
    if (length(rows) == 0L) return(NULL)
    data.table::data.table(row_id = rows, reason = reason)
  }
  invalid[[length(invalid) + 1L]] <- add_invalid(dt[is.na(chrom) | chrom == "", row_id], "empty chromosome")
  invalid[[length(invalid) + 1L]] <- add_invalid(dt[is.na(start) | is.na(end) | start > end, row_id], "invalid coordinates")
  invalid[[length(invalid) + 1L]] <- add_invalid(dt[!strand %in% c("+", "-", "*"), row_id], "invalid strand")
  invalid_records <- data.table::rbindlist(invalid, fill = TRUE)
  if (nrow(invalid_records) == 0L) invalid_records <- data.table::data.table(row_id = integer(), reason = character())
  list(
    invalid_records = invalid_records,
    invalid_summary = invalid_records[, .N, by = reason][order(-N)],
    warnings = character()
  )
}

# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.0
# Function: Validate standardized Feature annotation objects and gene-model tables
# Input: GenePred object or transcript/exon tables
# Output: Validation reports

#' Validate a GenePred object
#'
#' @param object A GenePred object.
#' @return A validation list with invalid records, summaries, and warnings.
#' @export
validate_genepred <- function(object) {
  stop_if_not(is_gene_model_feature(object), "`object` must contain gene-model tables.")
  tx <- data.table::copy(object$transcripts)
  ex <- data.table::copy(object$exons)
  if (!"row_id" %in% names(tx)) {
    tx[, row_id := .I]
    ex <- merge(ex, tx[, .(row_id, transcript_id)], by = "transcript_id", all.x = TRUE)
  }
  validate_genepred_tables(tx, ex)
}

validate_genepred_tables <- function(transcripts, exons) {
  tx <- data.table::copy(transcripts)
  ex <- data.table::copy(exons)
  invalid <- list()

  add_invalid <- function(rows, reason) {
    rows <- unique(rows[!is.na(rows)])
    if (length(rows) == 0L) {
      return(NULL)
    }
    data.table::data.table(row_id = rows, reason = reason)
  }

  invalid[[length(invalid) + 1L]] <- add_invalid(tx[is.na(chrom) | chrom == "", row_id], "empty chromosome")
  invalid[[length(invalid) + 1L]] <- add_invalid(tx[is.na(transcript_id) | transcript_id == "", row_id], "empty transcript_id")
  invalid[[length(invalid) + 1L]] <- add_invalid(tx[is.na(gene_id) | gene_id == "", row_id], "empty gene_id")
  invalid[[length(invalid) + 1L]] <- add_invalid(tx[!strand %in% c("+", "-"), row_id], "invalid strand")
  invalid[[length(invalid) + 1L]] <- add_invalid(tx[is.na(tx_start) | is.na(tx_end) | tx_start > tx_end, row_id], "invalid transcript coordinates")

  coding_bad <- tx[gene_type == "coding" & (is.na(cds_start) | is.na(cds_end) | cds_start > cds_end | cds_start < tx_start | cds_end > tx_end), row_id]
  invalid[[length(invalid) + 1L]] <- add_invalid(coding_bad, "invalid CDS coordinates")

  exon_n <- ex[, .N, by = row_id]
  tx_n <- merge(tx[, .(row_id, exon_count)], exon_n, by = "row_id", all.x = TRUE)
  tx_n[is.na(N), N := 0L]
  invalid[[length(invalid) + 1L]] <- add_invalid(tx_n[exon_count != N, row_id], "exonCount does not match exonStarts/exonEnds")

  invalid[[length(invalid) + 1L]] <- add_invalid(ex[is.na(exon_start) | is.na(exon_end) | exon_start > exon_end, unique(row_id)], "invalid exon coordinates")

  ex_tx <- merge(
    ex,
    tx[, .(row_id, tx_start, tx_end)],
    by = "row_id",
    all.x = TRUE
  )
  invalid[[length(invalid) + 1L]] <- add_invalid(ex_tx[exon_start < tx_start | exon_end > tx_end, unique(row_id)], "exon outside transcript region")

  overlap_bad <- ex[order(row_id, exon_start, exon_end), {
    bad <- any(exon_start[-1L] <= exon_end[-.N])
    .(bad = bad)
  }, by = row_id][bad == TRUE, row_id]
  invalid[[length(invalid) + 1L]] <- add_invalid(overlap_bad, "overlapping exons")

  invalid_records <- data.table::rbindlist(invalid, fill = TRUE)
  if (nrow(invalid_records) == 0L) {
    invalid_records <- data.table::data.table(row_id = integer(), reason = character())
  }

  invalid_summary <- invalid_records[, .N, by = reason][order(-N)]
  warnings <- character()

  multi_gene <- tx[, data.table::uniqueN(gene_id), by = transcript_id][V1 > 1L]
  if (nrow(multi_gene) > 0L) {
    warnings <- c(warnings, "Some transcript IDs are assigned to multiple gene IDs.")
  }

  gene_multi_chrom <- tx[, data.table::uniqueN(chrom), by = gene_id][V1 > 1L]
  if (nrow(gene_multi_chrom) > 0L) {
    warnings <- c(warnings, "Some gene IDs are assigned to multiple chromosomes.")
  }

  gene_multi_strand <- tx[, data.table::uniqueN(strand), by = gene_id][V1 > 1L]
  if (nrow(gene_multi_strand) > 0L) {
    warnings <- c(warnings, "Some gene IDs are assigned to multiple strands.")
  }

  list(
    invalid_records = invalid_records,
    invalid_summary = invalid_summary,
    warnings = warnings
  )
}
