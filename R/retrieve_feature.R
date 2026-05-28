# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.21
# Function: Retrieve Feature/GenePred-compatible annotation sub-objects or tables
# Input: Feature-compatible annotation object
# Output: Retrieved Feature/GenePred-compatible sub-object or data.table

#' Retrieve annotation features or sub-objects
#'
#' @description
#' `retrieve_feature()` extracts a true sub-object from a Feature/GenePred-compatible
#' annotation object. By default it returns a Feature/GenePred-compatible object.
#' Use `as = "data.table"` only when a plain table is required.
#'
#' @param object A Feature-compatible annotation object.
#' @param pattern Optional pattern used for ID/name matching.
#' @param level Output table level when `as = "data.table"`. One of feature, gene,
#' transcript, or exon.
#' @param chrom Optional chromosome.
#' @param start Optional start coordinate.
#' @param end Optional end coordinate.
#' @param mode Region matching mode. `overlap` keeps overlapping transcripts/features,
#' `within` keeps fully contained records, and `trim` clips gene-model records to
#' the query interval.
#' @param gene_id Optional gene ID filter. Exact matching is used.
#' @param transcript_id Optional transcript ID filter. Exact matching is used.
#' @param type Optional feature type filter for feature-table retrieval.
#' @param ignore_case Logical. Whether pattern matching ignores case.
#' @param fixed Logical. Whether pattern is fixed text.
#' @param as Output type. Default `Feature` returns a sub-object. Use `data.table`
#' for a plain table.
#'
#' @return A Feature/GenePred-compatible object or a data.table.
#' @export
retrieve_feature <- function(object,
                             pattern = NULL,
                             level = c("feature", "gene", "transcript", "exon"),
                             chrom = NULL,
                             start = NULL,
                             end = NULL,
                             mode = c("overlap", "within", "trim"),
                             gene_id = NULL,
                             transcript_id = NULL,
                             type = NULL,
                             ignore_case = TRUE,
                             fixed = FALSE,
                             as = c("Feature", "data.table")) {
  stop_if_not(
    inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"),
    "`object` must be a Feature-compatible annotation object."
  )

  level <- match.arg(level)
  mode <- match.arg(mode)
  as <- match.arg(as)

  has_region <- !is.null(chrom) && !is.null(start) && !is.null(end)

  normalize_id <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    unique(x)
  }

  gene_id <- normalize_id(gene_id)
  transcript_id <- normalize_id(transcript_id)
  type <- normalize_id(type)

  pattern_filter <- function(dt, fields) {
    dt <- data.table::as.data.table(dt)
    if (is.null(pattern)) {
      return(dt)
    }
    stop_if_not(
      length(pattern) == 1L && !is.na(pattern) && nzchar(pattern),
      "`pattern` must be a non-empty character string."
    )
    fields <- intersect(fields, names(dt))
    if (length(fields) == 0L || nrow(dt) == 0L) {
      return(dt[0])
    }
    hit <- rep(FALSE, nrow(dt))
    for (field in fields) {
      hit <- hit | grepl(
        pattern = pattern,
        x = as.character(dt[[field]]),
        ignore.case = ignore_case,
        fixed = fixed
      )
    }
    dt[hit]
  }

  exact_filter <- function(dt, field, values) {
    dt <- data.table::as.data.table(dt)
    if (length(values) == 0L) {
      return(dt)
    }
    if (!field %in% names(dt)) {
      return(dt[0])
    }
    dt[as.character(dt[[field]]) %in% values]
  }

  region_filter <- function(dt, start_col, end_col) {
    dt <- data.table::as.data.table(dt)
    if (!is.null(chrom) && "chrom" %in% names(dt)) {
      dt <- dt[as.character(dt[["chrom"]]) %in% as.character(chrom)]
    }
    if (has_region) {
      stop_if_not(
        all(c(start_col, end_col) %in% names(dt)),
        "Region filtering requires valid start and end columns."
      )
      query_start <- as.integer(start)[1L]
      query_end <- as.integer(end)[1L]
      stop_if_not(
        !is.na(query_start) && !is.na(query_end) && query_start <= query_end,
        "`start` and `end` must be valid integers with start <= end."
      )
      if (identical(mode, "within")) {
        dt <- dt[as.integer(dt[[start_col]]) >= query_start & as.integer(dt[[end_col]]) <= query_end]
      } else {
        dt <- dt[as.integer(dt[[start_col]]) <= query_end & as.integer(dt[[end_col]]) >= query_start]
      }
    }
    dt
  }

  make_gene_model_subobject <- function(tx) {
    tx <- data.table::copy(data.table::as.data.table(tx))
    ex <- data.table::copy(object$exons)

    if (nrow(tx) == 0L) {
      ex <- ex[0]
    } else {
      ex <- ex[as.character(ex[["transcript_id"]]) %in% unique(as.character(tx[["transcript_id"]]))]
    }

    if (has_region && identical(mode, "trim") && nrow(tx) > 0L) {
      query_start <- as.integer(start)[1L]
      query_end <- as.integer(end)[1L]

      tx[, "tx_start" := as.integer(pmax(as.integer(tx[["tx_start"]]), query_start))]
      tx[, "tx_end" := as.integer(pmin(as.integer(tx[["tx_end"]]), query_end))]

      if ("cds_start" %in% names(tx)) {
        tx[, "cds_start" := as.integer(pmax(as.integer(tx[["cds_start"]]), query_start))]
      }
      if ("cds_end" %in% names(tx)) {
        tx[, "cds_end" := as.integer(pmin(as.integer(tx[["cds_end"]]), query_end))]
      }

      ex[, "exon_start" := as.integer(pmax(as.integer(ex[["exon_start"]]), query_start))]
      ex[, "exon_end" := as.integer(pmin(as.integer(ex[["exon_end"]]), query_end))]
      ex <- ex[as.integer(ex[["exon_start"]]) <= as.integer(ex[["exon_end"]])]
      tx <- tx[as.character(tx[["transcript_id"]]) %in% unique(as.character(ex[["transcript_id"]]))]

      if (nrow(tx) > 0L) {
        exon_counts <- ex[, .(exon_count_new = as.integer(.N)), by = "transcript_id"]
        tx <- merge(
          tx[, setdiff(names(tx), "exon_count"), with = FALSE],
          exon_counts,
          by = "transcript_id",
          all.x = TRUE
        )
        data.table::setnames(tx, "exon_count_new", "exon_count")
      }
    }

    if (nrow(tx) > 0L) {
      data.table::setorderv(
        tx,
        intersect(c("chrom", "tx_start", "tx_end", "gene_id", "transcript_id"), names(tx))
      )
    }
    if (nrow(ex) > 0L) {
      data.table::setorderv(
        ex,
        intersect(c("chrom", "exon_start", "exon_end", "gene_id", "transcript_id", "exon_number"), names(ex))
      )
    }

    genes <- build_gene_table(tx)
    data <- genepred_to_feature_table(tx, ex, genes)

    out <- Feature(
      data = data,
      genes = genes,
      transcripts = tx,
      exons = ex,
      meta = modifyList(object$meta %||% list(), list(format = "retrieved", coordinate_internal = "1-based closed")),
      validation = make_empty_validation()
    )
    class(out) <- unique(c(class(out), "GenePred"))
    out
  }

  if (is_gene_model_feature(object)) {
    tx <- data.table::copy(object$transcripts)

    tx <- region_filter(tx, "tx_start", "tx_end")
    tx <- exact_filter(tx, "gene_id", gene_id)
    tx <- exact_filter(tx, "transcript_id", transcript_id)
    tx <- pattern_filter(tx, c("transcript_id", "gene_id", "gene_type"))

    if (identical(as, "Feature")) {
      return(make_gene_model_subobject(tx))
    }

    if (identical(level, "transcript")) {
      return(tx[])
    }

    if (identical(level, "gene")) {
      genes <- build_gene_table(tx)
      return(genes[])
    }

    if (identical(level, "exon")) {
      ex <- data.table::copy(object$exons)
      if (nrow(tx) == 0L) {
        ex <- ex[0]
      } else {
        ex <- ex[as.character(ex[["transcript_id"]]) %in% unique(as.character(tx[["transcript_id"]]))]
      }
      ex <- exact_filter(ex, "gene_id", gene_id)
      ex <- exact_filter(ex, "transcript_id", transcript_id)
      if (!is.null(chrom) || has_region) {
        ex <- region_filter(ex, "exon_start", "exon_end")
      }
      ex <- pattern_filter(ex, c("transcript_id", "gene_id"))
      return(ex[])
    }

    sub_obj <- make_gene_model_subobject(tx)
    dt <- data.table::copy(sub_obj$data)
    dt <- exact_filter(dt, "type", type)
    dt <- pattern_filter(dt, c("feature_id", "name", "gene_id", "transcript_id", "parent_id", "type"))
    return(dt[])
  }

  dt <- data.table::copy(as_feature_table(object))
  dt <- region_filter(dt, "start", "end")
  dt <- exact_filter(dt, "gene_id", gene_id)
  dt <- exact_filter(dt, "transcript_id", transcript_id)
  dt <- exact_filter(dt, "type", type)
  dt <- pattern_filter(dt, c("feature_id", "name", "gene_id", "transcript_id", "parent_id", "type"))

  if (identical(as, "data.table")) {
    return(dt[])
  }

  Feature(
    data = dt,
    meta = modifyList(object$meta %||% list(), list(format = "retrieved", coordinate_internal = "1-based closed")),
    validation = make_empty_validation()
  )
}
