# Author: Rensc
# Date: 2026-08-31
# Version: dev004
# Function: BwgTrack S3 class constructors, schema normalization, and print methods
# Input: Signal sample metadata, optional in-memory signal table, and sequence metadata
# Output: Schema-v2 BwgTrack object

normalize_bwg_seqinfo <- function(seqinfo, sample_ids) {
  if (is.null(seqinfo)) {
    return(NULL)
  }

  x <- data.table::copy(data.table::as.data.table(seqinfo))
  stop_if_not(
    all(c("chrom", "length") %in% names(x)),
    "`seqinfo` must contain `chrom` and `length` columns."
  )

  sample_ids <- as.character(sample_ids)
  if (!"sample_id" %in% names(x)) {
    x <- data.table::rbindlist(lapply(sample_ids, function(sample_id) {
      y <- data.table::copy(x)
      y[, sample_id := sample_id]
      y
    }), fill = TRUE)
  }

  x[, `:=`(
    sample_id = as.character(sample_id),
    chrom = as.character(chrom)
  )]
  length_value <- suppressWarnings(as.numeric(x[["length"]]))

  stop_if_not(
    all(!is.na(x$sample_id) & nzchar(x$sample_id)),
    "`seqinfo` contains missing or empty sample IDs."
  )
  stop_if_not(
    all(!is.na(x$chrom) & nzchar(x$chrom)),
    "`seqinfo` contains missing or empty chromosome names."
  )

  unknown_samples <- setdiff(unique(x$sample_id), sample_ids)
  stop_if_not(
    length(unknown_samples) == 0L,
    paste0(
      "`seqinfo` contains unknown sample IDs: ",
      paste(unknown_samples, collapse = ", "),
      "."
    )
  )

  invalid_length <- !is.na(length_value) & (
    !is.finite(length_value) |
      length_value < 1 |
      length_value != floor(length_value)
  )
  stop_if_not(
    !any(invalid_length),
    "`seqinfo$length` must contain positive integer chromosome lengths or NA."
  )
  length_value[!is.na(length_value) & length_value > .Machine$integer.max] <- NA_real_
  x[, length := as.integer(length_value)]

  key <- paste(x$sample_id, x$chrom, sep = "\r")
  conflicting <- vapply(unique(key), function(k) {
    values <- unique(x$length[key == k])
    values <- values[!is.na(values)]
    length(values) > 1L
  }, logical(1))
  stop_if_not(
    !any(conflicting),
    "`seqinfo` contains conflicting chromosome lengths for a sample/chromosome pair."
  )

  sample_order <- match(x$sample_id, sample_ids)
  data.table::set(x, j = "sample_order__", value = sample_order)
  data.table::setorderv(x, c("sample_order__", "chrom"))
  x[, "sample_order__" := NULL]

  core <- c("sample_id", "chrom", "length")
  data.table::setcolorder(x, c(core, setdiff(names(x), core)))
  unique(x)
}

normalize_bwg_meta <- function(meta, data = NULL) {
  stop_if_not(is.list(meta), "`meta` must be a list.")
  meta <- as.list(meta)

  coordinate <- meta$coordinate %||% "1-based closed"
  stop_if_not(
    identical(as.character(coordinate)[1L], "1-based closed"),
    "BwgTrack coordinates must use the '1-based closed' convention."
  )

  mode <- meta$mode %||% if (is.null(data)) "lazy" else "memory"
  mode <- as.character(mode)[1L]
  stop_if_not(
    mode %in% c("lazy", "memory"),
    "`meta$mode` must be 'lazy' or 'memory'."
  )

  backend <- meta$backend %||% "R"
  backend <- as.character(backend)[1L]
  stop_if_not(
    !is.na(backend) && nzchar(backend),
    "`meta$backend` must be a non-empty character value."
  )

  meta$coordinate <- "1-based closed"
  meta$schema_version <- "2"
  meta$backend <- backend
  meta$mode <- mode
  meta
}

subset_bwg_seqinfo <- function(object, sample_ids = NULL, chrom = NULL) {
  seqinfo <- object$seqinfo
  if (is.null(seqinfo)) {
    return(NULL)
  }

  x <- data.table::copy(data.table::as.data.table(seqinfo))
  if (!is.null(sample_ids)) {
    sample_ids <- as.character(sample_ids)
    x <- x[x[["sample_id"]] %in% sample_ids]
  }
  if (!is.null(chrom)) {
    chrom_value <- as.character(chrom)
    x <- x[x[["chrom"]] %in% chrom_value]
  }
  x[]
}

#' Create a BwgTrack object
#'
#' @description
#' Creates a standardized `BwgTrack` object. GeneTrackR schema version 2 keeps
#' sample metadata, optional in-memory signal records, optional chromosome
#' metadata, object metadata, and validation results in stable named slots.
#' Public signal coordinates use 1-based closed intervals.
#'
#' @param samples Sample metadata table containing `sample_id`.
#' @param data Optional in-memory signal table.
#' @param meta Metadata list. `coordinate`, `schema_version`, `backend`, and
#' `mode` are normalized by the constructor.
#' @param validation Validation result list.
#' @param seqinfo Optional chromosome metadata containing `chrom` and `length`.
#' If `sample_id` is absent, the chromosome metadata are shared across all
#' samples. Unknown chromosome lengths can be represented by `NA`.
#' @return A BwgTrack object with `samples`, `data`, `seqinfo`, `meta`, and
#' `validation` slots.
#' @examples
#' samples <- data.frame(sample_id = "sample1")
#' BwgTrack(samples)
#' @export
BwgTrack <- function(samples,
                     data = NULL,
                     meta = list(),
                     validation = make_empty_validation(),
                     seqinfo = NULL) {
  samples <- data.table::copy(data.table::as.data.table(samples))
  stop_if_not("sample_id" %in% names(samples), "`samples` must contain a `sample_id` column.")
  samples[, sample_id := as.character(sample_id)]
  stop_if_not(
    all(!is.na(samples$sample_id) & nzchar(samples$sample_id)),
    "`samples$sample_id` values must be non-missing and non-empty."
  )

  data <- if (is.null(data)) {
    NULL
  } else {
    data.table::copy(data.table::as.data.table(data))
  }
  seqinfo <- normalize_bwg_seqinfo(seqinfo, samples$sample_id)
  meta <- normalize_bwg_meta(meta, data = data)

  structure(
    list(
      samples = samples,
      data = data,
      seqinfo = seqinfo,
      meta = meta,
      validation = validation
    ),
    class = "BwgTrack"
  )
}

#' @export
print.BwgTrack <- function(x, ...) {
  cat("<BwgTrack>\n")
  cat("  samples: ", nrow(x$samples), "\n", sep = "")
  cat("  mode   : ", x$meta$mode %||% "unknown", "\n", sep = "")
  cat("  backend: ", x$meta$backend %||% "R", "\n", sep = "")
  cat("  coordinate: ", x$meta$coordinate %||% "unknown", "\n", sep = "")
  cat("  schema : ", x$meta$schema_version %||% "unknown", "\n", sep = "")

  if ("has_strand" %in% names(x$samples)) {
    n_unstranded <- sum(!x$samples$has_strand, na.rm = TRUE)
    if (n_unstranded > 0L) {
      cat("  strand : ", n_unstranded, " unstranded sample(s)\n", sep = "")
    }
  }

  if (is.null(x$seqinfo)) {
    cat("  seqinfo: unavailable\n")
  } else {
    cat("  seqinfo: ", nrow(x$seqinfo), " sample/chromosome record(s)\n", sep = "")
  }

  if (is.null(x$data)) {
    cat("  records: not loaded in memory\n")
    cat("  note   : use retrieve_bwg() to read a genomic region\n")
  } else {
    cat("  records: ", nrow(x$data), "\n", sep = "")
  }

  invisible(x)
}

#' Show the first rows of a BwgTrack object
#'
#' @param x A BwgTrack object.
#' @param n Number of rows to show.
#' @param ... Additional arguments.
#' @return Invisibly returns the input object.
#' @export
head.BwgTrack <- function(x, n = 6L, ...) {
  print(x)
  cat("\nSamples:\n")
  print(utils::head(x$samples, n = n))

  if ("has_strand" %in% names(x$samples)) {
    n_unstranded <- sum(!x$samples$has_strand, na.rm = TRUE)
    if (n_unstranded > 0L) {
      cat("  strand : ", n_unstranded, " unstranded sample(s)\n", sep = "")
    }
  }

  if (!is.null(x$seqinfo)) {
    cat("\nSequence metadata:\n")
    print(utils::head(x$seqinfo, n = n))
  }

  if (is.null(x$data)) {
    cat("\nSignal data: NULL because this object is in lazy mode.\n")
    cat("Use retrieve_bwg(object, chrom, start, end) to read signal from a specific region.\n")
  } else {
    cat("\nSignal data:\n")
    print(utils::head(x$data, n = n))
  }

  invisible(x)
}

#' @export
summary.BwgTrack <- function(object, ...) {
  summary_bwg(object, ...)
}
