# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.5
# Function: BwgTrack S3 class constructors and print methods
# Input: Signal sample metadata and optional in-memory signal table
# Output: BwgTrack object

#' Create a BwgTrack object
#'
#' @param samples Sample metadata table.
#' @param data Optional in-memory signal table.
#' @param meta Metadata list.
#' @param validation Validation result list.
#' @return A BwgTrack object.
#' @export
BwgTrack <- function(samples, data = NULL, meta = list(), validation = make_empty_validation()) {
  structure(
    list(
      samples = data.table::as.data.table(samples),
      data = if (is.null(data)) NULL else data.table::as.data.table(data),
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

  if ("has_strand" %in% names(x$samples)) {
    n_unstranded <- sum(!x$samples$has_strand, na.rm = TRUE)
    if (n_unstranded > 0L) {
      cat("  strand : ", n_unstranded, " unstranded sample(s)
", sep = "")
    }
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
      cat("  strand : ", n_unstranded, " unstranded sample(s)
", sep = "")
    }
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
