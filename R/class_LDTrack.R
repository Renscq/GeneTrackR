# Author: Rensc
# Date: 2026-08-31
# Version: dev002
# Function: LDTrack class definition
# Input: LD calculation outputs
# Output: LDTrack S3 objects

#' Create an LDTrack object
#'
#' @description
#' `LDTrack()` stores pairwise linkage disequilibrium results, retained variant
#' metadata, the selected genomic region, and optionally the dosage matrix and
#' the generated figure. It follows the same lightweight S3 list-class style used
#' by `VariantTrack` and `HapVariant`.
#'
#' @param data A data.frame/data.table of pairwise LD statistics.
#' @param matrix A symmetric LD matrix. Diagonal values are usually 1.
#' @param variants Variant metadata used in the calculation.
#' @param region A list with `chrom`, `start`, and `end`.
#' @param genotype Optional numeric genotype dosage matrix.
#' @param figure Optional figure object produced by `plot_ld_block()`.
#' @param meta A metadata list.
#' @param plot Deprecated alias for `figure`, kept only for compatibility with
#' older development versions.
#' @return An LDTrack object.
#' @examples
#' LDTrack(region = list(chrom = "chr1", start = 1L, end = 100L))
#' @export
LDTrack <- function(data = NULL,
                    matrix = NULL,
                    variants = NULL,
                    region = list(),
                    genotype = NULL,
                    figure = NULL,
                    meta = list(),
                    plot = NULL) {
  if (!is.null(plot) && is.null(figure)) {
    figure <- plot
  }

  ld_data <- if (is.null(data)) data.table::data.table() else data.table::as.data.table(data)
  var_data <- if (is.null(variants)) data.table::data.table() else data.table::as.data.table(variants)

  if (nrow(var_data) > 0L) {
    stop_if_not(all(c("chrom", "pos", "variant_id") %in% names(var_data)),
                "`variants` must contain `chrom`, `pos`, and `variant_id` columns.")
    if (!"variant_index" %in% names(var_data)) {
      var_data[, "variant_index" := seq_len(.N)]
    }
    data.table::setorderv(var_data, intersect(c("chrom", "pos", "variant_id"), names(var_data)))
    var_data[, "variant_index" := seq_len(.N)]
  }

  structure(
    list(
      data = ld_data[],
      matrix = matrix,
      variants = var_data[],
      region = region,
      genotype = genotype,
      figure = figure,
      meta = meta
    ),
    class = "LDTrack"
  )
}

#' @export
print.LDTrack <- function(x, ...) {
  cat("<LDTrack>\n")
  region <- x$region %||% list()
  if (!is.null(region$chrom) && !is.null(region$start) && !is.null(region$end)) {
    cat("  region    : ", region$chrom, ":", region$start, "-", region$end, "\n", sep = "")
  }
  cat("  method    : ", x$meta$method %||% "unknown", "\n", sep = "")
  cat("  variants  : ", format(nrow(x$variants), big.mark = ","), "\n", sep = "")
  cat("  pairs     : ", format(nrow(x$data), big.mark = ","), "\n", sep = "")
  cat("  samples   : ", format(x$meta$sample_n %||% 0L, big.mark = ","), "\n", sep = "")
  cat("  figure    : ", if (!is.null(x$figure)) "yes" else "no", "\n", sep = "")
  invisible(x)
}
