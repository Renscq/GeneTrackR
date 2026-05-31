# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.2
# Function: Define VariantTrack class for genome-level variant data
# Input: VCF-like variant tables
# Output: VariantTrack S3 objects

#' Create a VariantTrack object
#'
#' @description
#' `VariantTrack()` stores VCF-like variant records as genomic point features.
#' Coordinates are stored internally as 1-based positions with `start == end`.
#'
#' @param data A data.frame or data.table containing at least `chrom` and `pos`.
#' @param meta A list of metadata.
#' @return A VariantTrack object.
#' @export
VariantTrack <- function(data = NULL, meta = list()) {
  if (is.null(data)) {
    meta$lazy <- isTRUE(meta$lazy)
    return(structure(
      list(data = NULL, meta = meta),
      class = "VariantTrack"
    ))
  }

  dt <- data.table::as.data.table(data)
  stop_if_not(all(c("chrom", "pos") %in% names(dt)), "`data` must contain `chrom` and `pos` columns.")

  if (!"variant_id" %in% names(dt)) dt[, "variant_id" := paste0("variant_", seq_len(.N))]
  if (!"ref" %in% names(dt)) dt[, "ref" := NA_character_]
  if (!"alt" %in% names(dt)) dt[, "alt" := NA_character_]
  if (!"qual" %in% names(dt)) dt[, "qual" := NA_real_]
  if (!"filter" %in% names(dt)) dt[, "filter" := NA_character_]
  if (!"info" %in% names(dt)) dt[, "info" := NA_character_]
  if (!"variant_type" %in% names(dt)) dt[, "variant_type" := infer_variant_type(dt[["ref"]], dt[["alt"]])]

  dt[, "chrom" := as.character(dt[["chrom"]])]
  dt[, "pos" := as.integer(dt[["pos"]])]
  dt[, "start" := as.integer(dt[["pos"]])]
  dt[, "end" := as.integer(dt[["pos"]])]
  dt[, "variant_id" := as.character(dt[["variant_id"]])]
  dt[, "ref" := as.character(dt[["ref"]])]
  dt[, "alt" := as.character(dt[["alt"]])]
  dt[, "qual" := suppressWarnings(as.numeric(dt[["qual"]]))]
  dt[, "filter" := as.character(dt[["filter"]])]
  dt[, "info" := as.character(dt[["info"]])]
  dt[, "variant_type" := as.character(dt[["variant_type"]])]

  data.table::setorderv(dt, c("chrom", "pos", "variant_id"))
  meta$lazy <- FALSE
  structure(
    list(data = dt, meta = meta),
    class = "VariantTrack"
  )
}

make_lazy_vcf_track <- function(file,
                                keep_genotype = TRUE,
                                verbose = TRUE,
                                progress = interactive() && isTRUE(verbose)) {
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  stop_if_not(has_vcf_tabix_index(file), "Lazy VCF mode requires a bgzip-compressed VCF with a .tbi index.")

  progress_msg <- vcf_progress_message(isTRUE(verbose) && isTRUE(progress))
  progress_msg(1L, 2L, "Reading VCF header.")
  header <- read_vcf_header_line(file)
  col_names <- parse_vcf_header_names(header)

  sample_names <- character()
  if (isTRUE(keep_genotype) && length(col_names) > 9L) {
    sample_names <- col_names[10L:length(col_names)]
  }

  progress_msg(2L, 2L, "Created lazy VariantTrack.")

  VariantTrack(
    data = NULL,
    meta = list(
      source_file = file,
      format = "VCF",
      coordinate_internal = "1-based position",
      lazy = TRUE,
      keep_genotype = isTRUE(keep_genotype),
      sample_names = sample_names,
      header_names = col_names,
      indexed = TRUE
    )
  )
}

is_lazy_variant_track <- function(x) {
  inherits(x, "VariantTrack") && is.null(x$data) && isTRUE(x$meta$lazy)
}

#' @export
print.VariantTrack <- function(x, ...) {
  cat("<VariantTrack>\n")
  if (is_lazy_variant_track(x)) {
    cat("  mode      : lazy\n")
    cat("  source    :", x$meta$source_file %||% "unknown", "\n")
    cat("  variants  : not loaded in memory\n")
    cat("  samples   :", length(x$meta$sample_names %||% character()), "\n")
    cat("  format    :", x$meta$format %||% "VCF", "\n")
    cat("  coordinate:", x$meta$coordinate_internal %||% "1-based position", "\n")
    cat("  note      : use retrieve_vcf(object, chrom, start, end) to read a region\n")
    return(invisible(x))
  }
  cat("  variants  :", format(nrow(x$data), big.mark = ","), "\n")
  cat("  format    :", x$meta$format %||% "VCF", "\n")
  cat("  coordinate:", x$meta$coordinate_internal %||% "1-based position", "\n")
  invisible(x)
}

infer_variant_type <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  first_alt <- sub(",.*$", "", alt)
  ref_n <- nchar(ref)
  alt_n <- nchar(first_alt)
  out <- rep("variant", length(ref))
  out[!is.na(ref_n) & !is.na(alt_n) & ref_n == 1L & alt_n == 1L] <- "SNP"
  out[!is.na(ref_n) & !is.na(alt_n) & ref_n < alt_n] <- "INS"
  out[!is.na(ref_n) & !is.na(alt_n) & ref_n > alt_n] <- "DEL"
  out[!is.na(ref_n) & !is.na(alt_n) & ref_n == alt_n & ref_n > 1L] <- "MNV"
  out
}

#' @export
summary.VariantTrack <- function(object, ...) {
  summary_vcf(object, ...)
}
