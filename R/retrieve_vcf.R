# Author: Rensc
# Date: 2026-05-29
# Version: 0.3.1
# Function: Retrieve variants from VariantTrack objects or indexed VCF files
# Input: VariantTrack object or VCF path and filters
# Output: data.table or VariantTrack object

#' Retrieve variants from a VariantTrack or an indexed VCF file
#'
#' @description
#' Retrieves VCF records by genomic region, variant ID, type, or text pattern.
#' For large bgzip-compressed VCF files with a tabix index, `object` can be a
#' file path and the function will query only the requested region using
#' `Rsamtools::scanTabix()`.
#'
#' @param object A VariantTrack object or a path to a VCF/VCF.GZ file.
#' @param pattern Optional text pattern matched against variant ID, REF, ALT, INFO, and variant type.
#' @param chrom Optional chromosome name. Required for indexed file region queries.
#' @param start Optional 1-based region start.
#' @param end Optional 1-based region end.
#' @param variant_id Optional variant ID vector.
#' @param variant_type Optional variant type vector, such as `SNP`, `INS`, `DEL`, or `MNV`.
#' @param ignore_case Logical. Whether pattern matching ignores case.
#' @param fixed Logical. Whether pattern is matched as a fixed string.
#' @param as Output type: `data.table` or `VariantTrack`.
#' @return A data.table or VariantTrack object.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' retrieve_vcf(vcf, chrom = "chr1", start = 100, end = 600)
#' retrieve_vcf(vcf, variant_type = "SNP")
#' retrieve_vcf(vcf, pattern = "rsA", fixed = TRUE)
#' vt <- retrieve_vcf(vcf, chrom = "chr1", start = 100, end = 600, as = "VariantTrack")
#' vt
#' @export
retrieve_vcf <- function(object,
                         pattern = NULL,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         variant_id = NULL,
                         variant_type = NULL,
                         ignore_case = TRUE,
                         fixed = FALSE,
                         as = c("data.table", "VariantTrack")) {
  as <- match.arg(as)

  if (is.character(object) && length(object) == 1L && file.exists(object)) {
    vt <- retrieve_vcf_file(object, chrom = chrom, start = start, end = end)
  } else {
    stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object or a VCF file path.")
    vt <- object
  }

  dt <- data.table::copy(vt$data)

  query_chrom <- as.character(chrom)
  query_chrom <- query_chrom[!is.na(query_chrom) & nzchar(query_chrom)]
  if (length(query_chrom) > 0L) {
    keep <- as.character(dt[["chrom"]]) %in% query_chrom
    dt <- dt[keep]
  }

  if (!is.null(start) && !is.null(end)) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    dt <- dt[dt[["pos"]] >= s & dt[["pos"]] <= e]
  }

  if (!is.null(variant_id)) {
    ids <- as.character(variant_id)
    dt <- dt[dt[["variant_id"]] %in% ids]
  }

  if (!is.null(variant_type)) {
    types <- as.character(variant_type)
    dt <- dt[dt[["variant_type"]] %in% types]
  }

  dt <- match_pattern_internal(
    dt,
    pattern = pattern,
    fields = c("variant_id", "ref", "alt", "info", "variant_type"),
    ignore_case = ignore_case,
    fixed = fixed
  )

  if (nrow(dt) > 0L) {
    data.table::setorderv(dt, intersect(c("chrom", "pos", "variant_id"), names(dt)))
  }

  if (as == "data.table") return(dt[])
  VariantTrack(dt, meta = vt$meta)
}

retrieve_vcf_file <- function(file, chrom = NULL, start = NULL, end = NULL) {
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  is_indexed <- has_vcf_tabix_index(file)

  if (isTRUE(is_indexed) && !is.null(chrom) && !is.null(start) && !is.null(end)) {
    stop_if_not(requireNamespace("Rsamtools", quietly = TRUE), "Package `Rsamtools` is required for indexed VCF queries.")
    region <- paste0(as.character(chrom)[1L], ":", as.integer(start)[1L], "-", as.integer(end)[1L])
    tbx <- Rsamtools::TabixFile(file)
    lines <- Rsamtools::scanTabix(tbx, param = region)[[1L]]
    header <- read_vcf_header_line(file)
    col_names <- parse_vcf_header_names(header)
    if (length(lines) == 0L) {
      dt <- data.table::data.table(matrix(ncol = length(col_names), nrow = 0L))
      data.table::setnames(dt, col_names)
    } else {
      dt <- data.table::fread(
        text = paste(c(paste(col_names, collapse = "\t"), lines), collapse = "\n"),
        sep = "\t",
        header = TRUE,
        data.table = TRUE,
        showProgress = FALSE
      )
    }
    return(vcf_table_to_variant_track(dt, file))
  }

  read_vcf(file, verbose = FALSE)
}

vcf_table_to_variant_track <- function(dt, source_file = NULL) {
  if ("#CHROM" %in% names(dt) && !"CHROM" %in% names(dt)) {
    data.table::setnames(dt, "#CHROM", "CHROM")
  }
  standard_in <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  stop_if_not(all(standard_in %in% names(dt)), "A VCF table must contain standard VCF columns.")
  sample_cols <- setdiff(names(dt), c(standard_in, "FORMAT"))
  has_format <- "FORMAT" %in% names(dt)
  data.table::setnames(
    dt,
    old = standard_in,
    new = c("chrom", "pos", "variant_id", "ref", "alt", "qual", "filter", "info")
  )
  dt[is.na(variant_id) | variant_id == "." | variant_id == "", "variant_id" := paste0(as.character(chrom), ":", as.integer(pos), ":", as.character(ref), ":", as.character(alt))]
  VariantTrack(
    dt,
    meta = list(
      source_file = source_file,
      format = "VCF",
      coordinate_input = "1-based position",
      coordinate_internal = "1-based position",
      indexed = if (!is.null(source_file)) has_vcf_tabix_index(source_file) else FALSE,
      has_genotype = length(sample_cols) > 0L,
      has_format = isTRUE(has_format),
      sample_names = sample_cols
    )
  )
}
