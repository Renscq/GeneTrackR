# Author: Rensc
# Date: 2026-05-29
# Version: 0.3.1
# Function: Read VCF files into genome-level VariantTrack objects
# Input: VCF files
# Output: VariantTrack objects

#' Read a VCF file as a VariantTrack
#'
#' @description
#' Reads VCF records and stores them as genomic point variants. Standard VCF
#' columns are always imported. If genotype sample columns are present, they are
#' preserved for haplotype analysis.
#'
#' @param file VCF file path. Plain VCF, gzip-compressed VCF, and bgzip VCF are supported.
#' @param keep_genotype Logical. Whether to keep FORMAT and sample genotype columns.
#' @param verbose Logical. Whether to print progress messages.
#' @return A VariantTrack object.
#' @examples
#' \dontrun{
#' variants <- read_vcf("variants.vcf.gz")
#' plot_variant(variants, chrom = "chr1", start = 1, end = 10000)
#' }
#' @export
read_vcf <- function(file, keep_genotype = TRUE, verbose = TRUE) {
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))

  file_path <- normalizePath(file, winslash = "/", mustWork = FALSE)
  if (isTRUE(verbose)) {
    message("[GeneTrackR] Reading VCF file: ", file_path)
  }

  header <- read_vcf_header_line(file)
  col_names <- parse_vcf_header_names(header)

  dt <- data.table::fread(
    file,
    skip = "#CHROM",
    header = TRUE,
    sep = "\t",
    data.table = TRUE,
    showProgress = FALSE
  )

  if (nrow(dt) == 0L && length(col_names) > 0L) {
    dt <- data.table::data.table(matrix(ncol = length(col_names), nrow = 0L))
    data.table::setnames(dt, col_names)
  }

  if (length(col_names) == ncol(dt)) {
    data.table::setnames(dt, col_names)
  } else {
    data.table::setnames(dt, seq_len(min(8L, ncol(dt))), c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")[seq_len(min(8L, ncol(dt)))])
  }

  if ("#CHROM" %in% names(dt) && !"CHROM" %in% names(dt)) {
    data.table::setnames(dt, "#CHROM", "CHROM")
  }
  standard_in <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  stop_if_not(all(standard_in %in% names(dt)), "A VCF file must contain standard VCF columns.")

  sample_cols <- setdiff(names(dt), c(standard_in, "FORMAT"))
  has_format <- "FORMAT" %in% names(dt)

  data.table::setnames(
    dt,
    old = standard_in,
    new = c("chrom", "pos", "variant_id", "ref", "alt", "qual", "filter", "info")
  )

  dt[is.na(variant_id) | variant_id == "." | variant_id == "", "variant_id" := paste0(as.character(chrom), ":", as.integer(pos), ":", as.character(ref), ":", as.character(alt))]

  if (!isTRUE(keep_genotype)) {
    keep_cols <- c("chrom", "pos", "variant_id", "ref", "alt", "qual", "filter", "info")
    dt <- dt[, ..keep_cols]
    sample_cols <- character()
    has_format <- FALSE
  }

  vt <- VariantTrack(
    dt,
    meta = list(
      source_file = file_path,
      format = "VCF",
      coordinate_input = "1-based position",
      coordinate_internal = "1-based position",
      indexed = has_vcf_tabix_index(file_path),
      has_genotype = length(sample_cols) > 0L,
      has_format = isTRUE(has_format),
      sample_names = sample_cols
    )
  )

  if (isTRUE(verbose)) {
    message("[GeneTrackR] Finished. variants: ", format(nrow(vt$data), big.mark = ","), "; samples: ", length(sample_cols), ".")
  }

  vt
}

#' Read a VCF file as a VariantTrack
#'
#' @description Backward-compatible alias of `read_vcf()`.
#' @inheritParams read_vcf
#' @return A VariantTrack object.
#' @export
read_vcf_track <- function(file, keep_genotype = TRUE, verbose = TRUE) {
  read_vcf(file = file, keep_genotype = keep_genotype, verbose = verbose)
}

read_vcf_header_line <- function(file) {
  con <- if (grepl("\\.gz$|\\.bgz$", file, ignore.case = TRUE)) gzfile(file, open = "rt") else file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    x <- readLines(con, n = 1L, warn = FALSE)
    if (length(x) == 0L) break
    if (startsWith(x, "#CHROM")) return(x)
  }
  stop("VCF header line starting with #CHROM was not found.", call. = FALSE)
}

parse_vcf_header_names <- function(header) {
  header <- sub("^#", "", header)
  strsplit(header, "\t", fixed = TRUE)[[1L]]
}

has_vcf_tabix_index <- function(file) {
  file.exists(paste0(file, ".tbi")) || file.exists(sub("\\.gz$", ".tbi", file, ignore.case = TRUE))
}
