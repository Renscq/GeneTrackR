# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.2
# Function: Read VCF files into genome-level VariantTrack objects
# Input: VCF files
# Output: VariantTrack objects

#' Read a VCF file as a VariantTrack
#'
#' @description
#' Reads VCF records and stores them as genomic point variants. Genotype sample
#' columns are not parsed; this function focuses on track visualization.
#'
#' @param file VCF file path. Gzip-compressed files are supported by `data.table::fread()`.
#' @param verbose Whether to print progress messages.
#' @return A VariantTrack object.
#' @examples
#' \dontrun{
#' variants <- read_vcf("variants.vcf.gz")
#' plot_variant_track(variants, chrom = "chr1", start = 1, end = 10000)
#' }
#' @export
read_vcf <- function(file, verbose = TRUE) {
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))
  if (isTRUE(verbose)) message("[GeneTrackR] Reading VCF file: ", normalizePath(file, winslash = "/", mustWork = FALSE))

  dt <- data.table::fread(
    file,
    header = FALSE,
    sep = "\t",
    data.table = TRUE,
    comment.char = "#",
    showProgress = isTRUE(verbose)
  )
  stop_if_not(ncol(dt) >= 8L, "A VCF file must contain at least 8 columns.")
  data.table::setnames(dt, seq_len(8L), c("chrom", "pos", "variant_id", "ref", "alt", "qual", "filter", "info"))
  dt[is.na(variant_id) | variant_id == "." | variant_id == "", "variant_id" := paste0(as.character(chrom), ":", as.integer(pos), ":", as.character(ref), ":", as.character(alt))]
  VariantTrack(
    dt[, .(
      chrom = as.character(chrom),
      pos = as.integer(pos),
      variant_id = as.character(variant_id),
      ref = as.character(ref),
      alt = as.character(alt),
      qual = suppressWarnings(as.numeric(qual)),
      filter = as.character(filter),
      info = as.character(info)
    )],
    meta = list(
      source_file = normalizePath(file, winslash = "/", mustWork = FALSE),
      format = "VCF",
      coordinate_input = "1-based position",
      coordinate_internal = "1-based position"
    )
  )
}

#' Read a VCF file as a VariantTrack
#'
#' @description Backward-compatible alias of `read_vcf()`.
#' @inheritParams read_vcf
#' @return A VariantTrack object.
#' @export
read_vcf_track <- function(file, verbose = TRUE) {
  read_vcf(file = file, verbose = verbose)
}
