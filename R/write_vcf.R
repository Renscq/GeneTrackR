# Author: Rensc
# Date: 2026-05-28
# Version: dev001
# Function: Write VariantTrack objects to VCF files
# Input: VariantTrack object
# Output: VCF file

#' Write a VariantTrack object to VCF
#'
#' @param object A `VariantTrack` object created by `read_vcf()`.
#' @param file Output VCF file path.
#' @param overwrite Whether to overwrite an existing file.
#' @param include_header Whether to include minimal VCF header lines.
#' @return Invisibly returns the output file path.
#' @examples
#' \dontrun{
#' vcf <- read_vcf("variants.vcf")
#' write_vcf(vcf, "variants.out.vcf")
#' }
#' @export
write_vcf <- function(object, file, overwrite = FALSE, include_header = TRUE) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  check_output_file(file, overwrite)
  dt <- data.table::copy(object$data)
  out <- dt[, .(
    `#CHROM` = as.character(chrom),
    POS = as.integer(pos),
    ID = ifelse(is.na(variant_id) | variant_id == "", ".", as.character(variant_id)),
    REF = ifelse(is.na(ref) | ref == "", ".", as.character(ref)),
    ALT = ifelse(is.na(alt) | alt == "", ".", as.character(alt)),
    QUAL = ifelse(is.na(qual), ".", as.character(qual)),
    FILTER = ifelse(is.na(filter) | filter == "", ".", as.character(filter)),
    INFO = ifelse(is.na(info) | info == "", ".", as.character(info))
  )]
  if (isTRUE(include_header)) {
    writeLines(c("##fileformat=VCFv4.2", "##source=GeneTrackR"), con = file)
    data.table::fwrite(out, file, sep = "\t", col.names = TRUE, append = TRUE)
  } else {
    data.table::fwrite(out, file, sep = "\t", col.names = TRUE)
  }
  invisible(file)
}

#' Write a VariantTrack object
#'
#' @description Backward-compatible wrapper around `write_vcf()`.
#' @inheritParams write_vcf
#' @export
write_variant_track <- function(object, file, overwrite = FALSE, include_header = TRUE) {
  write_vcf(object = object, file = file, overwrite = overwrite, include_header = include_header)
}
