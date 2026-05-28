# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.4
# Function: Read GTF files into unified Feature annotation objects
# Input: GTF annotation files
# Output: Feature objects

#' Read a GTF file as a FeatureTrack
#'
#' @param file GTF file path.
#' @param feature_types Optional feature types to keep, such as `gene`, `transcript`, `exon`, or `CDS`.
#' @param verbose Whether to print progress messages.
#' @param progress Whether to show a stage progress bar. The parser reports major stages using the same style as `read_genepred()`.
#' @return A FeatureTrack object.
#' @export
read_gtf <- function(file, feature_types = NULL, verbose = TRUE, progress = interactive() && isTRUE(verbose)) {
  read_gff_gtf(file, format = "GTF", feature_types = feature_types, verbose = verbose, progress = progress)
}
