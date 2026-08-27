# Author: Rensc
# Date: 2026-05-28
# Version: dev002
# Function: Read GFF3 files into unified Feature annotation objects
# Input: GFF3 annotation files
# Output: Feature objects

#' Read a GFF3 file as a FeatureTrack
#'
#' @param file GFF3 file path.
#' @param feature_types Optional feature types to keep, such as `gene`, `mRNA`, `exon`, or `CDS`.
#' @param verbose Whether to print progress messages.
#' @param progress Whether to show a stage progress bar. The parser reports major stages using the same style as `read_genepred()`.
#' @return A FeatureTrack object.
#' @examples
#' gff_file <- system.file("extdata", "gtr_demo.gff3", package = "GeneTrackR")
#' gff <- read_gff(gff_file, verbose = FALSE, progress = FALSE)
#' gff
#' head(gff$genes)
#' @export
read_gff <- function(file, feature_types = NULL, verbose = TRUE, progress = interactive() && isTRUE(verbose)) {
  read_gff_gtf(file, format = "GFF", feature_types = feature_types, verbose = verbose, progress = progress)
}
