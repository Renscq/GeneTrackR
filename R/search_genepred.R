# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.26
# Function: Search genes and transcripts in GenePred objects
# Input: GenePred object and search pattern
# Output: data.table search results

#' Search genes in a GenePred object
#'
#' @param object A GenePred object.
#' @param pattern Character pattern used to search gene IDs. Regular expressions are supported by default.
#' @param chrom Optional chromosome filter.
#' @param ignore_case Logical. Whether to ignore case when matching.
#' @param fixed Logical. Whether `pattern` should be treated as a fixed string instead of a regular expression.
#'
#' @return A data.table with gene coordinates and transcript counts.
#'
#' @examples
#' \dontrun{
#' gp <- read_genepred(system.file("extdata", "example.genePredExt", package = "GeneTrackR"), format = "genePredExt")
#' search_gene(gp, "GeneA")
#' }
#' @export
search_gene <- function(object, pattern, chrom = NULL, ignore_case = TRUE, fixed = FALSE) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  stop_if_not(length(pattern) == 1L && !is.na(pattern) && nzchar(pattern), "`pattern` must be a non-empty character string.")
  genes <- data.table::copy(data.table::as.data.table(object$genes))
  if (!is.null(chrom)) {
    genes <- genes[genes[["chrom"]] %in% as.character(chrom)]
  }
  if (nrow(genes) == 0L) {
    return(genes)
  }
  hit <- grepl(pattern, genes[["gene_id"]], ignore.case = ignore_case, fixed = fixed)
  out <- genes[hit]
  data.table::setorderv(out, c("chrom", "gene_start", "gene_end", "gene_id"))
  out[]
}

#' Search transcripts in a GenePred object
#'
#' @param object A GenePred object.
#' @param pattern Character pattern used to search transcript IDs. Regular expressions are supported by default.
#' @param chrom Optional chromosome filter.
#' @param gene_id Optional gene ID filter.
#' @param ignore_case Logical. Whether to ignore case when matching.
#' @param fixed Logical. Whether `pattern` should be treated as a fixed string instead of a regular expression.
#'
#' @return A data.table with transcript coordinates and annotation fields.
#'
#' @examples
#' \dontrun{
#' gp <- read_genepred(system.file("extdata", "example.genePredExt", package = "GeneTrackR"), format = "genePredExt")
#' search_transcript(gp, "TxA")
#' }
#' @export
search_transcript <- function(object, pattern, chrom = NULL, gene_id = NULL, ignore_case = TRUE, fixed = FALSE) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  stop_if_not(length(pattern) == 1L && !is.na(pattern) && nzchar(pattern), "`pattern` must be a non-empty character string.")
  tx <- data.table::copy(data.table::as.data.table(object$transcripts))
  if (!is.null(chrom)) {
    tx <- tx[tx[["chrom"]] %in% as.character(chrom)]
  }
  if (!is.null(gene_id)) {
    tx <- tx[tx[["gene_id"]] %in% as.character(gene_id)]
  }
  if (nrow(tx) == 0L) {
    return(tx)
  }
  hit <- grepl(pattern, tx[["transcript_id"]], ignore.case = ignore_case, fixed = fixed)
  out <- tx[hit]
  data.table::setorderv(out, c("chrom", "tx_start", "tx_end", "transcript_id"))
  out[]
}
