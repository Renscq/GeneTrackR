# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.26
# Function: Search genes, transcripts, and standardized features
# Input: GenePred object and search pattern
# Output: data.table search results

#' Search genes in a GenePred object
#'
#' @param object A Feature-compatible annotation object.
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
  stop_if_not(is_gene_model_feature(object), "`object` must be a Feature-compatible annotation object.")
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
#' @param object A Feature-compatible annotation object.
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
  stop_if_not(is_gene_model_feature(object), "`object` must be a Feature-compatible annotation object.")
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


#' Search standardized annotation features
#'
#' @param object A Feature-compatible annotation object.
#' @param pattern Character pattern used to search `feature_id`, `name`, `gene_id`, or `transcript_id`.
#' @param chrom Optional chromosome filter.
#' @param type Optional feature type filter.
#' @param level Optional standardized feature level filter.
#' @param ignore_case Logical. Whether to ignore case when matching.
#' @param fixed Logical. Whether `pattern` should be treated as a fixed string.
#' @return A data.table with matching standardized feature records.
#' @export
search_feature <- function(object, pattern, chrom = NULL, type = NULL, level = NULL, ignore_case = TRUE, fixed = FALSE) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"), "`object` must be a Feature-compatible annotation object.")
  stop_if_not(length(pattern) == 1L && !is.na(pattern) && nzchar(pattern), "`pattern` must be a non-empty character string.")
  dt <- as_feature_table(object)
  if (!is.null(chrom)) dt <- dt[dt[["chrom"]] %in% as.character(chrom)]
  if (!is.null(type)) dt <- dt[dt[["type"]] %in% as.character(type)]
  if (!is.null(level)) dt <- dt[dt[["level"]] %in% as.character(level)]
  if (nrow(dt) == 0L) return(dt)
  fields <- c("feature_id", "name", "gene_id", "transcript_id")
  fields <- intersect(fields, names(dt))
  hit <- rep(FALSE, nrow(dt))
  for (field in fields) {
    hit <- hit | grepl(pattern, as.character(dt[[field]]), ignore.case = ignore_case, fixed = fixed)
  }
  out <- dt[hit]
  data.table::setorderv(out, c("chrom", "start", "end", "feature_id"))
  out[]
}
