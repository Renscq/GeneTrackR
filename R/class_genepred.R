# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.0
# Function: GenePred S3 class constructors and coercion helpers
# Input: Gene, transcript, and exon tables
# Output: GenePred object and derived tables

#' Create a GenePred object
#'
#' @param transcripts Transcript-level annotation table.
#' @param exons Exon-level annotation table.
#' @param genes Gene-level annotation table.
#' @param meta Metadata list.
#' @param validation Validation result list.
#' @return A GenePred object.
#' @export
GenePred <- function(transcripts, exons, genes = NULL, meta = list(), validation = make_empty_validation()) {
  if (is.null(genes)) {
    genes <- build_gene_table(transcripts)
  }
  structure(
    list(
      transcripts = data.table::as.data.table(transcripts),
      exons = data.table::as.data.table(exons),
      genes = data.table::as.data.table(genes),
      meta = meta,
      validation = validation
    ),
    class = "GenePred"
  )
}

#' @export
print.GenePred <- function(x, ...) {
  cat("<GenePred>\n")
  cat("  genes      : ", nrow(x$genes), "\n", sep = "")
  cat("  transcripts: ", nrow(x$transcripts), "\n", sep = "")
  cat("  exons      : ", nrow(x$exons), "\n", sep = "")
  cat("  coordinate : ", x$meta$coordinate_internal %||% "1-based closed", "\n", sep = "")
  invisible(x)
}

#' @export
summary.GenePred <- function(object, ...) {
  summary_genepred(object, ...)
}

#' Extract transcript table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with transcript-level annotation.
#' @export
as_transcript_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$transcripts)
}

#' Extract exon table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with exon-level annotation.
#' @export
as_exon_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$exons)
}

#' Extract gene table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with gene-level annotation.
#' @export
as_gene_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$genes)
}

#' Convert a GenePred object to GRanges
#'
#' @param object A GenePred object.
#' @param level Feature level. One of gene, transcript, or exon.
#' @return A GRanges object.
#' @export
as_granges <- function(object, level = c("gene", "transcript", "exon")) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  level <- match.arg(level)

  if (level == "gene") {
    dt <- as_gene_table(object)
    gr <- GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(start = dt$gene_start, end = dt$gene_end),
      strand = dt$strand,
      gene_id = dt$gene_id,
      gene_type = dt$gene_type
    )
    return(gr)
  }

  if (level == "transcript") {
    dt <- as_transcript_table(object)
    gr <- GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(start = dt$tx_start, end = dt$tx_end),
      strand = dt$strand,
      transcript_id = dt$transcript_id,
      gene_id = dt$gene_id,
      gene_type = dt$gene_type
    )
    return(gr)
  }

  dt <- as_exon_table(object)
  GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges = IRanges::IRanges(start = dt$exon_start, end = dt$exon_end),
    strand = dt$strand,
    transcript_id = dt$transcript_id,
    gene_id = dt$gene_id,
    exon_number = dt$exon_number
  )
}

build_gene_table <- function(transcripts) {
  dt <- data.table::as.data.table(transcripts)
  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      gene_id = character(), chrom = character(), strand = character(),
      gene_start = integer(), gene_end = integer(), n_transcripts = integer(),
      gene_type = character()
    ))
  }
  dt[, .(
    chrom = chrom[1],
    strand = strand[1],
    gene_start = min(tx_start),
    gene_end = max(tx_end),
    n_transcripts = data.table::uniqueN(transcript_id),
    gene_type = if (any(gene_type == "coding")) "coding" else "non-coding"
  ), by = gene_id][order(chrom, gene_start, gene_end, gene_id)]
}
