# Author: Rensc
# Date: 2026-08-31
# Version: dev001
# Function: Convert GeneTrackR genomic track objects to GenomicRanges
# Input: Annotation, signal, or variant track objects
# Output: GenomicRanges::GRanges objects

set_granges_metadata <- function(gr, metadata) {
  metadata <- data.table::copy(data.table::as.data.table(metadata))
  if (ncol(metadata) == 0L) {
    return(gr)
  }

  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    metadata,
    check.names = FALSE
  )
  gr
}

bwg_granges_seqlengths <- function(object, chrom) {
  if (is.null(object$seqinfo)) {
    return(NULL)
  }

  x <- data.table::copy(data.table::as.data.table(object$seqinfo))
  chrom <- unique(as.character(chrom))
  chrom <- chrom[!is.na(chrom) & nzchar(chrom)]
  if (length(chrom) == 0L) {
    return(NULL)
  }

  x <- x[x[["chrom"]] %in% chrom]
  if (nrow(x) == 0L) {
    return(NULL)
  }

  out <- vapply(chrom, function(chrom_value) {
    values <- unique(x[x[["chrom"]] == chrom_value][["length"]])
    values <- values[!is.na(values)]
    if (length(values) == 1L) as.integer(values[[1L]]) else NA_integer_
  }, integer(1L))
  names(out) <- chrom
  out
}

annotation_as_granges <- function(object, level) {
  if (!inherits(object, "GenePred") && level != "feature") {
    object <- as_genepred(object)
  }

  if (level == "feature") {
    dt <- as_feature_table(object)
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(start = dt$start, end = dt$end),
      strand = dt$strand,
      feature_id = dt$feature_id,
      name = dt$name,
      type = dt$type,
      gene_id = dt$gene_id,
      transcript_id = dt$transcript_id
    ))
  }

  if (level == "gene") {
    dt <- as_gene_table(object)
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(start = dt$gene_start, end = dt$gene_end),
      strand = dt$strand,
      gene_id = dt$gene_id,
      gene_type = dt$gene_type
    ))
  }

  if (level == "transcript") {
    dt <- as_transcript_table(object)
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(start = dt$tx_start, end = dt$tx_end),
      strand = dt$strand,
      transcript_id = dt$transcript_id,
      gene_id = dt$gene_id,
      gene_type = dt$gene_type
    ))
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

bwg_as_granges <- function(object) {
  stop_if_not(
    !is.null(object$data),
    "Cannot convert a lazy BwgTrack directly to GRanges; use retrieve_bwg(..., as = 'GRanges') for a genomic region first."
  )

  dt <- data.table::copy(data.table::as.data.table(object$data))
  stop_if_not(
    all(c("chrom", "start", "end", "value") %in% names(dt)),
    "BwgTrack data must contain chrom, start, end, and value columns."
  )
  if (!"sample_id" %in% names(dt)) dt[, sample_id := NA_character_]
  if (!"strand" %in% names(dt)) dt[, strand := "*"]

  seqlengths <- bwg_granges_seqlengths(object, dt$chrom)
  gr <- GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges = IRanges::IRanges(start = dt$start, end = dt$end),
    strand = dt$strand,
    seqlengths = seqlengths
  )
  set_granges_metadata(
    gr,
    data.table::data.table(
      sample_id = as.character(dt$sample_id),
      score = as.numeric(dt$value)
    )
  )
}

variant_as_granges <- function(object) {
  stop_if_not(
    !is.null(object$data),
    "Cannot convert a lazy VariantTrack directly to GRanges; use retrieve_vcf(..., as = 'GRanges') for a genomic region first."
  )

  dt <- data.table::copy(data.table::as.data.table(object$data))
  stop_if_not(
    all(c("chrom", "pos") %in% names(dt)),
    "VariantTrack data must contain chrom and pos columns."
  )

  gr <- GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges = IRanges::IRanges(start = dt$pos, end = dt$pos),
    strand = rep("*", nrow(dt))
  )
  metadata_columns <- setdiff(names(dt), c("chrom", "start", "end"))
  set_granges_metadata(gr, dt[, metadata_columns, with = FALSE])
}

#' Convert a GeneTrackR genomic track object to GRanges
#'
#' @description
#' Converts GeneTrackR annotation, signal, and variant objects to
#' `GenomicRanges::GRanges` so they can be used directly with Bioconductor
#' genomic interval infrastructure. Annotation objects support gene,
#' transcript, exon, and feature levels. In-memory `BwgTrack` objects expose
#' `sample_id` and signal `score` metadata, while `VariantTrack` objects retain
#' their VCF-derived metadata columns. Lazy signal or variant objects should be
#' retrieved by genomic region before conversion.
#'
#' @param object A GenePred, Feature, FeatureTrack, BwgTrack, or VariantTrack object.
#' @param level Annotation feature level. One of `gene`, `transcript`, `exon`,
#' or `feature`. Ignored for BwgTrack and VariantTrack objects.
#' @return A `GenomicRanges::GRanges` object.
#' @examples
#' gp_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
#' gp <- read_genepred(gp_file, format = "genePredExt", verbose = FALSE)
#' as_granges(gp, level = "gene")
#'
#' signal <- BwgTrack(
#'   samples = data.frame(sample_id = "sampleA"),
#'   data = data.frame(
#'     sample_id = "sampleA", chrom = "chr1", start = 10L, end = 12L,
#'     value = 1.5, strand = "*"
#'   )
#' )
#' as_granges(signal)
#'
#' vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
#' as_granges(vcf)
#' @export
as_granges <- function(object, level = c("gene", "transcript", "exon", "feature")) {
  if (inherits(object, "BwgTrack")) {
    return(bwg_as_granges(object))
  }
  if (inherits(object, "VariantTrack")) {
    return(variant_as_granges(object))
  }

  stop_if_not(
    inherits(object, "GenePred") || inherits(object, "Feature") || inherits(object, "FeatureTrack"),
    "`object` must be a GenePred, Feature, FeatureTrack, BwgTrack, or VariantTrack object."
  )
  level <- match.arg(level)
  annotation_as_granges(object, level = level)
}
