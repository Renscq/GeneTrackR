# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Adapt typed Bioconductor VCFs without replacing native variant I/O.
# Input: A CollapsedVCF or a VariantTrack created from one source VCF.
# Output: A POS-based analysis track or a typed VCF with matching rows and samples.

vcf_row_key <- ".gtr_vcf_row"

require_variant_annotation <- function() {
  stop_if_not(requireNamespace("VariantAnnotation", quietly = TRUE),
    "Install the Bioconductor package 'VariantAnnotation' to use the VCF adapters.")
}

vcf_with_samples <- function(object, sample_names) {
  header <- VariantAnnotation::header(object)
  if (!identical(VariantAnnotation::samples(header), sample_names)) {
    VariantAnnotation::header(object) <- VariantAnnotation::VCFHeader(
      reference = VariantAnnotation::reference(header), samples = sample_names,
      header = VariantAnnotation::header(header)
    )
  }
  if (!length(sample_names) && length(VariantAnnotation::geno(object))) {
    VariantAnnotation::geno(object) <- S4Vectors::SimpleList()
  }
  object
}

vcf_analysis_table <- function(object, fields, sample_names, row_keys, ranges) {
  if (nrow(object) == 0L) {
    columns <- c(fields, if (length(sample_names)) c("FORMAT", sample_names))
    return(data.table::as.data.table(stats::setNames(rep(list(character()), length(columns)), columns)))
  }

  # Project fields with the standard writer; keep typed values in the source.
  ids <- rownames(object)
  rownames(object) <- row_keys
  file <- tempfile(fileext = ".vcf")
  on.exit(unlink(file), add = TRUE)
  VariantAnnotation::writeVcf(object, file)
  header_line <- read_vcf_header_line(file)
  table <- data.table::fread(file = file, skip = header_line, header = TRUE,
    sep = "\t", quote = "", colClasses = "character", na.strings = NULL,
    check.names = FALSE, showProgress = FALSE)
  expected <- c("#CHROM", fields[-1L], if (length(sample_names)) c("FORMAT", sample_names))
  stop_if_not(identical(names(table), expected),
    "The standard VCF writer did not produce the expected site and sample columns.")
  index <- match(row_keys, table[["ID"]])
  stop_if_not(nrow(table) == length(row_keys) && !anyNA(index),
    "The standard VCF writer did not preserve the source rows.")
  table <- table[index]
  stop_if_not(identical(table[["#CHROM"]], as.character(GenomicRanges::seqnames(ranges))) &&
    identical(table[["POS"]], as.character(IRanges::start(ranges))),
    "The standard VCF writer did not preserve source coordinates.")
  if (is.null(ids)) ids <- rep.int(".", nrow(table))
  data.table::set(table, j = "ID", value = ids)
  table
}

#' Convert between a standard VCF and a variant track
#'
#' @param object For `as_variant_track()`, a `VariantAnnotation::CollapsedVCF`
#'   with fixed fields, a `fileformat` header entry and declared INFO and
#'   FORMAT fields.
#'   For `as_vcf()`, an in-memory track produced by `as_variant_track()`.
#' @param samples Unique sample names in the requested order. `NULL` retains
#'   all source samples. Use `character()` to return a sites-only VCF.
#'
#' @return `as_variant_track()` returns a [VariantTrack] with a POS-based
#'   analysis table and its source VCF. `as_vcf()` returns a `CollapsedVCF`
#'   with the current track row order and the requested sample order.
#' @details
#' These optional adapters require VariantAnnotation. The analysis table uses
#' one-based point coordinates while the source retains its ranges, alleles,
#' header, typed INFO and genotype data. Conversion projects the standard
#' writer's field representations into the analysis table using a temporary
#' file. A track retains the full input VCF and the analysis values needed to
#' verify a return conversion. Use `VariantAnnotation::ScanVcfParam()` to
#' read a region first when the input is large. `VariantTrack` sorts the input
#' by chromosome, POS and display ID.
#'
#' [retrieve_vcf()] with `as = "VariantTrack"`, row filtering and reordering
#' preserve the internal row keys. Keys are independent of variant IDs and
#' positions. Keep the keys and source metadata intact when using
#' `saveRDS()` and `readRDS()`. Return conversion requires unedited VCF fields
#' and one source record per row. Use `samples` to select or reorder samples.
#' Native tracks, merged tracks and edits to projected VCF fields cannot be
#' converted back to typed VCFs by this adapter. Sample names must be unique
#' and distinct from the track's reserved column names.
#'
#' Empty ALT elements are rejected because the standard reader can use the
#' same representation for an omitted allele (`.`) and a spanning deletion
#' (`*`). Conversion preserves the data represented in the supplied object.
#'
#' Pass the result to `VariantAnnotation::writeVcf()` for genotype-aware file
#' export. Its output reflects the standard object's fields and the writer's
#' formatting. The writer can update `fileDate`, round floating-point FORMAT
#' values, pad missing array entries or reorder fields. Use RDS for exact
#' typed values and keep the original VCF for its exact file representation.
#' [write_vcf()] and [write_variant_track()] still write eight site columns.
#' @seealso [read_vcf()], [retrieve_vcf()], [write_vcf()]
#' @examples
#' if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'   file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#'   vcf <- VariantAnnotation::readVcf(file, genome = NA_character_)
#'   track <- as_variant_track(vcf)
#'   region <- retrieve_vcf(track, chrom = "chr1", start = 12340000L, end = 12346000L,
#'     as = "VariantTrack", verbose = FALSE)
#'   selected <- as_vcf(region, samples = c("S02", "S01"))
#'   VariantAnnotation::geno(selected)$GT
#' }
#' @export
as_variant_track <- function(object) {
  require_variant_annotation()
  stop_if_not(methods::is(object, "CollapsedVCF"), "`object` must be a CollapsedVCF.")
  header <- VariantAnnotation::header(object)
  stop_if_not(methods::is(header, "VCFHeader"), "The CollapsedVCF must have a VCF header.")
  fileformat <- VariantAnnotation::meta(header)$fileformat
  stop_if_not(NROW(fileformat) == 1L && length(fileformat) == 1L &&
    isTRUE(grepl("^VCFv[0-9]+[.][0-9]+$", fileformat[[1L]])),
    "The VCF header must contain a valid fileformat declaration.")
  fixed_fields <- c("REF", "ALT", "QUAL", "FILTER")
  stop_if_not(all(fixed_fields %in% names(VariantAnnotation::fixed(object))),
    "The CollapsedVCF must contain REF, ALT, QUAL and FILTER.")
  stop_if_not(all(names(VariantAnnotation::info(object)) %in% rownames(VariantAnnotation::info(header))) &&
    all(names(VariantAnnotation::geno(object)) %in% rownames(VariantAnnotation::geno(header))),
    "The VCF header must declare every projected INFO and FORMAT field.")
  methods::validObject(object)
  alleles <- as.character(unlist(VariantAnnotation::alt(object), use.names = FALSE))
  stop_if_not(!anyNA(alleles) && all(nzchar(alleles)),
    "The CollapsedVCF contains an empty or missing ALT allele element.")

  fields <- c("CHROM", "POS", "ID", fixed_fields, "INFO")
  sample_names <- as.character(colnames(object))
  reserved <- c(fields, "#CHROM", "FORMAT", "chrom", "pos", "start", "end", "variant_id",
    "ref", "alt", "qual", "filter", "info", "variant_type", "track_source", vcf_row_key)
  stop_if_not(!any(sample_names %in% reserved), "A VCF sample name matches a reserved track column.")
  stop_if_not(length(sample_names) == ncol(object) && !anyNA(sample_names) &&
    all(nzchar(sample_names)) && !anyDuplicated(sample_names) &&
    !any(grepl("[[:space:]]", sample_names)), "The CollapsedVCF must have unique, nonempty sample names.")
  stop_if_not(all(sample_names %in% VariantAnnotation::samples(header)),
    "The VCF header must declare every sample column.")
  stop_if_not(!length(sample_names) || length(VariantAnnotation::geno(object)) > 0L,
    "A CollapsedVCF with samples must contain FORMAT data.")

  ranges <- GenomicRanges::granges(object)
  stop_if_not(all(IRanges::start(ranges) >= 1L), "The VCF positions must be positive integers.")
  object <- vcf_with_samples(object, sample_names)
  row_keys <- if (nrow(object)) paste0(basename(tempfile("gtr_vcf_")), "_", seq_len(nrow(object))) else character()
  table <- vcf_analysis_table(object, fields, sample_names, row_keys, ranges)
  track <- vcf_table_to_variant_track(table, row_keys = row_keys)
  track$meta$vcf_context <- list(source = object, row_keys = row_keys,
    projection = data.table::copy(track$data))
  track
}

#' @rdname as_variant_track
#' @export
as_vcf <- function(object, samples = NULL) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack.")
  context <- object$meta$vcf_context
  stop_if_not(is.list(context) && !is.null(context$source) && !is.null(object$data),
    "The VariantTrack has no compatible VCF provenance. Use as_variant_track() first.")
  require_variant_annotation()
  stop_if_not(methods::is(context$source, "CollapsedVCF"), "The VCF provenance is invalid.")

  table <- object$data
  keys <- table[[vcf_row_key]]
  source_rows <- match(keys, context$row_keys)
  projection_rows <- match(keys, context$projection[[vcf_row_key]])
  stop_if_not(is.character(keys) && !anyDuplicated(keys) && !anyNA(source_rows) &&
    !anyNA(projection_rows) && !anyDuplicated(names(table)),
    "The VCF row provenance is missing, duplicated or from a different source.")
  for (field in setdiff(names(context$projection), vcf_row_key)) {
    stop_if_not(identical(table[[field]], context$projection[[field]][projection_rows]),
      paste0("The projected VCF field was removed or edited: ", field, "."))
  }

  source <- context$source
  sample_names <- as.character(colnames(source))
  stop_if_not(identical(sample_names, object$meta$sample_names), "The VCF sample metadata was edited.")
  if (is.null(samples)) samples <- sample_names
  stop_if_not(is.character(samples) && !anyNA(samples) && !anyDuplicated(samples),
    "`samples` must contain unique sample names.")
  sample_cols <- match(samples, sample_names)
  stop_if_not(!anyNA(sample_cols), "Unknown VCF sample requested.")
  vcf_with_samples(source[source_rows, sample_cols], samples)
}
