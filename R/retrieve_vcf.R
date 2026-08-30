# Author: Rensc
# Date: 2026-08-31
# Version: dev004
# Function: Retrieve variants from VariantTrack objects or indexed VCF files
# Input: VariantTrack object or VCF path and filters
# Output: data.table, VariantTrack, or GRanges object

#' Retrieve variants from a VariantTrack or an indexed VCF file
#'
#' @description
#' Retrieves VCF records by optional genomic region, gene/transcript locator,
#' variant ID, type, or text pattern. With an in-memory `VariantTrack`, all
#' location arguments may be omitted to query the full object. For large
#' bgzip-compressed VCF files with a tabix index, a complete genomic region is
#' queried with `Rsamtools::scanTabix()`; non-regional filters require reading
#' the full file before filtering.
#'
#' @param object A VariantTrack object or a path to a VCF/VCF.GZ file.
#' @param pattern Optional text pattern matched against variant ID, REF, ALT, INFO, and variant type.
#' @param chrom Optional chromosome name or names. May be used alone as a chromosome filter. Required when `start` and `end` are supplied.
#' @param start Optional 1-based region start. Must be supplied together with `end`.
#' @param end Optional 1-based region end. Must be supplied together with `start`.
#' @param annotation Optional GenePred/Feature annotation object used for gene/transcript-aware retrieval.
#' @param gene_id Optional gene ID. When supplied, `annotation` is used to resolve the gene range.
#' @param transcript_id Optional transcript ID. When supplied, `annotation` is used to resolve the transcript range.
#' @param upstream Upstream flanking length in bp for gene/transcript queries.
#' @param downstream Downstream flanking length in bp for gene/transcript queries.
#' @param strand_aware Logical. Whether upstream/downstream should follow gene/transcript strand direction.
#' @param variant_id Optional variant ID vector.
#' @param variant_type Optional variant type vector, such as `SNP`, `INS`, `DEL`, or `MNV`.
#' @param keep_genotype Logical. Whether to keep FORMAT and sample genotype columns when `object` is a VCF file path.
#' @param ignore_case Logical. Whether pattern matching ignores case.
#' @param fixed Logical. Whether pattern is matched as a fixed string.
#' @param as Output type: `data.table`, `VariantTrack`, or `GRanges`.
#' @param verbose Logical. Whether to print progress messages.
#' @param progress Logical. Whether to print a compact stage-level progress indicator.
#' @return A data.table, VariantTrack, or `GenomicRanges::GRanges` object.
#' @examples
#' vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
#' retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200)
#' retrieve_vcf(vcf, variant_type = "SNP")
#' retrieve_vcf(vcf, pattern = "varA", fixed = TRUE)
#' vt <- retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200, as = "VariantTrack")
#' vt
#' retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200, as = "GRanges")
#' @export
retrieve_vcf <- function(object,
                         pattern = NULL,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         annotation = NULL,
                         gene_id = NULL,
                         transcript_id = NULL,
                         upstream = 0L,
                         downstream = 0L,
                         strand_aware = TRUE,
                         variant_id = NULL,
                         variant_type = NULL,
                         keep_genotype = TRUE,
                         ignore_case = TRUE,
                         fixed = FALSE,
                         as = c("data.table", "VariantTrack", "GRanges"),
                         verbose = TRUE,
                         progress = interactive() && isTRUE(verbose)) {
  as <- match.arg(as)
  verbose <- isTRUE(verbose)

  location <- resolve_optional_vcf_location(
    annotation = annotation,
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    start = start,
    end = end,
    upstream = upstream,
    downstream = downstream,
    strand_aware = strand_aware
  )
  chrom <- location$chrom
  start <- location$start
  end <- location$end

  if (is.character(object) && length(object) == 1L && file.exists(object)) {
    vt <- retrieve_vcf_file(
      object,
      chrom = chrom,
      start = start,
      end = end,
      keep_genotype = keep_genotype,
      verbose = verbose,
      progress = progress
    )
  } else {
    stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object or a VCF file path.")
    if (is_lazy_variant_track(object)) {
      vt <- retrieve_vcf_file(
        file = object$meta$source_file,
        chrom = chrom,
        start = start,
        end = end,
        keep_genotype = keep_genotype %||% object$meta$keep_genotype %||% TRUE,
        verbose = verbose,
        progress = progress
      )
    } else {
      vt <- object
    }
  }

  dt <- data.table::copy(vt$data)

  query_chrom <- as.character(chrom)
  query_chrom <- query_chrom[!is.na(query_chrom) & nzchar(query_chrom)]
  if (length(query_chrom) > 0L) {
    keep <- as.character(dt[["chrom"]]) %in% query_chrom
    dt <- dt[keep]
  }

  if (!is.null(start) && !is.null(end)) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    dt <- dt[dt[["pos"]] >= s & dt[["pos"]] <= e]
  }

  if (!is.null(variant_id)) {
    ids <- as.character(variant_id)
    dt <- dt[dt[["variant_id"]] %in% ids]
  }

  if (!is.null(variant_type)) {
    types <- as.character(variant_type)
    dt <- dt[dt[["variant_type"]] %in% types]
  }

  dt <- match_pattern_internal(
    dt,
    pattern = pattern,
    fields = c("variant_id", "ref", "alt", "info", "variant_type"),
    ignore_case = ignore_case,
    fixed = fixed
  )

  if (nrow(dt) > 0L) {
    data.table::setorderv(dt, intersect(c("chrom", "pos", "variant_id"), names(dt)))
  }

  if (verbose && !(is.character(object) && length(object) == 1L && file.exists(object))) {
    message("[GeneTrackR] Retrieved variants: ", format(nrow(dt), big.mark = ","), ".")
  }

  if (as == "data.table") return(dt[])
  result <- VariantTrack(dt, meta = vt$meta)
  if (as == "GRanges") return(as_granges(result))
  result
}


resolve_optional_vcf_location <- function(annotation = NULL,
                                          gene_id = NULL,
                                          transcript_id = NULL,
                                          chrom = NULL,
                                          start = NULL,
                                          end = NULL,
                                          upstream = 0L,
                                          downstream = 0L,
                                          strand_aware = TRUE) {
  has_gene <- !is.null(gene_id)
  has_tx <- !is.null(transcript_id)
  has_direct <- !is.null(chrom) || !is.null(start) || !is.null(end)

  if (has_gene || has_tx) {
    stop_if_not(!has_direct, "Do not mix `gene_id`/`transcript_id` with direct `chrom`/`start`/`end` region arguments.")
    return(resolve_haplotype_gene_region(
      annotation = annotation,
      gene_id = gene_id,
      transcript_id = transcript_id,
      upstream = upstream,
      downstream = downstream,
      strand_aware = strand_aware
    ))
  }

  upstream <- as.integer(upstream)[1L]
  downstream <- as.integer(downstream)[1L]
  stop_if_not(
    !is.na(upstream) && !is.na(downstream) && upstream >= 0L && downstream >= 0L,
    "`upstream` and `downstream` must be non-negative integers."
  )
  stop_if_not(
    upstream == 0L && downstream == 0L,
    "`upstream` and `downstream` require `gene_id` or `transcript_id`."
  )

  if (!has_direct) {
    return(list(
      locator = "all",
      chrom = NULL,
      start = NULL,
      end = NULL
    ))
  }

  stop_if_not(!is.null(chrom), "`chrom` is required when `start` or `end` is supplied.")
  has_start <- !is.null(start)
  has_end <- !is.null(end)
  stop_if_not(
    identical(has_start, has_end),
    "Supply both `start` and `end`, or neither, for a direct chromosome query."
  )

  chrom <- as.character(chrom)
  chrom <- chrom[!is.na(chrom) & nzchar(chrom)]
  stop_if_not(length(chrom) > 0L, "`chrom` must contain at least one non-empty chromosome name.")

  if (!has_start) {
    return(list(
      locator = "chrom",
      chrom = chrom,
      start = NULL,
      end = NULL
    ))
  }

  stop_if_not(length(chrom) == 1L, "A coordinate range can be queried for exactly one chromosome at a time.")
  check_region(chrom, start, end)
  list(
    locator = "region",
    chrom = chrom[1L],
    start = as.integer(start)[1L],
    end = as.integer(end)[1L]
  )
}

retrieve_vcf_file <- function(file,
                              chrom = NULL,
                              start = NULL,
                              end = NULL,
                              keep_genotype = TRUE,
                              verbose = TRUE,
                              progress = interactive() && isTRUE(verbose)) {
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  is_indexed <- has_vcf_tabix_index(file)
  verbose <- isTRUE(verbose)
  progress_msg <- vcf_progress_message(verbose && isTRUE(progress))

  if (isTRUE(is_indexed) && !is.null(chrom) && !is.null(start) && !is.null(end)) {
    stop_if_not(requireNamespace("Rsamtools", quietly = TRUE), "Package `Rsamtools` is required for indexed VCF queries.")
    region <- paste0(as.character(chrom)[1L], ":", as.integer(start)[1L], "-", as.integer(end)[1L])
    if (verbose) {
      message("[GeneTrackR] Querying indexed VCF region: ", region)
    }

    progress_msg(1L, 4L, "Reading VCF header.")
    header <- read_vcf_header_line(file)
    col_names <- parse_vcf_header_names(header)

    progress_msg(2L, 4L, "Scanning tabix index.")
    query_gr <- GenomicRanges::GRanges(
      seqnames = as.character(chrom)[1L],
      ranges = IRanges::IRanges(
        start = as.integer(start)[1L],
        end = as.integer(end)[1L]
      )
    )

    # Older Rsamtools versions do not dispatch scanTabix(TabixFile, character).
    # Using a GRanges param is compatible with TabixFile across Bioconductor
    # versions and avoids the "signature file = TabixFile, param = character" error.
    tbx <- Rsamtools::TabixFile(file)
    on.exit(try(close(tbx), silent = TRUE), add = TRUE)
    lines <- Rsamtools::scanTabix(tbx, param = query_gr)[[1L]]

    progress_msg(3L, 4L, paste0("Parsing ", format(length(lines), big.mark = ","), " VCF records."))
    if (length(lines) == 0L) {
      dt <- data.table::data.table(matrix(ncol = length(col_names), nrow = 0L))
      data.table::setnames(dt, col_names)
    } else {
      dt <- data.table::fread(
        text = paste(c(paste(col_names, collapse = "\t"), lines), collapse = "\n"),
        sep = "\t",
        header = TRUE,
        data.table = TRUE,
        showProgress = FALSE
      )
    }

    vt <- vcf_table_to_variant_track(dt, source_file = file, keep_genotype = keep_genotype)
    progress_msg(4L, 4L, "Finished indexed VCF query.")
    if (verbose) {
      sample_n <- length(vt$meta$sample_names %||% character())
      message("[GeneTrackR] Finished. variants: ", format(nrow(vt$data), big.mark = ","), "; samples: ", sample_n, ".")
    }
    return(vt)
  }

  if (verbose && !isTRUE(is_indexed)) {
    message("[GeneTrackR] No tabix index detected. Reading full VCF file: ", file)
  } else if (verbose) {
    message("[GeneTrackR] Indexed VCF detected, but no complete region was supplied. Reading full VCF file: ", file)
  }
  read_vcf(file, keep_genotype = keep_genotype, mode = "memory", verbose = verbose, progress = progress)
}

vcf_table_to_variant_track <- function(dt, source_file = NULL, keep_genotype = TRUE) {
  if (nrow(dt) == 0L && ncol(dt) == 0L) {
    standard_empty <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    dt <- data.table::data.table(matrix(ncol = length(standard_empty), nrow = 0L))
    data.table::setnames(dt, standard_empty)
  }

  if ("#CHROM" %in% names(dt) && !"CHROM" %in% names(dt)) {
    data.table::setnames(dt, "#CHROM", "CHROM")
  }
  standard_in <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  stop_if_not(all(standard_in %in% names(dt)), "A VCF table must contain standard VCF columns.")
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

  VariantTrack(
    dt,
    meta = list(
      source_file = source_file,
      format = "VCF",
      coordinate_input = "1-based position",
      coordinate_internal = "1-based position",
      indexed = if (!is.null(source_file)) has_vcf_tabix_index(source_file) else FALSE,
      has_genotype = length(sample_cols) > 0L,
      has_format = isTRUE(has_format),
      sample_names = sample_cols
    )
  )
}
