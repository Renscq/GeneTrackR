# Author: Rensc
# Date: 2026-05-31
# Version: 0.5.0
# Function: Read VCF files into genome-level VariantTrack objects
# Input: VCF files
# Output: VariantTrack objects

#' Read a VCF file as a VariantTrack
#'
#' @description
#' Reads VCF records and stores them as genomic point variants. Standard VCF
#' columns are always imported. If genotype sample columns are present, they are
#' preserved for haplotype analysis. For bgzip-compressed VCF files with a tabix
#' index, `chrom`, `start`, and `end` can be supplied to load only a genomic
#' interval instead of parsing the whole VCF file.
#'
#' @param file VCF file path. Plain VCF, gzip-compressed VCF, and bgzip VCF are supported.
#' @param keep_genotype Logical. Whether to keep FORMAT and sample genotype columns.
#' @param mode Reading mode. `memory` parses records immediately, `lazy` stores indexed VCF metadata only, and `auto` uses lazy mode for indexed VCF files when no region is supplied.
#' @param chrom Optional chromosome name for indexed regional reading.
#' @param start Optional 1-based region start for indexed regional reading.
#' @param end Optional 1-based region end for indexed regional reading.
#' @param verbose Logical. Whether to print progress messages.
#' @param progress Logical. Whether to print a compact stage-level progress indicator.
#' @return A VariantTrack object.
#' @examples
#' \dontrun{
#' variants <- read_vcf("variants.vcf.gz")
#' region_variants <- read_vcf("variants.vcf.gz", chrom = "chr1", start = 1, end = 10000)
#' plot_variant(region_variants, chrom = "chr1", start = 1, end = 10000)
#' }
#' @export
read_vcf <- function(file,
                     keep_genotype = TRUE,
                     mode = c("auto", "memory", "lazy"),
                     chrom = NULL,
                     start = NULL,
                     end = NULL,
                     verbose = TRUE,
                     progress = interactive() && isTRUE(verbose)) {
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))
  mode <- match.arg(mode)

  file_path <- normalizePath(file, winslash = "/", mustWork = FALSE)
  verbose <- isTRUE(verbose)
  progress_msg <- vcf_progress_message(verbose && isTRUE(progress))

  has_region <- !is.null(chrom) || !is.null(start) || !is.null(end)
  is_indexed <- has_vcf_tabix_index(file_path)
  if (!isTRUE(has_region) && mode %in% c("auto", "lazy") && isTRUE(is_indexed)) {
    if (verbose) {
      message("[GeneTrackR] Indexed VCF detected. Creating lazy VariantTrack: ", file_path)
    }
    return(make_lazy_vcf_track(
      file = file_path,
      keep_genotype = keep_genotype,
      verbose = verbose,
      progress = progress
    ))
  }
  if (!isTRUE(has_region) && mode == "lazy" && !isTRUE(is_indexed)) {
    warning("`mode = 'lazy'` requires a bgzip-compressed VCF with a .tbi index. Falling back to memory mode.", call. = FALSE)
  }

  if (isTRUE(has_region)) {
    stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end), "`chrom`, `start`, and `end` must all be supplied for regional VCF reading.")
    if (is_indexed) {
      if (verbose) {
        message("[GeneTrackR] Indexed VCF detected. Reading region with tabix: ", file_path)
      }
      vt <- retrieve_vcf_file(
        file = file_path,
        chrom = chrom,
        start = start,
        end = end,
        keep_genotype = keep_genotype,
        verbose = verbose,
        progress = progress
      )
      return(vt)
    }
    if (verbose) {
      message("[GeneTrackR] No tabix index detected. Falling back to full VCF reading before region filtering.")
    }
  }

  if (verbose) {
    file_size <- tryCatch(file.info(file_path)$size, error = function(e) NA_real_)
    size_label <- if (is.na(file_size)) "unknown size" else format_file_size(file_size)
    message("[GeneTrackR] Reading VCF file: ", file_path, " (", size_label, ")")
    if (has_vcf_tabix_index(file_path) && !has_region) {
      message("[GeneTrackR] Tabix index detected. For large files, use `read_vcf(file, chrom=, start=, end=)` or `retrieve_vcf(file, ...)` to avoid full-file parsing.")
    }
  }

  progress_msg(1L, 4L, "Reading VCF header.")
  header <- read_vcf_header_line(file_path)
  col_names <- parse_vcf_header_names(header)

  progress_msg(2L, 4L, "Loading VCF records.")
  dt <- data.table::fread(
    file_path,
    skip = "#CHROM",
    header = TRUE,
    sep = "\t",
    data.table = TRUE,
    showProgress = FALSE
  )

  progress_msg(3L, 4L, "Standardizing VCF table.")
  vt <- vcf_table_to_variant_track(dt, source_file = file_path, keep_genotype = keep_genotype)

  if (isTRUE(has_region)) {
    query_chrom <- as.character(chrom)[1L]
    query_start <- as.integer(start)[1L]
    query_end <- as.integer(end)[1L]
    vt$data <- vt$data[as.character(vt$data[["chrom"]]) == query_chrom & as.integer(vt$data[["pos"]]) >= query_start & as.integer(vt$data[["pos"]]) <= query_end]
  }

  progress_msg(4L, 4L, "Finished VCF reading.")
  if (verbose) {
    sample_n <- length(vt$meta$sample_names %||% character())
    message("[GeneTrackR] Finished. variants: ", format(nrow(vt$data), big.mark = ","), "; samples: ", sample_n, ".")
  }

  vt
}

#' Read a VCF file as a VariantTrack
#'
#' @description Backward-compatible alias of `read_vcf()`.
#' @inheritParams read_vcf
#' @return A VariantTrack object.
#' @export
read_vcf_track <- function(file,
                           keep_genotype = TRUE,
                           mode = c("auto", "memory", "lazy"),
                           chrom = NULL,
                           start = NULL,
                           end = NULL,
                           verbose = TRUE,
                           progress = interactive() && isTRUE(verbose)) {
  mode <- match.arg(mode)
  read_vcf(
    file = file,
    keep_genotype = keep_genotype,
    mode = mode,
    chrom = chrom,
    start = start,
    end = end,
    verbose = verbose,
    progress = progress
  )
}

read_vcf_header_line <- function(file) {
  con <- if (grepl("\\.gz$|\\.bgz$", file, ignore.case = TRUE)) gzfile(file, open = "rt") else file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    x <- readLines(con, n = 1L, warn = FALSE)
    if (length(x) == 0L) break
    if (startsWith(x, "#CHROM")) return(x)
  }
  stop("VCF header line starting with #CHROM was not found.", call. = FALSE)
}

parse_vcf_header_names <- function(header) {
  header <- sub("^#", "", header)
  strsplit(header, "\t", fixed = TRUE)[[1L]]
}

has_vcf_tabix_index <- function(file) {
  file <- as.character(file)[1L]
  if (is.na(file) || !nzchar(file)) {
    return(FALSE)
  }

  # A valid tabix index for bgzip-compressed VCF is normally
  # <file>.tbi. CSI is also accepted for very large contigs.
  # Do not infer <prefix>.tbi from <prefix>.vcf.gz because it can
  # incorrectly match unrelated stale indexes and make non-indexed
  # files look indexed.
  file.exists(paste0(file, ".tbi")) || file.exists(paste0(file, ".csi"))
}

vcf_progress_message <- function(enabled) {
  force(enabled)
  function(step, total, label = NULL) {
    if (!isTRUE(enabled)) return(invisible(NULL))
    pct <- if (total > 0L) step / total * 100 else 100
    width <- 30L
    filled <- max(0L, min(width, round(width * pct / 100)))
    bar <- paste0(strrep("=", filled), strrep(".", width - filled))
    msg <- sprintf("[GeneTrackR] Progress [%s] %s/%s (%3.0f%%)", bar, step, total, pct)
    if (!is.null(label) && nzchar(label)) msg <- paste(msg, label)
    message(msg)
    invisible(NULL)
  }
}

