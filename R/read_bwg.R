# Author: Rensc
# Date: 2026-08-30
# Version: dev007
# Function: Read signal tracks through the unified native R I/O layer
# Input: Signal track file paths
# Output: Schema-v2 BwgTrack object

#' Read bedGraph, wig, or bigWig signal tracks
#'
#' @param files Signal file paths.
#' @param format Input format. Use auto, bedgraph, bigwig, or wig.
#' @param sample_names Optional sample names. If NULL, names are inferred from file basenames after removing compression suffixes and signal-track suffixes such as `.bedgraph.gz`, `.bigwig`, `.bw`, or `.wig`.
#' @param strand Optional strand labels for files.
#' @param mode lazy stores file paths, memory reads data into R.
#' @param genome Optional genome label.
#' @param check_chrom Whether to check chromosome names for in-memory files.
#' @param use_tabix Whether to use tabix/Rsamtools for indexed bedGraph-like files in lazy mode. Accepts `"auto"`, `"yes"`, `"no"`, TRUE, or FALSE. `TRUE` is equivalent to `"yes"`, and `FALSE` is equivalent to `"no"`.
#' @param verbose Logical. Whether to print read/setup progress messages.
#' @param tabix_empty_fallback Logical. Whether an empty tabix result should be verified by a full-file fread query. Default FALSE for performance. Set TRUE only when you suspect the tabix index coordinate convention is incompatible with the queried region.
#' @details
#' `format = "auto"` infers the input type from file extensions after removing
#' compression suffixes such as `.gz` or `.bgz`. `sample_names = NULL` infers
#' names from file basenames, removing suffixes such as `.bedgraph.gz`, `.bw`,
#' `.bigwig`, and `.wig`.
#'
#' For bigWig and wig files, strand information is not stored in the file, so
#' the sample is treated as unstranded. For paired plus/minus bedGraph files,
#' provide `strand = c("+", "-")` and matching `sample_names` if strand-specific
#' filtering is required.
#'
#' `use_tabix = "auto"` uses indexed querying only when a `.tbi` index and an available backend are detected. GeneTrackR first checks the system `tabix` command and then the R package `Rsamtools`. Otherwise, lazy bedGraph queries fall back to full-file reading.
#'
#' BigWig metadata, full-memory loading, and lazy regional queries use the
#' GeneTrackR native R BigWig reader. The signal subsystem does not require
#' Rcpp, a compiled BigWig library, or an external BigWig executable.
#'
#' `BwgTrack$seqinfo` follows the schema-v2 sequence metadata contract. bigWig
#' inputs record chromosome lengths from the file header. In-memory bedGraph and
#' wig inputs record observed chromosome names with unknown lengths (`NA`).
#'
#' bedGraph, WIG, and BigWig memory reads are dispatched through one internal
#' native R I/O layer so all formats return the same canonical signal schema.
#' @return A BwgTrack object.
#' @examples
#' \dontrun{
#' rnaseq_files <- system.file(
#'   "extdata",
#'   c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
#'   package = "GeneTrackR"
#' )
#' riboseq_files <- system.file(
#'   "extdata",
#'   c("gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"),
#'   package = "GeneTrackR"
#' )
#'
#' rnaseq <- read_bwg(
#'   rnaseq_files,
#'   format = "bedgraph",
#'   sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#' riboseq <- read_bwg(
#'   riboseq_files,
#'   format = "bedgraph",
#'   sample_names = c("Ribo_seq_plus", "Ribo_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#' }
#' @export
read_bwg <- function(files,
                     format = c("auto", "bedgraph", "bigwig", "wig"),
                     sample_names = NULL,
                     strand = NULL,
                     mode = c("lazy", "memory"),
                     genome = NULL,
                     check_chrom = TRUE,
                     use_tabix = c("auto", "yes", "no"),
                     tabix_empty_fallback = FALSE,
                     verbose = TRUE) {
  format <- match.arg(format)
  mode <- match.arg(mode)
  use_tabix <- normalize_use_tabix_arg(use_tabix)
  tabix_empty_fallback <- isTRUE(tabix_empty_fallback)
  verbose <- isTRUE(verbose)
  files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  stop_if_not(all(file.exists(files)), "One or more signal files do not exist.")

  if (verbose) {
    message(sprintf("[GeneTrackR] Preparing %s signal file(s).", length(files)))
  }
  sample_names <- safe_sample_names(files, sample_names)
  formats <- infer_signal_format(files, format)
  if (is.null(strand)) {
    strand <- rep("*", length(files))
  }
  stop_if_not(length(strand) == length(files), "`strand` must have the same length as `files`.")

  strand <- normalize_sample_strand(strand, formats)
  has_strand <- infer_sample_has_strand(strand, formats)

  has_tabix <- vapply(files, has_tabix_index, logical(1))
  tabix_backend <- rep(NA_character_, length(files))
  available_backend <- get_tabix_backend()

  if (use_tabix == "no") {
    tabix_enabled <- rep(FALSE, length(files))
  } else {
    tabix_enabled <- formats == "bedgraph" & has_tabix & !is.na(available_backend)
    tabix_backend[tabix_enabled] <- available_backend
  }

  if (use_tabix == "yes" && any(formats == "bedgraph" & !tabix_enabled)) {
    missing_index <- formats == "bedgraph" & !has_tabix
    missing_backend <- formats == "bedgraph" & has_tabix & is.na(available_backend)
    msg <- c()
    if (any(missing_index)) {
      msg <- c(msg, sprintf("%s file(s) do not have a .tbi index", sum(missing_index)))
    }
    if (any(missing_backend)) {
      msg <- c(msg, "no tabix backend was found on this system. Install htslib/tabix or the R package Rsamtools")
    }
    warning(
      paste0("Some bedGraph files cannot use indexed tabix queries: ", paste(msg, collapse = "; "), ". Falling back to full-file fread for those files."),
      call. = FALSE
    )
  }

  samples <- data.table::data.table(
    sample_id = sample_names,
    file = files,
    format = formats,
    strand = strand,
    has_strand = has_strand,
    has_tabix = has_tabix,
    use_tabix = tabix_enabled,
    tabix_backend = tabix_backend,
    tabix_empty_fallback = tabix_empty_fallback,
    library_size = NA_real_,
    norm_method = "none",
    scale_factor = 1
  )

  data <- NULL
  if (mode == "memory") {
    if (verbose) {
      message("[GeneTrackR] Loading signal files into memory.")
    }
    if (any(formats == "bigwig")) {
      warning(
        "Full-memory bigWig loading may use a large amount of memory. ",
        "For routine analysis, use `mode = 'lazy'` and retrieve_bwg().",
        call. = FALSE
      )
    }
    data <- data.table::rbindlist(lapply(seq_along(files), function(i) {
      if (verbose) {
        message(sprintf("[GeneTrackR] Loading %s/%s: %s", i, length(files), sample_names[i]))
      }
      read_signal_file_memory(
        file = files[i],
        format = formats[i],
        sample_id = sample_names[i],
        strand = strand[i]
      )
    }), fill = TRUE)
  }

  seqinfo <- build_signal_seqinfo(
    files = files,
    formats = formats,
    sample_ids = sample_names,
    data = data
  )

  meta <- list(
    mode = mode,
    genome = genome,
    coordinate = "1-based closed",
    backend = "GeneTrackR-native-R"
  )

  if (verbose) {
    message(sprintf("[GeneTrackR] Finished reading signal metadata. samples: %s; mode: %s.", nrow(samples), mode))
  }
  obj <- BwgTrack(
    samples = samples,
    data = data,
    meta = meta,
    validation = make_empty_validation(),
    seqinfo = seqinfo
  )
  if (check_chrom) {
    obj$validation <- validate_bwg(obj)
  }
  obj
}


normalize_use_tabix_arg <- function(use_tabix) {
  if (is.logical(use_tabix)) {
    stop_if_not(length(use_tabix) == 1L && !is.na(use_tabix), "`use_tabix` must be TRUE, FALSE, 'auto', 'yes', or 'no'.")
    return(if (isTRUE(use_tabix)) "yes" else "no")
  }
  use_tabix <- as.character(use_tabix)[1L]
  if (is.na(use_tabix) || !nzchar(use_tabix)) {
    stop("`use_tabix` must be TRUE, FALSE, 'auto', 'yes', or 'no'.", call. = FALSE)
  }
  use_tabix <- tolower(use_tabix)
  use_tabix <- switch(
    use_tabix,
    "true" = "yes",
    "false" = "no",
    "t" = "yes",
    "f" = "no",
    "1" = "yes",
    "0" = "no",
    use_tabix
  )
  match.arg(use_tabix, choices = c("auto", "yes", "no"))
}

normalize_sample_strand <- function(strand, formats) {
  strand <- as.character(strand)
  strand[is.na(strand) | strand == ""] <- "*"
  valid <- strand %in% c("+", "-", "*", ".", "both", "ignore")
  stop_if_not(all(valid), "`strand` must be one of '+', '-', '*', '.', 'both', or 'ignore'.")
  strand[strand %in% c(".", "both", "ignore")] <- "*"

  # bigWig and wig do not contain strand information in the file itself.
  # A user-provided strand can describe the sample/source file, but it must
  # not be interpreted as per-record strand information for filtering.
  strand[formats %in% c("bigwig", "wig")] <- "*"
  strand
}

infer_sample_has_strand <- function(strand, formats) {
  # Standard bigWig/wig files are unstranded signal tracks. BedGraph is also
  # formally unstranded, but paired plus/minus bedGraph files are common, so a
  # user-provided '+' or '-' is respected for bedGraph only.
  formats == "bedgraph" & strand %in% c("+", "-")
}

#' Validate a BwgTrack object
#'
#' @param object A BwgTrack object.
#' @return A validation list.
#' @export
validate_bwg <- function(object) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  invalid <- data.table::data.table(sample_id = character(), reason = character())
  warnings <- character()

  if (anyDuplicated(object$samples$sample_id)) {
    invalid <- rbind(invalid, data.table::data.table(sample_id = object$samples$sample_id[duplicated(object$samples$sample_id)], reason = "duplicated sample_id"))
  }
  if (!is.null(object$data)) {
    bad <- object$data[is.na(chrom) | is.na(start) | is.na(end) | start > end | is.na(value), unique(sample_id)]
    if (length(bad) > 0L) {
      invalid <- rbind(invalid, data.table::data.table(sample_id = bad, reason = "invalid signal record"))
    }
  }

  list(
    invalid_records = invalid,
    invalid_summary = invalid[, .N, by = reason],
    warnings = warnings
  )
}
