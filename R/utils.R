# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.0
# Function: Internal utility functions
# Input: Generic R objects and coordinate tables
# Output: Validated and transformed helper objects

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

arg_match <- function(x, choices) {
  x <- match.arg(x, choices)
  x
}

stop_if_not <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

parse_comma_integer <- function(x) {
  if (length(x) != 1L || is.na(x) || x == "") {
    return(integer())
  }
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  vals <- vals[vals != ""]
  as.integer(vals)
}

paste_comma_integer <- function(x) {
  paste0(paste(as.integer(x), collapse = ","), ",")
}

normalize_chrom_filter <- function(chrom) {
  if (is.null(chrom)) {
    return(NULL)
  }
  as.character(chrom)
}

check_region <- function(chrom = NULL, start = NULL, end = NULL) {
  if (!is.null(start) || !is.null(end)) {
    stop_if_not(!is.null(chrom), "`chrom` is required when `start` or `end` is provided.")
    stop_if_not(!is.null(start) && !is.null(end), "Both `start` and `end` are required for region queries.")
    stop_if_not(is.numeric(start) && is.numeric(end), "`start` and `end` must be numeric.")
    stop_if_not(length(start) == 1L && length(end) == 1L, "`start` and `end` must be length-one values.")
    stop_if_not(start <= end, "`start` must be less than or equal to `end`.")
  }
  invisible(TRUE)
}

strip_signal_file_extensions <- function(x) {
  x <- basename(as.character(x))

  # Remove compression extensions first, then remove signal-track extensions.
  # This keeps biological/sample suffixes such as mRNA2 while removing
  # file-format suffixes such as .bedgraph.gz, .bigwig, .bw, or .wig.
  x <- sub("\\.(gz|bgz|bz2|xz|zip)$", "", x, ignore.case = TRUE)
  x <- sub("\\.(bedgraph|bdg|bg|bigwig|bw|wig)$", "", x, ignore.case = TRUE)
  x
}

safe_sample_names <- function(files, sample_names = NULL) {
  if (is.null(sample_names)) {
    sample_names <- strip_signal_file_extensions(files)
  }
  sample_names <- as.character(sample_names)
  stop_if_not(length(sample_names) == length(files), "`sample_names` must have the same length as `files`.")
  stop_if_not(!anyDuplicated(sample_names), "`sample_names` must be unique.")
  sample_names
}

infer_signal_format <- function(files, format) {
  format <- match.arg(format, c("auto", "bedgraph", "bigwig", "wig"))
  if (format != "auto") {
    return(rep(format, length(files)))
  }
  ext <- tolower(tools::file_ext(files))
  out <- ifelse(ext %in% c("bw", "bigwig"), "bigwig",
    ifelse(ext %in% c("wig"), "wig", "bedgraph")
  )
  out
}

make_empty_validation <- function() {
  list(
    invalid_records = data.frame(),
    invalid_summary = data.frame(),
    warnings = character()
  )
}
