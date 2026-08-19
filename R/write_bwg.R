# Author: Rensc
# Date: 2026-08-19
# Version: dev001
# Function: Write BwgTrack objects to bedGraph, wig, or bigWig files
# Input: BwgTrack object
# Output: Signal track files

#' Write BwgTrack signal data
#'
#' @description
#' Writes a `BwgTrack` object to bedGraph, wig, or bigWig files. bedGraph and
#' wig output are implemented directly in R. bigWig output requires the external
#' UCSC `bedGraphToBigWig` command and chromosome sizes.
#'
#' For in-memory `BwgTrack` objects, signal records are written from `object$data`.
#' For lazy objects, direct conversion is not possible because signal records are
#' not loaded; if the requested format matches the original file format, the
#' original files are copied to the output directory.
#'
#' @param object A `BwgTrack` object.
#' @param outdir Output directory.
#' @param format Output format. One of `bedgraph`, `wig`, or `bigwig`.
#' @param samples Optional sample IDs to write. Default writes all samples.
#' @param chrom_sizes Chromosome sizes for bigWig output. Can be a file path or
#' a data.frame/data.table with two columns: chromosome and size.
#' @param overwrite Whether to overwrite existing files.
#' @param compress Whether to gzip-compress bedGraph or wig output.
#' @param bedGraphToBigWig Path to UCSC `bedGraphToBigWig`. By default the
#' command is searched from `PATH`.
#' @return Invisibly returns a data.table containing sample IDs and output files.
#' @examples
#' \dontrun{
#' write_bwg(bg, outdir = "tracks", format = "bedgraph")
#' write_bwg(bg, outdir = "tracks", format = "wig")
#' write_bwg(bg, outdir = "tracks", format = "bigwig", chrom_sizes = "chrom.sizes")
#' }
#' @export
write_bwg <- function(object,
                      outdir,
                      format = c("bedgraph", "wig", "bigwig"),
                      samples = NULL,
                      chrom_sizes = NULL,
                      overwrite = FALSE,
                      compress = FALSE,
                      bedGraphToBigWig = Sys.which("bedGraphToBigWig")) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  format <- match.arg(format)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sample_tbl <- data.table::copy(object$samples)
  if (!is.null(samples)) {
    samples <- as.character(samples)
    sample_tbl <- sample_tbl[sample_id %in% samples]
  }
  stop_if_not(nrow(sample_tbl) > 0L, "No samples were selected for writing.")

  if (is.null(object$data)) {
    return(write_bwg_lazy_copy(object, sample_tbl, outdir, format, overwrite))
  }

  dt <- data.table::copy(object$data)
  dt <- dt[sample_id %in% sample_tbl$sample_id]
  stop_if_not(nrow(dt) > 0L, "No in-memory signal records are available for selected samples.")
  data.table::setorderv(dt, c("sample_id", "chrom", "start", "end"))

  if (format == "bigwig") {
    stop_if_not(!is.null(chrom_sizes), "`chrom_sizes` is required when `format = 'bigwig'`.")
    stop_if_not(nzchar(bedGraphToBigWig), "`bedGraphToBigWig` was not found. Provide its path or add it to PATH.")
    chrom_sizes_file <- prepare_chrom_sizes_file(chrom_sizes)
  }

  out <- vector("list", nrow(sample_tbl))
  for (i in seq_len(nrow(sample_tbl))) {
    sid <- sample_tbl$sample_id[i]
    x <- dt[sample_id == sid]
    if (nrow(x) == 0L) next
    if (format == "bedgraph") {
      file <- file.path(outdir, paste0(sid, ".bedgraph", if (compress) ".gz" else ""))
      check_output_file(file, overwrite)
      write_bedgraph_table(x, file, compress = compress)
    } else if (format == "wig") {
      file <- file.path(outdir, paste0(sid, ".wig", if (compress) ".gz" else ""))
      check_output_file(file, overwrite)
      write_wig_table(x, file, compress = compress)
    } else {
      file <- file.path(outdir, paste0(sid, ".bigwig"))
      check_output_file(file, overwrite)
      tmp_bg <- tempfile(pattern = paste0(sid, "_"), fileext = ".bedgraph")
      on.exit(unlink(tmp_bg), add = TRUE)
      write_bedgraph_table(x, tmp_bg, compress = FALSE)
      cmd_status <- system2(bedGraphToBigWig, args = c(normalizePath(tmp_bg, winslash = "/"), normalizePath(chrom_sizes_file, winslash = "/"), normalizePath(file, winslash = "/", mustWork = FALSE)))
      stop_if_not(identical(cmd_status, 0L), paste0("bedGraphToBigWig failed for sample: ", sid))
    }
    out[[i]] <- data.table::data.table(sample_id = sid, file = normalizePath(file, winslash = "/", mustWork = FALSE), format = format)
  }
  invisible(data.table::rbindlist(out, fill = TRUE))
}

write_bwg_lazy_copy <- function(object, sample_tbl, outdir, format, overwrite) {
  fmt_map <- c(bedgraph = "bedgraph", wig = "wig", bigwig = "bigwig")
  out <- vector("list", nrow(sample_tbl))
  for (i in seq_len(nrow(sample_tbl))) {
    sid <- sample_tbl$sample_id[i]
    src <- sample_tbl$file[i]
    src_format <- tolower(sample_tbl$format[i])
    requested_format <- fmt_map[[format]]
    stop_if_not(identical(src_format, requested_format), paste0(
      "Lazy BwgTrack objects can only be copied when output format matches the original format. ",
      "Sample '", sid, "' is ", src_format, " but requested ", requested_format, ". Use retrieve_bwg(as = \"BwgTrack\") to create an in-memory object first."
    ))
    stop_if_not(file.exists(src), paste0("Source signal file does not exist: ", src))

    compression_ext <- if (format %in% c("bedgraph", "wig") && grepl("\\.bgz$", src, ignore.case = TRUE)) {
      ".bgz"
    } else if (format %in% c("bedgraph", "wig") && grepl("\\.gz$", src, ignore.case = TRUE)) {
      ".gz"
    } else {
      ""
    }
    ext <- switch(format, bedgraph = ".bedgraph", wig = ".wig", bigwig = ".bigwig")
    dst <- file.path(outdir, paste0(sid, ext, compression_ext))
    check_output_file(dst, overwrite)
    copied <- file.copy(src, dst, overwrite = overwrite)
    stop_if_not(isTRUE(copied), paste0("Failed to copy signal file for sample '", sid, "': ", src))
    out[[i]] <- data.table::data.table(sample_id = sid, file = normalizePath(dst, winslash = "/", mustWork = TRUE), format = format)
  }
  invisible(data.table::rbindlist(out, fill = TRUE))
}

write_bedgraph_table <- function(dt, file, compress = FALSE) {
  x <- dt[, .(chrom = as.character(chrom), start = as.integer(start) - 1L, end = as.integer(end), value = as.numeric(value))]
  data.table::fwrite(x, file, sep = "\t", col.names = FALSE, compress = if (compress) "gzip" else "auto")
}

write_wig_table <- function(dt, file, compress = FALSE) {
  con <- if (compress) gzfile(file, open = "wt") else file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  x <- data.table::copy(dt)
  x[, "span" := as.integer(end) - as.integer(start) + 1L]
  stop_if_not(all(!is.na(x$span) & x$span >= 1L), "WIG output requires valid intervals with end >= start.")
  data.table::setorderv(x, c("chrom", "start", "end"))
  split_x <- split(x, x$chrom)
  for (chr in names(split_x)) {
    y <- split_x[[chr]]
    block_id <- cumsum(c(TRUE, y$span[-1L] != y$span[-nrow(y)]))
    split_y <- split(y, block_id)
    for (block in split_y) {
      span <- unique(block$span)
      writeLines(paste0("variableStep chrom=", chr, " span=", span), con)
      writeLines(paste(as.integer(block$start), as.numeric(block$value), sep = "\t"), con)
    }
  }
}

prepare_chrom_sizes_file <- function(chrom_sizes) {
  if (is.character(chrom_sizes) && length(chrom_sizes) == 1L) {
    stop_if_not(file.exists(chrom_sizes), paste0("Chromosome size file does not exist: ", chrom_sizes))
    return(chrom_sizes)
  }
  cs <- data.table::as.data.table(chrom_sizes)
  stop_if_not(ncol(cs) >= 2L, "`chrom_sizes` must have at least two columns: chrom and size.")
  out <- tempfile(fileext = ".chrom.sizes")
  data.table::fwrite(cs[, .(chrom = as.character(cs[[1L]]), size = as.integer(cs[[2L]]))], out, sep = "\t", col.names = FALSE)
  out
}

check_output_file <- function(file, overwrite = FALSE) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop(paste0("File exists: ", file), call. = FALSE)
  }
  invisible(TRUE)
}
