# Author: Rensc
# Date: 2026-08-30
# Version: dev007
# Function: Write BwgTrack objects to bedGraph, wig, or bigWig files
# Input: BwgTrack object
# Output: Signal track files

#' Write BwgTrack signal data
#'
#' @description
#' Writes a `BwgTrack` object to bedGraph, wig, or bigWig files. All three
#' formats are written directly in R. bigWig output uses GeneTrackR's native R
#' binary writer and requires chromosome sizes; no external conversion program
#' or compiled BigWig backend is used by `write_bwg()`.
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
#' @return Invisibly returns a data.table containing sample IDs and output files.
#' @examples
#' \dontrun{
#' rnaseq <- read_bwg(
#'   system.file(
#'     "extdata",
#'     c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
#'     package = "GeneTrackR"
#'   ),
#'   sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#' write_bwg(rnaseq, outdir = tempdir(), format = "bedgraph", overwrite = TRUE)
#' }
#' @export
write_bwg <- function(object,
                      outdir,
                      format = c("bedgraph", "wig", "bigwig"),
                      samples = NULL,
                      chrom_sizes = NULL,
                      overwrite = FALSE,
                      compress = FALSE) {
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
    chrom_sizes_table <- prepare_chrom_sizes_table(chrom_sizes)
  }

  out <- vector("list", nrow(sample_tbl))
  for (i in seq_len(nrow(sample_tbl))) {
    sid <- sample_tbl$sample_id[i]
    x <- dt[sample_id == sid]
    if (nrow(x) == 0L) next

    if (format == "bedgraph") {
      file <- file.path(
        outdir,
        paste0(sid, ".bedgraph", if (compress) ".gz" else "")
      )
      check_output_file(file, overwrite)
      write_bedgraph_table(x, file, compress = compress)
    } else if (format == "wig") {
      file <- file.path(
        outdir,
        paste0(sid, ".wig", if (compress) ".gz" else "")
      )
      check_output_file(file, overwrite)
      write_wig_table(x, file, compress = compress)
    } else {
      file <- file.path(outdir, paste0(sid, ".bigwig"))
      check_output_file(file, overwrite)
      x <- prepare_bigwig_signal(x, chrom_sizes_table, sample_id = sid)
      write_bigwig_native(file, x, chrom_sizes_table)
    }

    out[[i]] <- data.table::data.table(
      sample_id = sid,
      file = normalizePath(file, winslash = "/", mustWork = FALSE),
      format = format
    )
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
      "Sample '", sid, "' is ", src_format, " but requested ", requested_format,
      ". Use retrieve_bwg(as = \"BwgTrack\") to create an in-memory object first."
    ))
    stop_if_not(
      file.exists(src),
      paste0("Source signal file does not exist: ", src)
    )

    compression_ext <- if (
      format %in% c("bedgraph", "wig") &&
        grepl("\\.bgz$", src, ignore.case = TRUE)
    ) {
      ".bgz"
    } else if (
      format %in% c("bedgraph", "wig") &&
        grepl("\\.gz$", src, ignore.case = TRUE)
    ) {
      ".gz"
    } else {
      ""
    }
    ext <- switch(
      format,
      bedgraph = ".bedgraph",
      wig = ".wig",
      bigwig = ".bigwig"
    )
    dst <- file.path(outdir, paste0(sid, ext, compression_ext))
    check_output_file(dst, overwrite)
    copied <- file.copy(src, dst, overwrite = overwrite)
    stop_if_not(
      isTRUE(copied),
      paste0("Failed to copy signal file for sample '", sid, "': ", src)
    )
    out[[i]] <- data.table::data.table(
      sample_id = sid,
      file = normalizePath(dst, winslash = "/", mustWork = TRUE),
      format = format
    )
  }
  invisible(data.table::rbindlist(out, fill = TRUE))
}

write_bedgraph_table <- function(dt, file, compress = FALSE) {
  x <- dt[, .(
    chrom = as.character(chrom),
    start = as.integer(start) - 1L,
    end = as.integer(end),
    value = as.numeric(value)
  )]
  data.table::fwrite(
    x,
    file,
    sep = "\t",
    col.names = FALSE,
    compress = if (compress) "gzip" else "auto"
  )
}

write_wig_table <- function(dt, file, compress = FALSE) {
  con <- if (compress) gzfile(file, open = "wt") else file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  x <- data.table::copy(dt)
  x[, "span" := as.integer(end) - as.integer(start) + 1L]
  stop_if_not(
    all(!is.na(x$span) & x$span >= 1L),
    "WIG output requires valid intervals with end >= start."
  )
  data.table::setorderv(x, c("chrom", "start", "end"))
  split_x <- split(x, x$chrom)
  for (chr in names(split_x)) {
    y <- split_x[[chr]]
    block_id <- cumsum(c(TRUE, y$span[-1L] != y$span[-nrow(y)]))
    split_y <- split(y, block_id)
    for (block in split_y) {
      span <- unique(block$span)
      writeLines(paste0("variableStep chrom=", chr, " span=", span), con)
      writeLines(
        paste(as.integer(block$start), as.numeric(block$value), sep = "\t"),
        con
      )
    }
  }
}

prepare_chrom_sizes_table <- function(chrom_sizes) {
  if (is.character(chrom_sizes) && length(chrom_sizes) == 1L) {
    stop_if_not(
      file.exists(chrom_sizes),
      paste0("Chromosome size file does not exist: ", chrom_sizes)
    )
    cs <- data.table::fread(
      chrom_sizes,
      header = FALSE,
      select = 1:2,
      col.names = c("chrom", "size"),
      showProgress = FALSE
    )
  } else {
    cs <- data.table::as.data.table(chrom_sizes)
    stop_if_not(
      ncol(cs) >= 2L,
      "`chrom_sizes` must have at least two columns: chrom and size."
    )
    cs <- cs[, .(
      chrom = as.character(cs[[1L]]),
      size = as.numeric(cs[[2L]])
    )]
  }

  cs[, `:=`(chrom = as.character(chrom), size = as.numeric(size))]
  stop_if_not(
    nrow(cs) > 0L,
    "`chrom_sizes` must contain at least one chromosome."
  )
  stop_if_not(
    all(!is.na(cs$chrom) & nzchar(cs$chrom)),
    "`chrom_sizes` contains missing or empty chromosome names."
  )
  stop_if_not(
    !anyDuplicated(cs$chrom),
    "`chrom_sizes` contains duplicated chromosome names."
  )
  stop_if_not(
    all(is.finite(cs$size) & cs$size >= 1 & cs$size <= 4294967295),
    "Chromosome sizes must be finite integers between 1 and 4,294,967,295."
  )
  stop_if_not(
    all(cs$size == floor(cs$size)),
    "Chromosome sizes must be integer values."
  )
  cs[]
}

prepare_bigwig_signal <- function(dt, chrom_sizes, sample_id = NULL) {
  x <- data.table::copy(dt)
  required <- c("chrom", "start", "end", "value")
  stop_if_not(
    all(required %in% names(x)),
    "Signal data must contain chrom, start, end, and value columns."
  )

  x[, `:=`(
    chrom = as.character(chrom),
    start = as.numeric(start),
    end = as.numeric(end),
    value = as.numeric(value)
  )]
  label <- if (is.null(sample_id)) "Signal" else paste0("Sample '", sample_id, "'")

  stop_if_not(
    all(!is.na(x$chrom) & nzchar(x$chrom)),
    paste0(label, " contains missing chromosome names.")
  )
  stop_if_not(
    all(is.finite(x$start) & is.finite(x$end)),
    paste0(label, " contains non-finite coordinates.")
  )
  stop_if_not(
    all(x$start == floor(x$start) & x$end == floor(x$end)),
    paste0(label, " contains non-integer coordinates.")
  )
  stop_if_not(
    all(x$start >= 1 & x$end >= x$start),
    paste0(label, " contains invalid 1-based closed intervals.")
  )
  stop_if_not(
    all(is.finite(x$value)),
    paste0(label, " contains non-finite signal values.")
  )
  stop_if_not(
    all(x$chrom %in% chrom_sizes$chrom),
    paste0(label, " contains chromosomes absent from `chrom_sizes`.")
  )

  x[, chrom_order__ := match(chrom, chrom_sizes$chrom)]
  x[, chrom_size__ := chrom_sizes$size[chrom_order__]]
  stop_if_not(
    all(x$end <= x$chrom_size__),
    paste0(label, " contains intervals extending beyond chromosome sizes.")
  )
  data.table::setorderv(x, c("chrom_order__", "start", "end"))

  x[, previous_end__ := data.table::shift(end), by = chrom]
  has_overlap <- any(!is.na(x$previous_end__) & x$start <= x$previous_end__)
  stop_if_not(
    !has_overlap,
    paste0(
      label,
      " contains overlapping intervals; bigWig output requires non-overlapping intervals within each chromosome."
    )
  )

  x[, c("chrom_order__", "chrom_size__", "previous_end__") := NULL]
  x[]
}

check_output_file <- function(file, overwrite = FALSE) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop(paste0("File exists: ", file), call. = FALSE)
  }
  invisible(TRUE)
}
