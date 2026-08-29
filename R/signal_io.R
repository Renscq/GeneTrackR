# Author: Rensc
# Date: 2026-08-30
# Version: dev001
# Function: Unified native R I/O helpers for genomic signal track formats
# Input: Local bedGraph, WIG, and BigWig files or canonical signal tables
# Output: Canonical 1-based closed signal tables and signal track files

normalize_signal_input_file <- function(file) {
  stop_if_not(
    is.character(file) && length(file) == 1L && !is.na(file) && nzchar(file),
    "`file` must be a single non-missing character path."
  )
  stop_if_not(!grepl("^(https?|ftp)://", file, ignore.case = TRUE), "Remote signal URLs are not supported.")
  stop_if_not(!grepl("^file://", file, ignore.case = TRUE), "Use a local signal file path instead of a file:// URI.")
  stop_if_not(file.exists(file), paste0("Signal file does not exist: ", file))
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

open_signal_text_connection <- function(file) {
  file <- normalize_signal_input_file(file)
  if (grepl("\\.(gz|bgz)$", file, ignore.case = TRUE)) {
    return(gzfile(file, open = "rt"))
  }
  base::file(file, open = "rt")
}

parse_bedgraph_lines_native <- function(lines) {
  lines <- trimws(as.character(lines))
  keep <- nzchar(lines) & !grepl("^(#|browser|track)", lines, ignore.case = TRUE)
  lines <- lines[keep]
  if (length(lines) == 0L) {
    return(data.table::data.table(
      chrom = character(), start = integer(), end = integer(), value = numeric()
    ))
  }

  fields <- strsplit(lines, "[[:space:]]+")
  stop_if_not(
    all(lengths(fields) >= 4L),
    "bedGraph contains a non-empty record with fewer than four fields."
  )

  chrom <- vapply(fields, `[[`, character(1L), 1L)
  start0 <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 2L)))
  end0 <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 3L)))
  value <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 4L)))

  stop_if_not(
    !anyNA(start0) && !anyNA(end0) && !anyNA(value),
    "bedGraph contains non-numeric start, end, or value fields."
  )
  stop_if_not(
    all(start0 >= 0 & end0 > start0),
    "bedGraph contains invalid 0-based half-open intervals."
  )
  stop_if_not(
    all(start0 + 1 <= .Machine$integer.max) && all(end0 <= .Machine$integer.max),
    "bedGraph coordinates exceed the supported integer range."
  )

  data.table::data.table(
    chrom = chrom,
    start = as.integer(start0 + 1),
    end = as.integer(end0),
    value = as.numeric(value)
  )
}

read_bedgraph_native <- function(file) {
  con <- open_signal_text_connection(file)
  on.exit(close(con), add = TRUE)

  out <- list()
  out_n <- 0L
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    chunk <- parse_bedgraph_lines_native(lines)
    if (nrow(chunk) == 0L) {
      next
    }
    out_n <- out_n + 1L
    out[[out_n]] <- chunk
  }

  if (out_n == 0L) {
    return(data.table::data.table(
      chrom = character(), start = integer(), end = integer(), value = numeric()
    ))
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("chrom", "start", "end"))
  ans[]
}

parse_wig_header_native <- function(line) {
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  mode <- tolower(fields[1L])
  args <- fields[-1L]
  parsed <- list()
  for (arg in args) {
    pair <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(pair) == 2L) {
      parsed[[tolower(pair[1L])]] <- pair[2L]
    }
  }
  list(mode = mode, args = parsed)
}

read_wig_native <- function(file) {
  con <- open_signal_text_connection(file)
  on.exit(close(con), add = TRUE)

  chrom <- NULL
  mode <- NULL
  span <- 1L
  step <- NA_integer_
  next_start <- NA_integer_
  out <- vector("list", 1024L)
  out_n <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      break
    }
    line <- trimws(line)
    if (!nzchar(line) || grepl("^(#|browser|track)", line, ignore.case = TRUE)) {
      next
    }

    if (grepl("^(fixedStep|variableStep)\\b", line, ignore.case = TRUE)) {
      header <- parse_wig_header_native(line)
      mode <- header$mode
      chrom <- header$args$chrom
      stop_if_not(!is.null(chrom) && nzchar(chrom), "WIG step header is missing `chrom`.")
      span <- suppressWarnings(as.integer(header$args$span %||% "1"))
      stop_if_not(!is.na(span) && span >= 1L, "WIG `span` must be a positive integer.")
      if (identical(mode, "fixedstep")) {
        next_start <- suppressWarnings(as.integer(header$args$start %||% NA_character_))
        step <- suppressWarnings(as.integer(header$args$step %||% "1"))
        stop_if_not(
          !is.na(next_start) && next_start >= 1L && !is.na(step) && step >= 1L,
          "fixedStep WIG headers require positive `start` and `step` values."
        )
      }
      next
    }

    fields <- strsplit(line, "[[:space:]]+")[[1L]]
    if (length(fields) >= 4L) {
      record <- parse_bedgraph_lines_native(line)
      stop_if_not(nrow(record) == 1L, "Invalid bedGraph-style record in WIG input.")
      chrom_value <- record$chrom[1L]
      start_value <- record$start[1L]
      end_value <- record$end[1L]
      signal_value <- record$value[1L]
    } else if (identical(mode, "variablestep")) {
      stop_if_not(length(fields) >= 2L, "variableStep WIG records require position and value.")
      pos <- suppressWarnings(as.numeric(fields[1L]))
      signal_value <- suppressWarnings(as.numeric(fields[2L]))
      stop_if_not(!is.na(pos) && pos >= 1 && !is.na(signal_value), "Invalid variableStep WIG record.")
      chrom_value <- chrom
      start_value <- pos
      end_value <- pos + span - 1L
    } else if (identical(mode, "fixedstep")) {
      signal_value <- suppressWarnings(as.numeric(fields[1L]))
      stop_if_not(!is.na(signal_value), "Invalid fixedStep WIG value.")
      chrom_value <- chrom
      start_value <- next_start
      end_value <- next_start + span - 1L
      next_start <- next_start + step
    } else {
      stop("WIG data were encountered before a fixedStep or variableStep header.", call. = FALSE)
    }

    stop_if_not(
      start_value <= .Machine$integer.max && end_value <= .Machine$integer.max,
      "WIG coordinates exceed the supported integer range."
    )
    out_n <- out_n + 1L
    if (out_n > length(out)) {
      length(out) <- length(out) * 2L
    }
    out[[out_n]] <- list(
      chrom = as.character(chrom_value),
      start = as.integer(start_value),
      end = as.integer(end_value),
      value = as.numeric(signal_value)
    )
  }

  if (out_n == 0L) {
    return(data.table::data.table(
      chrom = character(), start = integer(), end = integer(), value = numeric()
    ))
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)])
  data.table::setorderv(ans, c("chrom", "start", "end"))
  ans[]
}

attach_signal_sample <- function(data, sample_id, strand = "*") {
  data <- data.table::copy(data.table::as.data.table(data))
  if (nrow(data) == 0L) {
    return(empty_signal_dt())
  }
  stop_if_not(
    all(c("chrom", "start", "end", "value") %in% names(data)),
    "Signal data must contain chrom, start, end, and value columns."
  )
  data[, `:=`(
    sample_id = as.character(sample_id)[1L],
    strand = as.character(strand)[1L]
  )]
  data.table::setcolorder(data, c("sample_id", "chrom", "start", "end", "value", "strand"))
  data[]
}

read_signal_file_memory <- function(file, format, sample_id, strand = "*") {
  format <- match.arg(tolower(as.character(format)[1L]), c("bedgraph", "wig", "bigwig"))
  if (format == "bigwig") {
    return(read_bigwig_whole_native(file, sample_id, strand))
  }
  data <- if (format == "bedgraph") {
    read_bedgraph_native(file)
  } else {
    read_wig_native(file)
  }
  attach_signal_sample(data, sample_id = sample_id, strand = strand)
}

subset_signal_region <- function(data, chrom, start, end) {
  data <- data.table::as.data.table(data)
  if (nrow(data) == 0L) {
    return(empty_signal_dt())
  }
  chrom_value <- as.character(chrom)[1L]
  start_value <- as.integer(start)[1L]
  end_value <- as.integer(end)[1L]
  out <- data[
    data[["chrom"]] == chrom_value &
      data[["start"]] <= end_value &
      data[["end"]] >= start_value
  ]
  if (nrow(out) == 0L) {
    return(empty_signal_dt())
  }
  out[, .(
    sample_id,
    chrom,
    start = pmax(start, start_value),
    end = pmin(end, end_value),
    value,
    strand
  )]
}

query_signal_file_region <- function(sample_row, chrom, start, end) {
  sample_row <- data.table::as.data.table(sample_row)
  stop_if_not(nrow(sample_row) == 1L, "`sample_row` must contain exactly one sample.")

  format <- tolower(as.character(sample_row$format[1L]))
  file <- sample_row$file[1L]
  sample_id <- sample_row$sample_id[1L]
  strand <- sample_row$strand[1L] %||% "*"

  if (identical(format, "bigwig")) {
    return(query_bigwig_native(file, sample_id, chrom, start, end, strand))
  }

  if (identical(format, "bedgraph") && isTRUE(sample_row$use_tabix[1L])) {
    return(query_bedgraph_tabix(
      file = file,
      sample_id = sample_id,
      chrom = chrom,
      start = start,
      end = end,
      strand = strand,
      backend = sample_row$tabix_backend[1L],
      empty_fallback = isTRUE(sample_row$tabix_empty_fallback[1L])
    ))
  }

  data <- read_signal_file_memory(
    file = file,
    format = format,
    sample_id = sample_id,
    strand = strand
  )
  subset_signal_region(data, chrom = chrom, start = start, end = end)
}

build_signal_seqinfo <- function(files, formats, sample_ids, data = NULL) {
  seqinfo_list <- lapply(seq_along(files), function(i) {
    if (!identical(formats[i], "bigwig")) {
      return(NULL)
    }
    si <- bigwig_seqinfo_native(files[i])
    if (nrow(si) == 0L) {
      return(NULL)
    }
    si[, sample_id := sample_ids[i]]
    data.table::setcolorder(si, c("sample_id", "chrom", "length"))
    si[]
  })

  if (!is.null(data) && nrow(data) > 0L) {
    text_sample_ids <- sample_ids[formats %in% c("bedgraph", "wig")]
    if (length(text_sample_ids) > 0L) {
      text_seqinfo <- unique(
        data[data[["sample_id"]] %in% text_sample_ids, .(sample_id, chrom)]
      )
      if (nrow(text_seqinfo) > 0L) {
        text_seqinfo[, length := NA_integer_]
        seqinfo_list[[length(seqinfo_list) + 1L]] <- text_seqinfo
      }
    }
  }

  seqinfo <- data.table::rbindlist(seqinfo_list, fill = TRUE)
  if (is.null(seqinfo) || nrow(seqinfo) == 0L) {
    return(NULL)
  }
  data.table::setorderv(seqinfo, c("sample_id", "chrom"))
  seqinfo[]
}

signal_output_path <- function(outdir, sample_id, format, compress = FALSE) {
  format <- match.arg(format, c("bedgraph", "wig", "bigwig"))
  ext <- switch(format, bedgraph = ".bedgraph", wig = ".wig", bigwig = ".bigwig")
  compression_ext <- if (isTRUE(compress) && format %in% c("bedgraph", "wig")) ".gz" else ""
  file.path(outdir, paste0(sample_id, ext, compression_ext))
}

write_bedgraph_native <- function(data, file, compress = FALSE) {
  x <- data[, .(
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
    compress = if (isTRUE(compress)) "gzip" else "auto"
  )
  invisible(file)
}

write_wig_native <- function(data, file, compress = FALSE) {
  con <- if (isTRUE(compress)) gzfile(file, open = "wt") else base::file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  x <- data.table::copy(data)
  x[, span := as.integer(end) - as.integer(start) + 1L]
  stop_if_not(
    all(!is.na(x$span) & x$span >= 1L),
    "WIG output requires valid intervals with end >= start."
  )
  data.table::setorderv(x, c("chrom", "start", "end"))
  split_x <- split(x, x$chrom)
  for (chrom_value in names(split_x)) {
    y <- split_x[[chrom_value]]
    block_id <- cumsum(c(TRUE, y$span[-1L] != y$span[-nrow(y)]))
    split_y <- split(y, block_id)
    for (block in split_y) {
      span_value <- unique(block$span)
      writeLines(paste0("variableStep chrom=", chrom_value, " span=", span_value), con)
      writeLines(
        paste(as.integer(block$start), as.numeric(block$value), sep = "\t"),
        con
      )
    }
  }
  invisible(file)
}

write_signal_file_memory <- function(data,
                                     file,
                                     format,
                                     compress = FALSE,
                                     chrom_sizes = NULL,
                                     sample_id = NULL) {
  format <- match.arg(format, c("bedgraph", "wig", "bigwig"))
  if (format == "bedgraph") {
    return(write_bedgraph_native(data, file, compress = compress))
  }
  if (format == "wig") {
    return(write_wig_native(data, file, compress = compress))
  }
  stop_if_not(!is.null(chrom_sizes), "`chrom_sizes` is required for bigWig output.")
  prepared <- prepare_bigwig_signal(data, chrom_sizes, sample_id = sample_id)
  write_bigwig_native(file, prepared, chrom_sizes)
}

# Compatibility wrappers retained for internal callers during the 0.7.x migration.
read_bedgraph_file <- function(file, sample_id, strand = "*") {
  read_signal_file_memory(file, "bedgraph", sample_id, strand)
}

read_wig_file <- function(file, sample_id, strand = "*") {
  read_signal_file_memory(file, "wig", sample_id, strand)
}

write_bedgraph_table <- function(dt, file, compress = FALSE) {
  write_bedgraph_native(dt, file, compress = compress)
}

write_wig_table <- function(dt, file, compress = FALSE) {
  write_wig_native(dt, file, compress = compress)
}
