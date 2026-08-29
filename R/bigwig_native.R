# Author: Rensc
# Date: 2026-08-29
# Version: dev002
# Function: Read local BigWig files using native R binary I/O
# Input: Local BigWig file paths and genomic intervals
# Output: Chromosome metadata and 1-based closed signal intervals
# Source: Adapted from the bwTools 0.8.10 native BigWig reader

.GTR_BIGWIG_MAGIC <- 0x888FFC26
.GTR_CIRTREE_MAGIC <- 0x78CA8C91
.GTR_BIGWIG_INDEX_MAGIC <- 0x2468ACE0
.GTR_BIGWIG_MAX_EXACT_UINT64 <- 2^53 - 1
.GTR_BIGWIG_HEADER_SIZE <- 64L
.GTR_BIGWIG_INDEX_HEADER_SIZE <- 48L
.GTR_BIGWIG_DATA_HEADER_SIZE <- 24L
.GTR_BIGWIG_SUMMARY_SIZE <- 40L

.gtr_bigwig_cache <- new.env(parent = emptyenv())

.gtr_bw_stop <- function(message) {
  stop(as.character(message)[1L], call. = FALSE)
}

.gtr_bw_validate_local_file <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    .gtr_bw_stop("`file` must be a single non-missing character string.")
  }
  if (grepl("^(https?|ftp)://", file, ignore.case = TRUE)) {
    .gtr_bw_stop("Remote BigWig files are not supported by the native R backend.")
  }
  if (grepl("^file://", file, ignore.case = TRUE)) {
    .gtr_bw_stop("Use a local file path instead of a file:// URI.")
  }
  if (!file.exists(file)) {
    .gtr_bw_stop(paste0("BigWig file does not exist: ", file))
  }
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

.gtr_bw_empty_signal <- function() {
  data.table::data.table(
    chrom = character(),
    start = integer(),
    end = integer(),
    value = numeric()
  )
}

.gtr_bw_cache_key <- function(file) {
  info <- file.info(file)
  paste(file, as.numeric(info$size[1L]), as.numeric(info$mtime[1L]), sep = "|")
}

.gtr_bw_read_exact <- function(con, offset, n) {
  offset <- as.numeric(offset)[1L]
  n <- as.numeric(n)[1L]

  if (
    !is.finite(offset) || offset < 0 ||
      offset > .GTR_BIGWIG_MAX_EXACT_UINT64
  ) {
    .gtr_bw_stop("Invalid or unsupported BigWig file offset.")
  }
  if (!is.finite(n) || n < 0 || n > .Machine$integer.max) {
    .gtr_bw_stop("Invalid or unsupported BigWig read size.")
  }

  seek(con, where = offset, origin = "start")
  out <- readBin(con, what = "raw", n = as.integer(n), size = 1L)
  if (length(out) != as.integer(n)) {
    .gtr_bw_stop(
      paste0(
        "Unexpected end of BigWig file at offset ",
        format(offset, scientific = FALSE),
        "."
      )
    )
  }
  out
}

.gtr_bw_uint <- function(x, offset, size, endian = "little") {
  offset <- as.integer(offset)[1L]
  size <- as.integer(size)[1L]
  idx <- seq.int(offset + 1L, offset + size)

  if (offset < 0L || size < 1L || max(idx) > length(x)) {
    .gtr_bw_stop("Binary integer read exceeds the available buffer.")
  }

  bytes <- as.numeric(x[idx])
  if (identical(endian, "big")) {
    bytes <- rev(bytes)
  }
  sum(bytes * 256^(seq_along(bytes) - 1L))
}

.gtr_bw_u8 <- function(x, offset) {
  as.integer(.gtr_bw_uint(x, offset, 1L, "little"))
}

.gtr_bw_u16 <- function(x, offset, endian) {
  as.integer(.gtr_bw_uint(x, offset, 2L, endian))
}

.gtr_bw_u32 <- function(x, offset, endian) {
  .gtr_bw_uint(x, offset, 4L, endian)
}

.gtr_bw_u64 <- function(x, offset, endian) {
  out <- .gtr_bw_uint(x, offset, 8L, endian)
  if (!is.finite(out) || out > .GTR_BIGWIG_MAX_EXACT_UINT64) {
    .gtr_bw_stop(
      paste0(
        "A 64-bit BigWig offset exceeds the exact integer range ",
        "supported by native R numeric values."
      )
    )
  }
  out
}

.gtr_bw_uint_vector <- function(x, offsets, size, endian = "little") {
  offsets <- as.numeric(offsets)
  size <- as.integer(size)[1L]

  if (length(offsets) == 0L) {
    return(numeric())
  }
  if (size < 1L || any(!is.finite(offsets)) || any(offsets < 0)) {
    .gtr_bw_stop("Invalid vectorized binary integer read.")
  }
  if (any(offsets + size > length(x))) {
    .gtr_bw_stop("Vectorized binary integer read exceeds the available buffer.")
  }

  idx <- rep(offsets, each = size) +
    rep(seq.int(0L, size - 1L), times = length(offsets)) + 1
  bytes <- matrix(
    as.numeric(x[idx]),
    nrow = length(offsets),
    ncol = size,
    byrow = TRUE
  )
  if (identical(endian, "big")) {
    bytes <- bytes[, seq.int(size, 1L), drop = FALSE]
  }
  drop(bytes %*% 256^(seq_len(size) - 1L))
}

.gtr_bw_u32_vector <- function(x, offsets, endian) {
  .gtr_bw_uint_vector(x, offsets, 4L, endian)
}

.gtr_bw_float32_vector <- function(x, offsets, endian) {
  bits <- .gtr_bw_u32_vector(x, offsets, endian)
  sign_bit <- floor(bits / 2^31)
  exponent <- floor((bits %% 2^31) / 2^23)
  fraction <- bits %% 2^23
  sign_value <- ifelse(sign_bit == 0, 1, -1)

  out <- sign_value * (1 + fraction / 2^23) * 2^(exponent - 127)
  subnormal <- exponent == 0
  out[subnormal] <- sign_value[subnormal] *
    (fraction[subnormal] / 2^23) * 2^-126

  infinite <- exponent == 255 & fraction == 0
  nan_value <- exponent == 255 & fraction != 0
  out[infinite] <- sign_value[infinite] * Inf
  out[nan_value] <- NaN
  out
}

.gtr_bw_double <- function(x, offset, endian) {
  offset <- as.integer(offset)[1L]
  idx <- seq.int(offset + 1L, offset + 8L)
  if (offset < 0L || max(idx) > length(x)) {
    .gtr_bw_stop("Binary double read exceeds the available buffer.")
  }

  con <- rawConnection(x[idx], open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = numeric(), n = 1L, size = 8L, endian = endian)
}

.gtr_bw_fixed_string <- function(x) {
  zero <- which(as.integer(x) == 0L)
  if (length(zero) > 0L) {
    end <- zero[1L] - 1L
    x <- if (end > 0L) x[seq_len(end)] else raw()
  }
  if (length(x) == 0L) {
    return("")
  }
  rawToChar(x)
}

.gtr_bw_detect_endian <- function(header_raw) {
  if (.gtr_bw_u32(header_raw, 0L, "little") == .GTR_BIGWIG_MAGIC) {
    return("little")
  }
  if (.gtr_bw_u32(header_raw, 0L, "big") == .GTR_BIGWIG_MAGIC) {
    return("big")
  }
  .gtr_bw_stop("The input file is not a valid BigWig file.")
}

.gtr_bw_read_header <- function(con) {
  x <- .gtr_bw_read_exact(con, 0, .GTR_BIGWIG_HEADER_SIZE)
  endian <- .gtr_bw_detect_endian(x)

  list(
    endian = endian,
    version = .gtr_bw_u16(x, 4L, endian),
    n_levels = .gtr_bw_u16(x, 6L, endian),
    chrom_tree_offset = .gtr_bw_u64(x, 8L, endian),
    data_offset = .gtr_bw_u64(x, 16L, endian),
    index_offset = .gtr_bw_u64(x, 24L, endian),
    field_count = .gtr_bw_u16(x, 32L, endian),
    defined_field_count = .gtr_bw_u16(x, 34L, endian),
    sql_offset = .gtr_bw_u64(x, 36L, endian),
    summary_offset = .gtr_bw_u64(x, 44L, endian),
    uncompress_buf_size = .gtr_bw_u32(x, 52L, endian),
    extension_offset = .gtr_bw_u64(x, 56L, endian)
  )
}

.gtr_bw_read_total_summary <- function(con, header) {
  if (header$summary_offset <= 0) {
    return(NULL)
  }

  x <- .gtr_bw_read_exact(
    con,
    header$summary_offset,
    .GTR_BIGWIG_SUMMARY_SIZE
  )
  list(
    n_bases_covered = .gtr_bw_u64(x, 0L, header$endian),
    min_value = .gtr_bw_double(x, 8L, header$endian),
    max_value = .gtr_bw_double(x, 16L, header$endian),
    sum_data = .gtr_bw_double(x, 24L, header$endian),
    sum_squared = .gtr_bw_double(x, 32L, header$endian)
  )
}

.gtr_bw_read_chrom_node <- function(
    con,
    offset,
    key_size,
    endian,
    chrom,
    chrom_length,
    depth = 0L) {
  if (depth > 128L) {
    .gtr_bw_stop("Chromosome B+ tree exceeds the supported recursion depth.")
  }

  node_header <- .gtr_bw_read_exact(con, offset, 4L)
  is_leaf <- .gtr_bw_u8(node_header, 0L) != 0L
  n_children <- .gtr_bw_u16(node_header, 2L, endian)
  entry_size <- key_size + 8L

  if (n_children == 0L) {
    return(list(chrom = chrom, length = chrom_length))
  }

  entries <- .gtr_bw_read_exact(
    con,
    offset + 4,
    as.numeric(n_children) * entry_size
  )
  for (i in seq_len(n_children)) {
    base <- (i - 1L) * entry_size
    key_raw <- entries[seq.int(base + 1L, base + key_size)]

    if (is_leaf) {
      tid <- .gtr_bw_u32(entries, base + key_size, endian)
      seq_len_value <- .gtr_bw_u32(
        entries,
        base + key_size + 4L,
        endian
      )
      r_index <- as.integer(tid + 1)
      if (r_index < 1L || r_index > length(chrom)) {
        .gtr_bw_stop("Invalid chromosome ID in BigWig chromosome tree.")
      }
      chrom[r_index] <- .gtr_bw_fixed_string(key_raw)
      chrom_length[r_index] <- seq_len_value
    } else {
      child_offset <- .gtr_bw_u64(entries, base + key_size, endian)
      child <- .gtr_bw_read_chrom_node(
        con,
        child_offset,
        key_size,
        endian,
        chrom,
        chrom_length,
        depth + 1L
      )
      chrom <- child$chrom
      chrom_length <- child$length
    }
  }

  list(chrom = chrom, length = chrom_length)
}

.gtr_bw_read_chrom_tree <- function(con, header) {
  x <- .gtr_bw_read_exact(con, header$chrom_tree_offset, 32L)
  endian <- header$endian

  if (.gtr_bw_u32(x, 0L, endian) != .GTR_CIRTREE_MAGIC) {
    .gtr_bw_stop("Invalid BigWig chromosome B+ tree magic number.")
  }

  key_size <- .gtr_bw_u32(x, 8L, endian)
  value_size <- .gtr_bw_u32(x, 12L, endian)
  item_count <- .gtr_bw_u64(x, 16L, endian)

  if (key_size < 1 || key_size > .Machine$integer.max) {
    .gtr_bw_stop("Invalid chromosome key size in BigWig file.")
  }
  if (value_size != 8) {
    .gtr_bw_stop("Unsupported chromosome value size in BigWig file.")
  }
  if (item_count > .Machine$integer.max) {
    .gtr_bw_stop("Unsupported number of chromosomes in BigWig file.")
  }

  n <- as.integer(item_count)
  out <- .gtr_bw_read_chrom_node(
    con = con,
    offset = header$chrom_tree_offset + 32,
    key_size = as.integer(key_size),
    endian = endian,
    chrom = rep(NA_character_, n),
    chrom_length = rep(NA_real_, n)
  )

  if (anyNA(out$chrom) || any(!nzchar(out$chrom))) {
    .gtr_bw_stop(
      "Failed to resolve all chromosomes from the BigWig chromosome tree."
    )
  }

  data.table::data.table(
    tid = seq_len(n) - 1L,
    chrom = out$chrom,
    length = out$length
  )
}

# Internal metadata reader for the native R BigWig backend.
bigwig_metadata_native <- function(file, use_cache = TRUE) {
  file <- .gtr_bw_validate_local_file(file)
  key <- .gtr_bw_cache_key(file)

  if (
    isTRUE(use_cache) &&
      exists(key, envir = .gtr_bigwig_cache, inherits = FALSE)
  ) {
    return(get(key, envir = .gtr_bigwig_cache, inherits = FALSE))
  }

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)

  header <- .gtr_bw_read_header(con)
  total_summary <- .gtr_bw_read_total_summary(con, header)
  chromosomes <- .gtr_bw_read_chrom_tree(con, header)

  out <- list(
    file = file,
    header = header,
    total_summary = total_summary,
    chromosomes = chromosomes
  )

  if (isTRUE(use_cache)) {
    assign(key, out, envir = .gtr_bigwig_cache)
  }
  out
}

.gtr_bw_read_index_header <- function(
    con,
    header,
    index_offset = header$index_offset) {
  index_offset <- as.numeric(index_offset)[1L]
  if (!is.finite(index_offset) || index_offset <= 0) {
    .gtr_bw_stop("The BigWig file does not contain the requested R-tree index.")
  }

  x <- .gtr_bw_read_exact(
    con,
    index_offset,
    .GTR_BIGWIG_INDEX_HEADER_SIZE
  )
  endian <- header$endian

  if (.gtr_bw_u32(x, 0L, endian) != .GTR_BIGWIG_INDEX_MAGIC) {
    .gtr_bw_stop("Invalid BigWig R-tree index magic number.")
  }

  list(
    block_size = .gtr_bw_u32(x, 4L, endian),
    n_items = .gtr_bw_u64(x, 8L, endian),
    chrom_start = .gtr_bw_u32(x, 16L, endian),
    base_start = .gtr_bw_u32(x, 20L, endian),
    chrom_end = .gtr_bw_u32(x, 24L, endian),
    base_end = .gtr_bw_u32(x, 28L, endian),
    index_size = .gtr_bw_u64(x, 32L, endian),
    items_per_slot = .gtr_bw_u32(x, 40L, endian),
    root_offset = index_offset + .GTR_BIGWIG_INDEX_HEADER_SIZE
  )
}

.gtr_bw_index_entry_overlaps <- function(
    chrom_start,
    base_start,
    chrom_end,
    base_end,
    tid,
    query_start,
    query_end) {
  if (tid < chrom_start || tid > chrom_end) {
    return(FALSE)
  }
  if (chrom_start != chrom_end) {
    if (tid == chrom_start && base_start >= query_end) {
      return(FALSE)
    }
    if (tid == chrom_end && base_end <= query_start) {
      return(FALSE)
    }
    return(TRUE)
  }
  !(base_start >= query_end || base_end <= query_start)
}

.gtr_bw_collect_blocks <- function(
    con,
    node_offset,
    endian,
    tid,
    query_start,
    query_end,
    depth = 0L) {
  if (depth > 128L) {
    .gtr_bw_stop("BigWig R-tree exceeds the supported recursion depth.")
  }

  node_header <- .gtr_bw_read_exact(con, node_offset, 4L)
  is_leaf <- .gtr_bw_u8(node_header, 0L) != 0L
  n_children <- .gtr_bw_u16(node_header, 2L, endian)

  if (n_children == 0L) {
    return(data.table::data.table(offset = numeric(), size = numeric()))
  }

  entry_size <- if (is_leaf) 32L else 24L
  entries <- .gtr_bw_read_exact(
    con,
    node_offset + 4,
    as.numeric(n_children) * entry_size
  )

  out <- vector("list", n_children)
  out_n <- 0L

  for (i in seq_len(n_children)) {
    base <- (i - 1L) * entry_size
    chrom_start <- .gtr_bw_u32(entries, base, endian)
    base_start <- .gtr_bw_u32(entries, base + 4L, endian)
    chrom_end <- .gtr_bw_u32(entries, base + 8L, endian)
    base_end <- .gtr_bw_u32(entries, base + 12L, endian)
    data_offset <- .gtr_bw_u64(entries, base + 16L, endian)

    if (tid < chrom_start) {
      break
    }
    if (!.gtr_bw_index_entry_overlaps(
      chrom_start,
      base_start,
      chrom_end,
      base_end,
      tid,
      query_start,
      query_end
    )) {
      next
    }

    if (is_leaf) {
      data_size <- .gtr_bw_u64(entries, base + 24L, endian)
      out_n <- out_n + 1L
      out[[out_n]] <- data.table::data.table(
        offset = data_offset,
        size = data_size
      )
    } else {
      child <- .gtr_bw_collect_blocks(
        con,
        data_offset,
        endian,
        tid,
        query_start,
        query_end,
        depth + 1L
      )
      if (nrow(child) > 0L) {
        out_n <- out_n + 1L
        out[[out_n]] <- child
      }
    }
  }

  if (out_n == 0L) {
    return(data.table::data.table(offset = numeric(), size = numeric()))
  }

  ans <- data.table::rbindlist(out[seq_len(out_n)])
  ans <- unique(ans, by = c("offset", "size"))
  ans[]
}

.gtr_bw_decompress_block <- function(x) {
  tryCatch(
    base::memDecompress(x, type = "gzip"),
    error = function(e) {
      .gtr_bw_stop(
        paste0(
          "Failed to decompress a zlib-compressed BigWig data block: ",
          conditionMessage(e)
        )
      )
    }
  )
}

.gtr_bw_decode_block_vectors <- function(
    x,
    endian,
    tid,
    query_start,
    query_end) {
  if (length(x) < .GTR_BIGWIG_DATA_HEADER_SIZE) {
    .gtr_bw_stop("Truncated BigWig data block header.")
  }

  block_tid <- .gtr_bw_u32(x, 0L, endian)
  block_start <- .gtr_bw_u32(x, 4L, endian)
  block_step <- .gtr_bw_u32(x, 12L, endian)
  block_span <- .gtr_bw_u32(x, 16L, endian)
  block_type <- .gtr_bw_u8(x, 20L)
  n_items <- .gtr_bw_u16(x, 22L, endian)

  if (block_tid != tid || n_items == 0L) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }

  item_size <- switch(
    as.character(block_type),
    "1" = 12L,
    "2" = 8L,
    "3" = 4L,
    .gtr_bw_stop("Unsupported BigWig data block type.")
  )

  required_size <- .GTR_BIGWIG_DATA_HEADER_SIZE +
    as.numeric(n_items) * item_size
  if (required_size > length(x)) {
    .gtr_bw_stop("Truncated BigWig data block payload.")
  }

  item_offsets <- .GTR_BIGWIG_DATA_HEADER_SIZE +
    seq.int(0, n_items - 1L) * item_size

  if (block_type == 1L) {
    start0 <- .gtr_bw_u32_vector(x, item_offsets, endian)
    end0 <- .gtr_bw_u32_vector(x, item_offsets + 4L, endian)
    value <- .gtr_bw_float32_vector(x, item_offsets + 8L, endian)
  } else if (block_type == 2L) {
    start0 <- .gtr_bw_u32_vector(x, item_offsets, endian)
    end0 <- start0 + block_span
    value <- .gtr_bw_float32_vector(x, item_offsets + 4L, endian)
  } else {
    start0 <- block_start + seq.int(0, n_items - 1L) * block_step
    end0 <- start0 + block_span
    value <- .gtr_bw_float32_vector(x, item_offsets, endian)
  }

  keep <- end0 > query_start & start0 < query_end
  if (!any(keep)) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }

  start0 <- pmax(start0[keep], query_start)
  end0 <- pmin(end0[keep], query_end)
  value <- as.numeric(value[keep])
  valid <- end0 > start0

  if (!any(valid)) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }

  list(
    start0 = as.numeric(start0[valid]),
    end0 = as.numeric(end0[valid]),
    value = value[valid]
  )
}

# Internal chromosome metadata reader for the native R BigWig backend.
bigwig_seqinfo_native <- function(file) {
  metadata <- bigwig_metadata_native(file)
  x <- metadata$chromosomes

  lengths <- rep(NA_integer_, nrow(x))
  safe <- is.finite(x$length) & x$length <= .Machine$integer.max
  lengths[safe] <- as.integer(x$length[safe])

  data.table::data.table(
    chrom = as.character(x$chrom),
    length = lengths
  )
}

# Internal interval query for the native R BigWig backend.
bigwig_query_native <- function(file, chrom, start, end) {
  file <- .gtr_bw_validate_local_file(file)
  if (
    !is.character(chrom) || length(chrom) != 1L ||
      is.na(chrom) || !nzchar(chrom)
  ) {
    .gtr_bw_stop("`chrom` must be a single non-missing character string.")
  }

  start <- suppressWarnings(as.integer(start)[1L])
  end <- suppressWarnings(as.integer(end)[1L])
  if (is.na(start) || start < 1L) {
    .gtr_bw_stop("`start` must be >= 1 in 1-based closed coordinates.")
  }
  if (is.na(end) || end < start) {
    .gtr_bw_stop("`end` must be >= `start`.")
  }

  metadata <- bigwig_metadata_native(file)
  chrom_idx <- match(chrom, metadata$chromosomes$chrom)
  if (is.na(chrom_idx)) {
    .gtr_bw_stop(
      paste0("Chromosome `", chrom, "` was not found in the BigWig file.")
    )
  }

  chrom_length <- metadata$chromosomes$length[chrom_idx]
  query_start <- as.numeric(start - 1L)
  query_end <- min(as.numeric(end), chrom_length)

  if (query_start >= chrom_length || query_start >= query_end) {
    return(.gtr_bw_empty_signal())
  }

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)

  index <- .gtr_bw_read_index_header(con, metadata$header)
  blocks <- .gtr_bw_collect_blocks(
    con = con,
    node_offset = index$root_offset,
    endian = metadata$header$endian,
    tid = metadata$chromosomes$tid[chrom_idx],
    query_start = query_start,
    query_end = query_end
  )

  if (nrow(blocks) == 0L) {
    return(.gtr_bw_empty_signal())
  }

  start_chunks <- vector("list", nrow(blocks))
  end_chunks <- vector("list", nrow(blocks))
  value_chunks <- vector("list", nrow(blocks))
  chunk_n <- 0L
  interval_n <- 0
  needs_sort <- FALSE
  previous_start <- -Inf
  previous_end <- -Inf

  for (i in seq_len(nrow(blocks))) {
    block <- .gtr_bw_read_exact(con, blocks$offset[i], blocks$size[i])
    if (metadata$header$uncompress_buf_size > 0) {
      block <- .gtr_bw_decompress_block(block)
    }

    decoded <- .gtr_bw_decode_block_vectors(
      block,
      metadata$header$endian,
      metadata$chromosomes$tid[chrom_idx],
      query_start,
      query_end
    )
    n_decoded <- length(decoded$start0)
    if (n_decoded < 1L) {
      next
    }

    start1 <- decoded$start0 + 1
    end1 <- decoded$end0
    if (
      any(start1 > .Machine$integer.max) ||
        any(end1 > .Machine$integer.max)
    ) {
      .gtr_bw_stop(
        "Queried coordinates exceed the supported integer coordinate range."
      )
    }

    if (n_decoded > 1L) {
      start_diff <- diff(start1)
      end_diff <- diff(end1)
      if (any(start_diff < 0 | (start_diff == 0 & end_diff < 0))) {
        needs_sort <- TRUE
      }
    }
    if (
      start1[1L] < previous_start ||
        (start1[1L] == previous_start && end1[1L] < previous_end)
    ) {
      needs_sort <- TRUE
    }

    previous_start <- start1[n_decoded]
    previous_end <- end1[n_decoded]
    interval_n <- interval_n + as.numeric(n_decoded)
    if (interval_n > .Machine$integer.max) {
      .gtr_bw_stop(
        "Queried BigWig interval count exceeds the supported R vector length."
      )
    }

    chunk_n <- chunk_n + 1L
    start_chunks[[chunk_n]] <- as.integer(start1)
    end_chunks[[chunk_n]] <- as.integer(end1)
    value_chunks[[chunk_n]] <- decoded$value
  }

  if (chunk_n < 1L) {
    return(.gtr_bw_empty_signal())
  }

  start_value <- unlist(start_chunks[seq_len(chunk_n)], use.names = FALSE)
  end_value <- unlist(end_chunks[seq_len(chunk_n)], use.names = FALSE)
  signal_value <- unlist(value_chunks[seq_len(chunk_n)], use.names = FALSE)

  ans <- data.table::data.table(
    chrom = rep(chrom, length(start_value)),
    start = start_value,
    end = end_value,
    value = signal_value
  )

  if (isTRUE(needs_sort)) {
    data.table::setorderv(ans, c("start", "end"))
  }
  ans[]
}

# Internal adapter returning the standard GeneTrackR signal schema.
query_bigwig_native <- function(file,
                                sample_id,
                                chrom,
                                start,
                                end,
                                strand = "*") {
  dt <- bigwig_query_native(
    file = file,
    chrom = chrom,
    start = start,
    end = end
  )

  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    ))
  }

  dt[, `:=`(
    sample_id = as.character(sample_id)[1L],
    strand = as.character(strand)[1L]
  )]
  data.table::setcolorder(
    dt,
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  dt[]
}


# Internal adapter for explicit full-memory BigWig loading.
read_bigwig_whole_native <- function(file, sample_id, strand = "*") {
  seqinfo <- bigwig_seqinfo_native(file)
  if (nrow(seqinfo) == 0L) {
    return(data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    ))
  }

  if (anyNA(seqinfo$length)) {
    .gtr_bw_stop(
      "Full-memory BigWig loading does not support chromosome lengths above the R integer range."
    )
  }

  out <- lapply(seq_len(nrow(seqinfo)), function(i) {
    query_bigwig_native(
      file = file,
      sample_id = sample_id,
      chrom = seqinfo$chrom[i],
      start = 1L,
      end = seqinfo$length[i],
      strand = strand
    )
  })

  data.table::rbindlist(out, fill = TRUE)
}
