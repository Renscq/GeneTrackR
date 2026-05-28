# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.4
# Function: R wrappers for the local libBigWig C++ backend
# Input: Local bigWig file path and genomic regions
# Output: Chromosome metadata and signal intervals

normalize_bigwig_path <- function(file) {
  stop_if_not(length(file) == 1L, "`file` must be a single path.")
  file <- as.character(file)

  if (grepl("^(https?|ftp)://", file, ignore.case = TRUE)) {
    stop(
      "Remote bigWig URLs are not supported by the local libBigWig backend. ",
      "Download the file first or compile libBigWig with curl support.",
      call. = FALSE
    )
  }

  if (grepl("^file://", file, ignore.case = TRUE)) {
    stop(
      "The local libBigWig backend expects a normal file path, not a file:// URI. ",
      "Use a path like E:/path/file.bigwig.",
      call. = FALSE
    )
  }

  normalizePath(file, winslash = "/", mustWork = FALSE)
}

bw_seqinfo_cpp <- function(file) {
  out <- .Call(`_GeneTrackR_bw_seqinfo_cpp`, normalize_bigwig_path(file))
  data.table::as.data.table(out)
}

bw_query_cpp <- function(file, chrom, start, end) {
  out <- .Call(
    `_GeneTrackR_bw_query_cpp`,
    normalize_bigwig_path(file),
    as.character(chrom),
    as.integer(start),
    as.integer(end)
  )
  data.table::as.data.table(out)
}

query_bigwig_cpp <- function(file, sample_id, chrom, start, end, strand = "*") {
  dt <- bw_query_cpp(file, chrom, start, end)
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
    sample_id = sample_id,
    strand = strand
  )]
  data.table::setcolorder(dt, c("sample_id", "chrom", "start", "end", "value", "strand"))
  dt[]
}

read_bigwig_whole_cpp <- function(file, sample_id, strand = "*") {
  si <- bw_seqinfo_cpp(file)
  if (nrow(si) == 0L) {
    return(data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    ))
  }

  out <- lapply(seq_len(nrow(si)), function(i) {
    query_bigwig_cpp(
      file = file,
      sample_id = sample_id,
      chrom = si$chrom[i],
      start = 1L,
      end = si$length[i],
      strand = strand
    )
  })

  data.table::rbindlist(out, fill = TRUE)
}
