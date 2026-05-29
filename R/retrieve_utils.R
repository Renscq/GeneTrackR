# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.20
# Function: Internal helpers for retrieve APIs
# Input: data.table-like genomic records
# Output: Filtered records

filter_by_region_internal <- function(dt, chrom = NULL, start = NULL, end = NULL,
                                      start_col = "start", end_col = "end",
                                      mode = c("overlap", "within")) {
  mode <- match.arg(mode)
  dt <- data.table::as.data.table(dt)
  query_chrom <- as.character(chrom)
  query_chrom <- query_chrom[!is.na(query_chrom) & nzchar(query_chrom)]
  if (length(query_chrom) > 0L && "chrom" %in% names(dt)) {
    keep <- as.character(dt[["chrom"]]) %in% query_chrom
    dt <- dt[keep]
  }
  if (!is.null(start) && !is.null(end) && all(c(start_col, end_col) %in% names(dt))) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    if (mode == "within") {
      dt <- dt[dt[[start_col]] >= s & dt[[end_col]] <= e]
    } else {
      dt <- dt[dt[[start_col]] <= e & dt[[end_col]] >= s]
    }
  }
  dt
}

match_pattern_internal <- function(dt, pattern = NULL, fields = NULL,
                                   ignore_case = TRUE, fixed = FALSE) {
  dt <- data.table::as.data.table(dt)
  if (is.null(pattern)) return(dt)
  stop_if_not(length(pattern) == 1L && !is.na(pattern) && nzchar(pattern),
              "`pattern` must be a non-empty character string.")
  fields <- intersect(fields, names(dt))
  if (length(fields) == 0L || nrow(dt) == 0L) return(dt[0])
  hit <- rep(FALSE, nrow(dt))
  for (field in fields) {
    hit <- hit | grepl(pattern, as.character(dt[[field]]), ignore.case = ignore_case, fixed = fixed)
  }
  dt[hit]
}

