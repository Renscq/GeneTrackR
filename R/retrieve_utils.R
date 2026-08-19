# Author: Rensc
# Date: 2026-05-28
# Version: dev001
# Function: Internal helpers for retrieve APIs
# Input: data.table-like genomic records
# Output: Filtered records

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

