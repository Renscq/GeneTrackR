# Author: Rensc
# Date: 2026-08-31
# Version: dev009
# Function: Query, normalize, summarize, slice, merge, and write signal tracks
# Input: BwgTrack objects and genomic regions
# Output: Signal tables, subset objects, and exported signal files

# Internal signal query function.
.query_bwg_internal <- function(object,
                      chrom,
                      start,
                      end,
                      samples = NULL,
                      strand = c("ignore", "+", "-", "both", "auto"),
                      strand_policy = c("ignore_unstranded", "strict"),
                      verbose = FALSE,
                      progress = interactive() && isTRUE(verbose),
                      keep_empty_samples = FALSE,
                      tabix_empty_fallback = NULL) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  strand <- match.arg(strand)
  strand_policy <- match.arg(strand_policy)
  verbose <- isTRUE(verbose)
  progress <- isTRUE(progress)
  keep_empty_samples <- isTRUE(keep_empty_samples)
  if (is.null(tabix_empty_fallback)) {
    tabix_empty_fallback <- NULL
  } else {
    tabix_empty_fallback <- isTRUE(tabix_empty_fallback)
  }
  check_region(chrom, start, end)
  chrom_value <- as.character(chrom)
  query_start <- as.integer(start)
  query_end <- as.integer(end)

  sample_tbl <- data.table::copy(object$samples)
  if (!"use_tabix" %in% names(sample_tbl)) {
    sample_tbl[, "use_tabix" := FALSE]
  }
  if (!"tabix_backend" %in% names(sample_tbl)) {
    sample_tbl[, "tabix_backend" := NA_character_]
  }
  if (!"has_strand" %in% names(sample_tbl)) {
    sample_tbl[, "has_strand" := FALSE]
  }
  if (!"tabix_empty_fallback" %in% names(sample_tbl)) {
    sample_tbl[, "tabix_empty_fallback" := FALSE]
  }
  if (!is.null(tabix_empty_fallback)) {
    sample_tbl[, "tabix_empty_fallback" := tabix_empty_fallback]
  }
  if (!is.null(samples)) {
    samples_value <- as.character(samples)
    sample_tbl <- sample_tbl[sample_tbl[["sample_id"]] %in% samples_value]
  }

  sample_tbl <- filter_sample_table_by_strand(
    sample_tbl = sample_tbl,
    strand = strand,
    strand_policy = strand_policy
  )

  if (nrow(sample_tbl) == 0L) {
    return(data.table::data.table(sample_id = character(), chrom = character(), start = integer(), end = integer(), value = numeric(), strand = character()))
  }
  expected_sample_ids <- unique(as.character(sample_tbl[["sample_id"]]))

  if (!is.null(object$data)) {
    dt <- object$data[object$data[["sample_id"]] %in% sample_tbl$sample_id & object$data[["chrom"]] == chrom_value & object$data[["start"]] <= query_end & object$data[["end"]] >= query_start]
    dt <- ensure_signal_strand_column(dt)
    dt <- dt[, .(
      sample_id, chrom,
      start = pmax(start, query_start),
      end = pmin(end, query_end),
      value,
      strand
    )]
    dt <- filter_signal_strand(dt, strand, sample_tbl, strand_policy = strand_policy)
    if (keep_empty_samples) {
      dt <- complete_empty_signal_tracks_for_query(
        dt,
        sample_ids = expected_sample_ids,
        chrom = chrom_value,
        start = query_start,
        end = query_end,
        strand = strand
      )
    }
    return(dt)
  }

  if (verbose) {
    message(sprintf(
      "[GeneTrackR] Querying %s sample(s): %s:%s-%s",
      nrow(sample_tbl), chrom_value, query_start, query_end
    ))
  }
  pb <- NULL
  if (progress && nrow(sample_tbl) > 1L) {
    pb <- utils::txtProgressBar(min = 0L, max = nrow(sample_tbl), initial = 0L, style = 3)
    on.exit({
      if (!is.null(pb)) close(pb)
    }, add = TRUE)
  }

  out <- lapply(seq_len(nrow(sample_tbl)), function(i) {
    if (verbose) {
      message(sprintf("[GeneTrackR] Querying %s/%s: %s", i, nrow(sample_tbl), sample_tbl$sample_id[i]))
    }
    res <- query_signal_file_region(
      sample_row = sample_tbl[i],
      chrom = chrom_value,
      start = query_start,
      end = query_end
    )
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    res
  })
  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }
  out <- data.table::rbindlist(out, fill = TRUE)
  if (nrow(out) > 0L) {
    scale_tbl <- sample_tbl[, .(sample_id, scale_factor)]
    out <- merge(out, scale_tbl, by = "sample_id", all.x = TRUE)
    out[is.na(scale_factor), scale_factor := 1]
    out[, value := value * scale_factor]
    out[, scale_factor := NULL]
  }
  out <- filter_signal_strand(out, strand, sample_tbl, strand_policy = strand_policy)
  if (keep_empty_samples) {
    out <- complete_empty_signal_tracks_for_query(
      out,
      sample_ids = expected_sample_ids,
      chrom = chrom_value,
      start = query_start,
      end = query_end,
      strand = strand
    )
  }
  out
}


ensure_signal_strand_column <- function(dt) {
  if (!"strand" %in% names(dt)) {
    dt[, strand := "*"]
  }
  dt[is.na(strand) | strand == "", strand := "*"]
  dt[]
}

filter_sample_table_by_strand <- function(sample_tbl, strand, strand_policy = c("ignore_unstranded", "strict")) {
  strand_policy <- match.arg(strand_policy)
  sample_tbl <- data.table::copy(data.table::as.data.table(sample_tbl))

  if (strand %in% c("ignore", "both", "auto") || nrow(sample_tbl) == 0L) {
    return(sample_tbl)
  }

  if (!"has_strand" %in% names(sample_tbl)) {
    sample_tbl[, "has_strand" := FALSE]
  }
  if (!"strand" %in% names(sample_tbl)) {
    sample_tbl[, "strand" := "*"]
  }

  target_strand <- as.character(strand)[1L]
  sample_strand <- as.character(sample_tbl[["strand"]])
  sample_has_strand <- isTRUE_VECTOR(sample_tbl[["has_strand"]])

  if (strand_policy == "strict") {
    return(sample_tbl[sample_has_strand & sample_strand == target_strand])
  }

  sample_tbl[(!sample_has_strand) | sample_strand == target_strand]
}

isTRUE_VECTOR <- function(x) {
  out <- rep(FALSE, length(x))
  out[!is.na(x)] <- as.logical(x[!is.na(x)])
  out
}

filter_signal_strand <- function(dt, strand, sample_tbl = NULL, strand_policy = c("ignore_unstranded", "strict")) {
  strand_policy <- match.arg(strand_policy)
  dt <- ensure_signal_strand_column(dt)

  if (strand %in% c("ignore", "both", "auto") || nrow(dt) == 0L) {
    return(dt)
  }

  target_strand <- as.character(strand)[1L]

  # bigWig/wig files are unstranded by definition in this package. If the
  # selected samples are all unstranded, a requested '+' or '-' is interpreted
  # as a display/query context rather than a true strand filter. This avoids
  # dropping all signal from bigWig/wig tracks unless strict filtering is used.
  if (!is.null(sample_tbl) && "has_strand" %in% names(sample_tbl)) {
    selected <- isTRUE_VECTOR(sample_tbl[["has_strand"]])
    if (length(selected) > 0L && all(!selected, na.rm = TRUE)) {
      if (strand_policy == "strict") {
        return(dt[0])
      }
      return(dt)
    }
  }

  dt_strand <- as.character(dt[["strand"]])
  if (strand_policy == "strict") {
    return(dt[dt_strand == target_strand])
  }

  dt[dt_strand == target_strand | dt_strand == "*"]
}



complete_empty_signal_tracks_for_query <- function(dt, sample_ids, chrom, start, end, strand = "*") {
  sample_ids <- unique(as.character(sample_ids))
  dt <- data.table::as.data.table(dt)
  if (length(sample_ids) == 0L) {
    return(dt)
  }

  chrom_value <- as.character(chrom)[1L]
  start_value <- as.integer(start)[1L]
  end_value <- as.integer(end)[1L]
  strand_value <- as.character(strand)[1L]
  if (is.na(strand_value) || strand_value %in% c("ignore", "both", "auto")) {
    strand_value <- "*"
  }

  if (nrow(dt) == 0L) {
    present <- character()
    dt <- data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    )
  } else {
    dt <- ensure_signal_strand_column(dt)
    present <- unique(as.character(dt[["sample_id"]]))
  }

  missing_samples <- setdiff(sample_ids, present)
  if (length(missing_samples) == 0L) {
    return(dt[])
  }

  empty_dt <- data.table::data.table(
    sample_id = missing_samples,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
    value = 0,
    strand = strand_value
  )
  data.table::rbindlist(list(dt, empty_dt), fill = TRUE)
}

tabix_command_available <- function() {
  nzchar(Sys.which("tabix"))
}

rsamtools_tabix_available <- function() {
  requireNamespace("Rsamtools", quietly = TRUE)
}

get_tabix_backend <- function() {
  if (tabix_command_available()) {
    return("system")
  }
  if (rsamtools_tabix_available()) {
    return("Rsamtools")
  }
  NA_character_
}

has_tabix_index <- function(file) {
  file.exists(paste0(file, ".tbi"))
}

empty_signal_dt <- function() {
  data.table::data.table(
    sample_id = character(), chrom = character(), start = integer(),
    end = integer(), value = numeric(), strand = character()
  )
}

read_tabix_lines <- function(file, region, backend = NA_character_) {
  backend <- as.character(backend)[1L]
  if (is.na(backend) || !nzchar(backend)) {
    backend <- get_tabix_backend()
  }

  if (identical(backend, "system")) {
    lines <- tryCatch(
      system2("tabix", args = c(file, region), stdout = TRUE, stderr = TRUE),
      error = function(e) structure(character(), status = 1L)
    )
    status <- attr(lines, "status")
    ok <- is.null(status) || identical(as.integer(status), 0L)
    return(list(lines = as.character(lines), ok = ok, backend = "system", region = region))
  }

  if (identical(backend, "Rsamtools")) {
    ans <- tryCatch(
      Rsamtools::scanTabix(file = file, param = region),
      error = function(e) e
    )
    if (inherits(ans, "error")) {
      return(list(lines = character(), ok = FALSE, backend = "Rsamtools", region = region))
    }
    return(list(lines = unlist(ans, use.names = FALSE), ok = TRUE, backend = "Rsamtools", region = region))
  }

  list(lines = character(), ok = FALSE, backend = NA_character_, region = region)
}

query_bedgraph_tabix <- function(file, sample_id, chrom, start, end, strand = "*", backend = NA_character_, empty_fallback = FALSE) {
  chrom_value <- as.character(chrom)
  query_start <- as.integer(start)
  query_end <- as.integer(end)
  region <- sprintf("%s:%s-%s", chrom_value, query_start, query_end)
  empty <- empty_signal_dt()

  parse_bedgraph_lines <- function(lines) {
    if (length(lines) == 0L) {
      return(empty)
    }
    out <- tryCatch(
      attach_signal_sample(
        parse_bedgraph_lines_native(lines),
        sample_id = sample_id,
        strand = strand
      ),
      error = function(e) empty
    )
    subset_signal_region(
      out,
      chrom = chrom_value,
      start = query_start,
      end = query_end
    )
  }

  # Tabix returns no lines when a region has no coverage. This is a valid and
  # very common result, especially for sparse bedGraph tracks. Therefore an
  # empty indexed result must NOT trigger a full-file fread by default. The
  # expensive full-file check is only used when the backend fails or when
  # `empty_fallback = TRUE` is explicitly requested for debugging.
  res <- read_tabix_lines(file = file, region = region, backend = backend)
  if (isTRUE(res$ok)) {
    out <- parse_bedgraph_lines(res$lines)
    if (nrow(out) > 0L || !isTRUE(empty_fallback)) {
      return(out)
    }
  }

  if (query_start > 1L && (!isTRUE(res$ok) || isTRUE(empty_fallback))) {
    region0 <- sprintf("%s:%s-%s", chrom_value, query_start - 1L, query_end)
    res0 <- read_tabix_lines(file = file, region = region0, backend = backend)
    if (isTRUE(res0$ok)) {
      out0 <- parse_bedgraph_lines(res0$lines)
      if (nrow(out0) > 0L || !isTRUE(empty_fallback)) {
        return(out0)
      }
    }
  }

  query_bedgraph_full_file(
    file = file,
    sample_id = sample_id,
    chrom = chrom_value,
    start = query_start,
    end = query_end,
    strand = strand
  )
}

query_bedgraph_full_file <- function(file, sample_id, chrom, start, end, strand = "*") {
  chrom_value <- as.character(chrom)
  query_start <- as.integer(start)
  query_end <- as.integer(end)
  empty <- empty_signal_dt()
  dt <- tryCatch(
    read_signal_file_memory(
      file = file,
      format = "bedgraph",
      sample_id = sample_id,
      strand = strand
    ),
    error = function(e) empty
  )
  subset_signal_region(
    dt,
    chrom = chrom_value,
    start = query_start,
    end = query_end
  )
}


#' Diagnose tabix availability for bedGraph tracks
#'
#' @param object A BwgTrack object created by `read_bwg()`.
#' @return A data.table reporting index presence, backend availability, and whether indexed querying is enabled.
#' @examples
#' rnaseq <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph",
#'   sample_names = "RNA_seq_plus",
#'   strand = "+",
#'   mode = "lazy",
#'   verbose = FALSE
#' )
#' diagnose_tabix(rnaseq)
#' @export
diagnose_tabix <- function(object) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  sample_tbl <- data.table::copy(object$samples)
  if (!"has_tabix" %in% names(sample_tbl)) {
    sample_tbl[, "has_tabix" := vapply(file, has_tabix_index, logical(1))]
  }
  if (!"use_tabix" %in% names(sample_tbl)) {
    sample_tbl[, "use_tabix" := FALSE]
  }
  if (!"tabix_backend" %in% names(sample_tbl)) {
    sample_tbl[, "tabix_backend" := NA_character_]
  }
  sample_tbl[, system_tabix := tabix_command_available()]
  sample_tbl[, rsamtools_tabix := rsamtools_tabix_available()]
  sample_tbl[, available_backend := get_tabix_backend()]
  sample_tbl[, .(
    sample_id,
    file,
    format,
    has_tabix,
    use_tabix,
    tabix_backend,
    tabix_empty_fallback = if ("tabix_empty_fallback" %in% names(sample_tbl)) tabix_empty_fallback else FALSE,
    system_tabix,
    rsamtools_tabix,
    available_backend
  )]
}

#' Test whether a bedGraph tabix query returns records without full-file fallback
#'
#' @param object A BwgTrack object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param samples Optional sample IDs.
#' @return A data.table with per-sample query timing and number of records.
#' @examples
#' signal_file <- system.file(
#'   "extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"
#' )
#' signal <- read_bwg(
#'   signal_file, format = "bedgraph", sample_names = "RNA_plus",
#'   strand = "+", mode = "lazy", verbose = FALSE
#' )
#' benchmark_tabix_query(signal, "chr1", 12339001L, 12340000L)
#' @export
benchmark_tabix_query <- function(object, chrom, start, end, samples = NULL) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  sample_tbl <- data.table::copy(object$samples)
  if (!is.null(samples)) {
    sample_tbl <- sample_tbl[sample_id %in% as.character(samples)]
  }
  sample_tbl <- sample_tbl[format == "bedgraph"]
  if (nrow(sample_tbl) == 0L) {
    return(data.table::data.table())
  }
  out <- lapply(seq_len(nrow(sample_tbl)), function(i) {
    t <- system.time({
      dt <- if (isTRUE(sample_tbl$use_tabix[i])) {
        query_bedgraph_tabix(
          file = sample_tbl$file[i],
          sample_id = sample_tbl$sample_id[i],
          chrom = chrom,
          start = start,
          end = end,
          strand = sample_tbl$strand[i],
          backend = sample_tbl$tabix_backend[i],
          empty_fallback = FALSE
        )
      } else {
        query_bedgraph_full_file(
          file = sample_tbl$file[i],
          sample_id = sample_tbl$sample_id[i],
          chrom = chrom,
          start = start,
          end = end,
          strand = sample_tbl$strand[i]
        )
      }
    })
    data.table::data.table(
      sample_id = sample_tbl$sample_id[i],
      use_tabix = isTRUE(sample_tbl$use_tabix[i]),
      tabix_backend = as.character(sample_tbl$tabix_backend[i]),
      n_records = nrow(dt),
      elapsed_sec = unname(t[["elapsed"]])
    )
  })
  data.table::rbindlist(out, fill = TRUE)
}

#' Normalize signal tracks
#'
#' @param object A BwgTrack object.
#' @param method Normalization method. Use `RPM`/`CPM` with `library_size`, `scale` with `scale_factor`, or `custom` with `custom_factor`.
#' @param library_size Optional named library sizes.
#' @param scale_factor Optional named scale factors.
#' @param custom_factor Optional named custom factors.
#' @return A normalized BwgTrack object.
#' @examples
#' rnaseq <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph",
#'   sample_names = "RNA_seq_plus",
#'   strand = "+",
#'   mode = "memory",
#'   verbose = FALSE
#' )
#' norm_bwg(rnaseq, method = "scale", scale_factor = c(RNA_seq_plus = 0.5))
#' @export
norm_bwg <- function(object, method = c("none", "RPM", "CPM", "scale", "custom"), library_size = NULL, scale_factor = NULL, custom_factor = NULL) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")
  method <- match.arg(method)
  obj <- object

  if (method == "none") {
    return(obj)
  }

  samples <- data.table::copy(obj$samples)
  factor <- rep(1, nrow(samples))
  names(factor) <- samples$sample_id

  if (method %in% c("RPM", "CPM")) {
    if (is.null(library_size)) {
      if (is.null(obj$data)) {
        stop("`library_size` is required for lazy signal normalization.", call. = FALSE)
      }
      library_size <- obj$data[, .(library_size = sum(value, na.rm = TRUE)), by = sample_id]
      library_size <- stats::setNames(library_size$library_size, library_size$sample_id)
    }
    factor[names(library_size)] <- 1e6 / as.numeric(library_size)
  }

  if (method == "scale") {
    stop_if_not(!is.null(scale_factor), "`scale_factor` is required when `method = 'scale'`.")
    factor[names(scale_factor)] <- as.numeric(scale_factor)
  }

  if (method == "custom") {
    stop_if_not(!is.null(custom_factor), "`custom_factor` is required when `method = 'custom'`.")
    factor[names(custom_factor)] <- as.numeric(custom_factor)
  }

  samples[, scale_factor := factor[sample_id]]
  samples[, norm_method := method]
  obj$samples <- samples

  if (!is.null(obj$data)) {
    scale_dt <- samples[, .(sample_id, scale_factor)]
    obj$data <- merge(obj$data, scale_dt, by = "sample_id", all.x = TRUE)
    obj$data[, value := value * scale_factor]
    obj$data[, scale_factor := NULL]
  }
  obj
}

#' Summarize signal tracks within a region
#'
#' @param object A BwgTrack object.
#' @param chrom Optional chromosome name.
#' @param start Optional region start.
#' @param end Optional region end.
#' @param samples Optional sample IDs.
#' @return A data.table summary.
#' @examples
#' signal_file <- system.file(
#'   "extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"
#' )
#' signal <- read_bwg(
#'   signal_file, format = "bedgraph", sample_names = "RNA_plus",
#'   strand = "+", mode = "memory", verbose = FALSE
#' )
#' summary_bwg(signal)
#' @export
summary_bwg <- function(object, chrom = NULL, start = NULL, end = NULL, samples = NULL) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")

  if (is.null(chrom)) {
    if (is.null(object$data)) {
      return(data.table::copy(object$samples))
    }
    dt <- data.table::copy(object$data)
  } else {
    dt <- .query_bwg_internal(object, chrom, start, end, samples = samples)
  }

  if (!is.null(samples)) {
    dt <- dt[sample_id %in% samples]
  }
  if (nrow(dt) == 0L) {
    return(data.table::data.table(sample_id = character(), sum_signal = numeric(), mean_signal = numeric(), max_signal = numeric(), covered_bases = integer()))
  }
  dt[, width := end - start + 1L]
  dt[, .(
    sum_signal = sum(value * width, na.rm = TRUE),
    mean_signal = stats::weighted.mean(value, width, na.rm = TRUE),
    max_signal = max(value, na.rm = TRUE),
    covered_bases = sum(width[value > 0], na.rm = TRUE)
  ), by = sample_id]
}

#' Slice signal tracks by genomic region
#'
#' @param object A BwgTrack object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param samples Optional sample IDs.
#' @param strand Strand selector. For unstranded bigWig/wig tracks, '+' and '-' do not filter records.
#' @param as Return type. Use `BwgTrack` to keep the object structure, `data.frame` for a plain table, or `GRanges` for genomic interval operations.
#' @return A BwgTrack, data.table, or GRanges object.
#' @examples
#' rnaseq <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph", sample_names = "RNA_seq_plus", strand = "+",
#'   mode = "memory", verbose = FALSE
#' )
#' slice_bwg(
#'   rnaseq, chrom = "chr1", start = 12339001, end = 12352000,
#'   strand = "+", as = "BwgTrack"
#' )
#' @export
slice_bwg <- function(object, chrom, start, end, samples = NULL, strand = "ignore", as = c("BwgTrack", "data.frame", "GRanges")) {
  as <- match.arg(as)
  dt <- .query_bwg_internal(object, chrom, start, end, samples = samples, strand = strand)

  if (as == "data.frame") {
    return(dt)
  }
  if (as == "GRanges") {
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(dt$start, dt$end),
      sample_id = dt$sample_id,
      score = dt$value,
      strand = dt$strand
    ))
  }

  samples_tbl <- object$samples
  if (!is.null(samples)) {
    samples_tbl <- samples_tbl[sample_id %in% samples]
  }
  seqinfo <- subset_bwg_seqinfo(
    object,
    sample_ids = samples_tbl[["sample_id"]],
    chrom = chrom
  )
  BwgTrack(
    samples = samples_tbl,
    data = dt,
    meta = modifyList(object$meta, list(mode = "memory")),
    validation = make_empty_validation(),
    seqinfo = seqinfo
  )
}

#' Merge signal track objects
#'
#' @param ... BwgTrack objects.
#' @param sample_conflict Conflict strategy for duplicated sample IDs. `error` stops, `rename` makes sample IDs unique, `sum`/`mean` combines records with identical genomic intervals, and `keep_first` keeps the first occurrence.
#' @param require_same_norm Require identical normalization labels before merging.
#' @return A merged BwgTrack object.
#' @examples
#' plus <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph", sample_names = "RNA_seq_plus", strand = "+",
#'   mode = "memory", verbose = FALSE
#' )
#' minus <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_minus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph", sample_names = "RNA_seq_minus", strand = "-",
#'   mode = "memory", verbose = FALSE
#' )
#' signal_all <- merge_bwg(plus, minus)
#' signal_all
#' @export
merge_bwg <- function(..., sample_conflict = c("error", "rename", "sum", "mean", "keep_first"), require_same_norm = TRUE) {
  objs <- list(...)
  stop_if_not(length(objs) >= 1L, "At least one BwgTrack object is required.")
  stop_if_not(all(vapply(objs, inherits, logical(1), "BwgTrack")), "All inputs must be BwgTrack objects.")
  sample_conflict <- match.arg(sample_conflict)

  samples <- data.table::rbindlist(lapply(seq_along(objs), function(i) {
    x <- data.table::copy(objs[[i]]$samples)
    x[, source_index := i]
    x
  }), fill = TRUE)

  if (require_same_norm && data.table::uniqueN(samples$norm_method) > 1L) {
    stop("Signal tracks have different normalization methods.", call. = FALSE)
  }

  dup <- samples[duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE), unique(sample_id)]
  if (length(dup) > 0L && sample_conflict == "error") {
    stop("Duplicated sample IDs found. Use `sample_conflict` to resolve them.", call. = FALSE)
  }
  if (length(dup) > 0L && sample_conflict == "rename") {
    samples[sample_id %in% dup, sample_id := paste0("set", source_index, "_", sample_id)]
  }

  data_list <- lapply(seq_along(objs), function(i) {
    dt <- data.table::copy(objs[[i]]$data)
    if (is.null(dt)) return(NULL)
    dt[, source_index := i]
    dt
  })
  data <- data.table::rbindlist(data_list, fill = TRUE)

  seqinfo_list <- lapply(seq_along(objs), function(i) {
    si <- objs[[i]]$seqinfo
    if (is.null(si)) return(NULL)
    si <- data.table::copy(data.table::as.data.table(si))
    si[, source_index := i]
    si
  })
  seqinfo <- data.table::rbindlist(seqinfo_list, fill = TRUE)
  if (is.null(seqinfo) || nrow(seqinfo) == 0L) {
    seqinfo <- NULL
  }

  if (length(dup) > 0L && sample_conflict == "keep_first") {
    keep_sources <- samples[, .(source_index = source_index[1L]), by = sample_id]
    keep_key <- paste(keep_sources$sample_id, keep_sources$source_index, sep = "\r")
    sample_key <- paste(samples$sample_id, samples$source_index, sep = "\r")
    samples <- samples[sample_key %in% keep_key]
    if (!is.null(data)) {
      data_key <- paste(data$sample_id, data$source_index, sep = "\r")
      data <- data[data_key %in% keep_key]
    }
    if (!is.null(seqinfo)) {
      seqinfo_key <- paste(seqinfo$sample_id, seqinfo$source_index, sep = "\r")
      seqinfo <- seqinfo[seqinfo_key %in% keep_key]
    }
  }

  if (!is.null(data) && length(dup) > 0L) {
    if (sample_conflict == "rename") {
      map <- samples[, .(source_index, old_sample_id = sub("^set[0-9]+_", "", sample_id), sample_id)]
      data <- merge(data, map, by.x = c("source_index", "sample_id"), by.y = c("source_index", "old_sample_id"), all.x = TRUE)
      data[!is.na(i.sample_id), sample_id := i.sample_id]
      data[, i.sample_id := NULL]
    }
    if (sample_conflict %in% c("sum", "mean")) {
      fun <- if (sample_conflict == "sum") sum else mean
      data <- data[, .(value = fun(value, na.rm = TRUE)), by = .(sample_id, chrom, start, end, strand)]
      samples <- samples[, .SD[1L], by = sample_id]
    }
  }

  if (!is.null(seqinfo) && length(dup) > 0L) {
    if (sample_conflict == "rename") {
      map <- samples[, .(source_index, old_sample_id = sub("^set[0-9]+_", "", sample_id), sample_id)]
      seqinfo <- merge(
        seqinfo,
        map,
        by.x = c("source_index", "sample_id"),
        by.y = c("source_index", "old_sample_id"),
        all.x = TRUE
      )
      seqinfo[!is.na(i.sample_id), sample_id := i.sample_id]
      seqinfo[, i.sample_id := NULL]
    }
    if (sample_conflict %in% c("sum", "mean", "keep_first")) {
      seqinfo <- unique(seqinfo[, .(sample_id, chrom, length)])
    }
  }

  if (!is.null(data) && "source_index" %in% names(data)) {
    data[, source_index := NULL]
  }
  if (!is.null(seqinfo) && "source_index" %in% names(seqinfo)) {
    seqinfo[, source_index := NULL]
  }
  if ("source_index" %in% names(samples)) {
    samples[, source_index := NULL]
  }

  BwgTrack(
    samples = samples,
    data = data,
    meta = list(
      mode = if (is.null(data)) "lazy" else "memory",
      coordinate = "1-based closed"
    ),
    validation = make_empty_validation(),
    seqinfo = seqinfo
  )
}

#' Bin signal tracks into fixed-width windows
#'
#' @param data Signal table returned by retrieve_bwg. The table must contain `sample_id`, `chrom`, `start`, `end`, and `value` columns.
#' @param bin_size Bin size in bases.
#' @return A binned signal table.
#' @examples
#' rnaseq <- read_bwg(
#'   system.file("extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"),
#'   format = "bedgraph", sample_names = "RNA_seq_plus", strand = "+",
#'   mode = "memory", verbose = FALSE
#' )
#' dt <- retrieve_bwg(
#'   rnaseq, chrom = "chr1", start = 12339001, end = 12352000, strand = "+"
#' )
#' bin_bwg(dt, bin_size = 50)
#' @export
bin_bwg <- function(data, bin_size = 50L) {
  dt <- data.table::copy(data.table::as.data.table(data))
  stop_if_not(all(c("sample_id", "chrom", "start", "end", "value") %in% names(dt)), "Input must be a signal table.")
  if (nrow(dt) == 0L) {
    out <- data.table::data.table(
      sample_id = character(), chrom = character(), start = integer(),
      end = integer(), value = numeric(), bin = integer()
    )
    return(out)
  }
  dt <- dt[is.finite(as.numeric(start)) & is.finite(as.numeric(end))]
  if (nrow(dt) == 0L) {
    out <- data.table::data.table(
      sample_id = character(), chrom = character(), start = integer(),
      end = integer(), value = numeric(), bin = integer()
    )
    return(out)
  }
  dt[, bin := floor((start - min(start)) / bin_size) + 1L, by = .(sample_id, chrom)]
  dt[, .(
    start = min(start),
    end = max(end),
    value = mean(value, na.rm = TRUE)
  ), by = .(sample_id, chrom, bin)]
}


#' Get chromosome information from BwgTrack objects
#'
#' @param object A BwgTrack object.
#' @param samples Optional sample IDs.
#' @return A data.table containing chromosome names and lengths for each sample.
#' @details
#' Schema-v2 objects return stored sequence metadata directly. Legacy objects
#' without a `seqinfo` slot query chromosome metadata from bigWig source files
#' through the built-in pure-R backend.
#' @examples
#' signal_file <- system.file(
#'   "extdata", "gtr_demo_rnaseq_plus.bedgraph", package = "GeneTrackR"
#' )
#' signal <- read_bwg(
#'   signal_file, format = "bedgraph", sample_names = "RNA_plus",
#'   strand = "+", mode = "memory", verbose = FALSE
#' )
#' seqinfo_bwg(signal)
#' @export
seqinfo_bwg <- function(object, samples = NULL) {
  stop_if_not(inherits(object, "BwgTrack"), "`object` must be a BwgTrack object.")

  if (!is.null(object$seqinfo)) {
    return(subset_bwg_seqinfo(object, sample_ids = samples))
  }

  sample_tbl <- data.table::copy(object$samples)
  if (!"use_tabix" %in% names(sample_tbl)) {
    sample_tbl[, "use_tabix" := FALSE]
  }
  if (!"has_strand" %in% names(sample_tbl)) {
    sample_tbl[, "has_strand" := FALSE]
  }
  if (!is.null(samples)) {
    sample_tbl <- sample_tbl[sample_id %in% samples]
  }
  if (nrow(sample_tbl) == 0L) {
    return(data.table::data.table(
      sample_id = character(),
      chrom = character(),
      length = integer()
    ))
  }

  out <- lapply(seq_len(nrow(sample_tbl)), function(i) {
    if (!identical(sample_tbl$format[i], "bigwig")) {
      return(NULL)
    }

    si <- bigwig_seqinfo_native(sample_tbl$file[i])
    if (nrow(si) == 0L) {
      return(NULL)
    }
    si[, sample_id := sample_tbl$sample_id[i]]
    data.table::setcolorder(si, c("sample_id", "chrom", "length"))
    si[]
  })

  data.table::rbindlist(out, fill = TRUE)
}
