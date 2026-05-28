# Author: Rensc
# Date: 2026-05-26
# Version: 0.2.5
# Function: Read GenePred or GenePredExt annotation files
# Input: GenePred-format text file
# Output: GenePred object with transcript, exon, and gene tables

#' Read a GenePred annotation file
#'
#' @param file GenePred or GenePredExt file path.
#' @param format Input format. Use auto, genePred, or genePredExt.
#' @param coordinate Input coordinate system. ucsc means 0-based half-open.
#' @param remove_invalid Remove invalid records after validation.
#' @param report_invalid Store invalid records in the object validation slot.
#' @param gene_col Gene identifier source. name2 is recommended for GenePredExt.
#' @param transcript_col Transcript identifier source.
#' @param verbose Logical. Whether to print step-level progress messages. Default TRUE.
#' @param progress Logical. Whether to show a stage-level progress bar. By default, a progress bar is shown only in interactive sessions when `verbose = TRUE`.
#' @details
#' Standard GenePred uses `name` as transcript ID. GenePredExt commonly uses
#' `name` as transcript ID and `name2` as gene ID, so the default `gene_col`
#' tries `name2` first and falls back to `name` when `name2` is unavailable.
#' `coordinate = "ucsc"` converts 0-based half-open GenePred coordinates to
#' the package's internal 1-based closed coordinate system.
#'
#' For large files, `verbose = TRUE` prints major processing stages, while
#' `progress = TRUE` additionally shows a stage-level text progress bar.
#' @return A GenePred object.
#' @examples
#' \dontrun{
#' gp <- read_genepred("annotation.genePredExt", format = "genePredExt")
#' gp_quiet <- read_genepred(
#'   "annotation.genePredExt",
#'   format = "genePredExt",
#'   verbose = FALSE
#' )
#' gp_custom <- read_genepred(
#'   "annotation.genePredExt",
#'   format = "genePredExt",
#'   gene_col = "name2",
#'   transcript_col = "name",
#'   remove_invalid = TRUE
#' )
#' }
#' @export
read_genepred <- function(file,
                          format = c("auto", "genePred", "genePredExt"),
                          coordinate = c("ucsc", "granges"),
                          remove_invalid = TRUE,
                          report_invalid = TRUE,
                          gene_col = c("name2", "name"),
                          transcript_col = "name",
                          verbose = TRUE,
                          progress = interactive() && isTRUE(verbose)) {
  format <- match.arg(format)
  coordinate <- match.arg(coordinate)
  gene_col <- match.arg(gene_col)
  verbose <- isTRUE(verbose)
  progress <- isTRUE(progress)
  old_dt_verbose <- getOption("datatable.verbose")
  options(datatable.verbose = FALSE)
  on.exit(options(datatable.verbose = old_dt_verbose), add = TRUE)

  stop_if_not(file.exists(file), paste0("File does not exist: ", file))

  progress_msg <- make_progress_message(verbose)
  stage <- make_stage_progress(total = 6L, progress = progress)
  on.exit(stage$close(), add = TRUE)

  input_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  progress_msg("%s", format_file_size_message("GenePred", file))

  dt <- data.table::fread(
    file,
    header = FALSE,
    sep = "\t",
    data.table = TRUE,
    showProgress = FALSE
  )
  stage$tick()
  n_col <- ncol(dt)
  progress_msg("[GeneTrackR] Loaded %s records with %s columns.", format(nrow(dt), big.mark = ","), n_col)

  if (format == "auto") {
    format <- if (n_col >= 15L) "genePredExt" else "genePred"
  }
  progress_msg("[GeneTrackR] Detected input format: %s", format)

  standard_cols <- c(
    "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
    "exonCount", "exonStarts", "exonEnds"
  )
  ext_cols <- c("score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

  if (format == "genePred") {
    stop_if_not(n_col >= 10L, "A genePred file must contain at least 10 columns.")
    data.table::setnames(dt, seq_len(10L), standard_cols)
  } else {
    stop_if_not(n_col >= 15L, "A genePredExt file must contain at least 15 columns.")
    data.table::setnames(dt, seq_len(15L), c(standard_cols, ext_cols))
  }

  stage$tick()

  if (!transcript_col %in% names(dt)) {
    stop("`transcript_col` is not available in the input file.", call. = FALSE)
  }

  if (!gene_col %in% names(dt)) {
    gene_col <- "name"
  }

  progress_msg("[GeneTrackR] Building transcript-level table.")
  dt[, transcript_id := as.character(get(transcript_col))]
  dt[, gene_id := as.character(get(gene_col))]
  dt[is.na(gene_id) | gene_id == "", gene_id := transcript_id]

  dt[, row_id := .I]
  dt[, gene_type := data.table::fifelse(as.integer(cdsStart) < as.integer(cdsEnd), "coding", "non-coding")]

  if (coordinate == "ucsc") {
    dt[, `:=`(
      tx_start = as.integer(txStart) + 1L,
      tx_end = as.integer(txEnd),
      cds_start = as.integer(cdsStart) + 1L,
      cds_end = as.integer(cdsEnd)
    )]
  } else {
    dt[, `:=`(
      tx_start = as.integer(txStart),
      tx_end = as.integer(txEnd),
      cds_start = as.integer(cdsStart),
      cds_end = as.integer(cdsEnd)
    )]
  }

  tx <- dt[, .(
    row_id,
    transcript_id,
    gene_id,
    chrom = as.character(chrom),
    strand = as.character(strand),
    tx_start,
    tx_end,
    cds_start,
    cds_end,
    exon_count = as.integer(exonCount),
    gene_type,
    score = if ("score" %in% names(dt)) suppressWarnings(as.numeric(score)) else NA_real_,
    cds_start_stat = if ("cdsStartStat" %in% names(dt)) as.character(cdsStartStat) else NA_character_,
    cds_end_stat = if ("cdsEndStat" %in% names(dt)) as.character(cdsEndStat) else NA_character_,
    exon_frames = if ("exonFrames" %in% names(dt)) as.character(exonFrames) else NA_character_
  )]
  stage$tick()

  progress_msg("[GeneTrackR] Parsing exon structures for %s transcripts.", format(nrow(dt), big.mark = ","))
  exons <- parse_genepred_exons_fast(dt = dt, coordinate = coordinate, has_frames = "exonFrames" %in% names(dt))
  stage$tick()
  progress_msg("[GeneTrackR] Parsed %s exons.", format(nrow(exons), big.mark = ","))

  progress_msg("[GeneTrackR] Validating transcript and exon tables.")
  validation <- validate_genepred_tables(tx, exons)
  stage$tick()

  n_invalid <- nrow(validation$invalid_records)
  if (n_invalid > 0L) {
    progress_msg("[GeneTrackR] Found %s invalid validation records.", format(n_invalid, big.mark = ","))
  } else {
    progress_msg("[GeneTrackR] No invalid records detected.")
  }

  if (remove_invalid && n_invalid > 0L) {
    bad_rows <- unique(validation$invalid_records$row_id)
    progress_msg("[GeneTrackR] Removing %s invalid transcripts.", format(length(bad_rows), big.mark = ","))
    tx <- tx[!row_id %in% bad_rows]
    exons <- exons[!row_id %in% bad_rows]
  }

  tx[, row_id := NULL]
  exons[, row_id := NULL]

  progress_msg("[GeneTrackR] Building gene-level table.")
  genes <- build_gene_table(tx)
  stage$tick()
  meta <- list(
    source_file = input_file,
    format = format,
    coordinate_input = coordinate,
    coordinate_internal = "1-based closed",
    gene_col = gene_col,
    transcript_col = transcript_col
  )

  if (!report_invalid) {
    validation <- make_empty_validation()
  }

  progress_msg(
    "[GeneTrackR] Finished. genes: %s; transcripts: %s; exons: %s.",
    format(nrow(genes), big.mark = ","),
    format(nrow(tx), big.mark = ","),
    format(nrow(exons), big.mark = ",")
  )

  GenePred(transcripts = tx, exons = exons, genes = genes, meta = meta, validation = validation)
}



parse_genepred_exons_fast <- function(dt, coordinate = c("ucsc", "granges"), has_frames = FALSE) {
  coordinate <- match.arg(coordinate)
  n_record <- nrow(dt)
  if (n_record == 0L) {
    return(data.table::data.table(
      row_id = integer(),
      transcript_id = character(),
      gene_id = character(),
      chrom = character(),
      strand = character(),
      exon_number = integer(),
      exon_start = integer(),
      exon_end = integer(),
      exon_frame = integer()
    ))
  }

  split_integer_fields <- function(x) {
    fields <- strsplit(as.character(x), ",", fixed = TRUE)
    lapply(fields, function(v) {
      v <- v[nzchar(v)]
      if (length(v) == 0L) integer() else as.integer(v)
    })
  }

  starts_list <- split_integer_fields(dt[["exonStarts"]])
  ends_list <- split_integer_fields(dt[["exonEnds"]])
  n_start <- lengths(starts_list)
  n_end <- lengths(ends_list)

  # Fast path for valid GenePred rows. Invalid rows with mismatched exonStarts
  # and exonEnds are rare; they are handled by the slow fallback so validation
  # can still report them cleanly.
  if (!identical(n_start, n_end)) {
    return(parse_genepred_exons_fallback(dt = dt, coordinate = coordinate, has_frames = has_frames))
  }

  total_exons <- sum(n_start)
  if (total_exons == 0L) {
    return(data.table::data.table(
      row_id = integer(),
      transcript_id = character(),
      gene_id = character(),
      chrom = character(),
      strand = character(),
      exon_number = integer(),
      exon_start = integer(),
      exon_end = integer(),
      exon_frame = integer()
    ))
  }

  exon_start <- unlist(starts_list, use.names = FALSE)
  exon_end <- unlist(ends_list, use.names = FALSE)
  if (coordinate == "ucsc") {
    exon_start <- exon_start + 1L
  }

  exon_frame <- rep.int(NA_integer_, total_exons)
  if (isTRUE(has_frames)) {
    frames_list <- split_integer_fields(dt[["exonFrames"]])
    n_frame <- lengths(frames_list)
    if (identical(n_frame, n_start)) {
      exon_frame <- unlist(frames_list, use.names = FALSE)
    } else {
      # Preserve valid frame rows and leave malformed rows as NA.
      offset <- cumsum(c(0L, n_start))
      valid <- which(n_frame == n_start & n_start > 0L)
      for (i in valid) {
        exon_frame[(offset[i] + 1L):offset[i + 1L]] <- frames_list[[i]]
      }
    }
  }

  strand_rep <- rep.int(as.character(dt[["strand"]]), n_start)
  exon_number_genomic <- sequence(n_start)
  exon_count_rep <- rep.int(n_start, n_start)
  exon_number_tx <- data.table::fifelse(
    strand_rep == "-",
    as.integer(exon_count_rep - exon_number_genomic + 1L),
    as.integer(exon_number_genomic)
  )

  data.table::data.table(
    row_id = rep.int(dt[["row_id"]], n_start),
    transcript_id = rep.int(dt[["transcript_id"]], n_start),
    gene_id = rep.int(dt[["gene_id"]], n_start),
    chrom = rep.int(as.character(dt[["chrom"]]), n_start),
    strand = strand_rep,
    exon_number = as.integer(exon_number_tx),
    exon_start = as.integer(exon_start),
    exon_end = as.integer(exon_end),
    exon_frame = as.integer(exon_frame)
  )
}

parse_genepred_exons_fallback <- function(dt, coordinate = c("ucsc", "granges"), has_frames = FALSE) {
  coordinate <- match.arg(coordinate)
  n_record <- nrow(dt)
  exon_list <- vector("list", n_record)

  for (i in seq_len(n_record)) {
    starts <- parse_comma_integer(dt$exonStarts[i])
    ends <- parse_comma_integer(dt$exonEnds[i])
    if (coordinate == "ucsc") {
      starts <- starts + 1L
    }
    frames <- if (isTRUE(has_frames)) parse_comma_integer(dt$exonFrames[i]) else rep(NA_integer_, length(starts))
    n <- max(length(starts), length(ends))
    exon_list[[i]] <- data.table::data.table(
      row_id = dt$row_id[i],
      transcript_id = dt$transcript_id[i],
      gene_id = dt$gene_id[i],
      chrom = as.character(dt$chrom[i]),
      strand = as.character(dt$strand[i]),
      exon_number = if (as.character(dt$strand[i]) == "-") rev(seq_len(n)) else seq_len(n),
      exon_start = if (length(starts) == n) starts else rep(NA_integer_, n),
      exon_end = if (length(ends) == n) ends else rep(NA_integer_, n),
      exon_frame = if (length(frames) == n) frames else rep(NA_integer_, n)
    )
  }

  data.table::rbindlist(exon_list, fill = TRUE)
}
