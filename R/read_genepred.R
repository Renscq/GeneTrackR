# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.29
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
#' @param progress Logical. Whether to show a progress bar for exon parsing. By default, a progress bar is shown only in interactive sessions when `verbose = TRUE`.
#' @details
#' Standard GenePred uses `name` as transcript ID. GenePredExt commonly uses
#' `name` as transcript ID and `name2` as gene ID, so the default `gene_col`
#' tries `name2` first and falls back to `name` when `name2` is unavailable.
#' `coordinate = "ucsc"` converts 0-based half-open GenePred coordinates to
#' the package's internal 1-based closed coordinate system.
#'
#' For large files, `verbose = TRUE` prints major processing stages, while
#' `progress = TRUE` additionally shows a text progress bar for exon parsing.
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

  stop_if_not(file.exists(file), paste0("File does not exist: ", file))

  progress_msg <- function(..., appendLF = TRUE) {
    if (verbose) {
      message(sprintf(...), appendLF = appendLF)
    }
  }

  input_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  file_size_mb <- suppressWarnings(file.info(file)$size / 1024^2)
  if (is.finite(file_size_mb)) {
    progress_msg("[GeneTrackR] Reading GenePred file: %s (%.2f MB)", input_file, file_size_mb)
  } else {
    progress_msg("[GeneTrackR] Reading GenePred file: %s", input_file)
  }

  dt <- data.table::fread(
    file,
    header = FALSE,
    sep = "\t",
    data.table = TRUE,
    showProgress = verbose
  )
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

  progress_msg("[GeneTrackR] Parsing exon structures for %s transcripts.", format(nrow(dt), big.mark = ","))

  n_record <- nrow(dt)
  exon_list <- vector("list", n_record)
  pb <- NULL
  if (progress && n_record > 0L) {
    pb <- utils::txtProgressBar(min = 0L, max = n_record, initial = 0L, style = 3)
    on.exit({
      if (!is.null(pb)) {
        close(pb)
      }
    }, add = TRUE)
  }
  progress_step <- max(1L, floor(n_record / 100L))

  for (i in seq_len(n_record)) {
    starts <- parse_comma_integer(dt$exonStarts[i])
    ends <- parse_comma_integer(dt$exonEnds[i])
    if (coordinate == "ucsc") {
      starts <- starts + 1L
    }
    frames <- if ("exonFrames" %in% names(dt)) parse_comma_integer(dt$exonFrames[i]) else rep(NA_integer_, length(starts))
    n <- max(length(starts), length(ends))
    exon_list[[i]] <- data.table::data.table(
      row_id = dt$row_id[i],
      transcript_id = dt$transcript_id[i],
      gene_id = dt$gene_id[i],
      chrom = as.character(dt$chrom[i]),
      strand = as.character(dt$strand[i]),
      exon_number = seq_len(n),
      exon_start = starts,
      exon_end = ends,
      exon_frame = if (length(frames) == n) frames else rep(NA_integer_, n)
    )
    if (!is.null(pb) && (i %% progress_step == 0L || i == n_record)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }

  progress_msg("[GeneTrackR] Combining exon records.")
  exons <- data.table::rbindlist(exon_list, fill = TRUE)
  progress_msg("[GeneTrackR] Parsed %s exons.", format(nrow(exons), big.mark = ","))

  progress_msg("[GeneTrackR] Validating transcript and exon tables.")
  validation <- validate_genepred_tables(tx, exons)

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
