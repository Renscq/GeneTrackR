# Author: Rensc
# Date: 2026-08-31
# Version: dev003
# Function: Read BED track files into unified Feature objects
# Input: BED interval files
# Output: Feature objects

#' Read a BED file as a FeatureTrack
#'
#' @description
#' Reads BED3-BED12 style files and converts BED coordinates from 0-based
#' half-open to 1-based closed intervals by default.
#'
#' @param file BED file path. Gzip-compressed files are supported by `data.table::fread()`.
#' @param coordinate Input coordinate system. `bed` means 0-based half-open;
#' `granges` means 1-based closed.
#' @param name_col BED column used as feature name. Default 4.
#' @param type Feature type assigned to BED records.
#' @param verbose Whether to print progress messages.
#' @param progress Whether to show a stage progress bar.
#' @return A FeatureTrack object.
#' @examples
#' bed_file <- system.file("extdata", "gtr_demo_features.bed", package = "GeneTrackR")
#' features <- read_bed(bed_file, verbose = FALSE, progress = FALSE)
#' features
#' plot_feature_track(features, chrom = "chr1", start = 12338201, end = 12374500)
#' @export
read_bed <- function(file,
                     coordinate = c("bed", "granges"),
                     name_col = 4L,
                     type = "BED",
                     verbose = TRUE,
                     progress = interactive() && isTRUE(verbose)) {
  coordinate <- match.arg(coordinate)
  verbose <- isTRUE(verbose)
  progress <- isTRUE(progress)
  old_dt_options <- options(datatable.verbose = FALSE)
  on.exit(options(old_dt_options), add = TRUE)
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))

  progress_msg <- make_progress_message(verbose)
  stage <- make_stage_progress(total = 3L, progress = progress)
  on.exit(stage$close(), add = TRUE)
  input_file <- normalizePath(file, winslash = "/", mustWork = FALSE)

  progress_msg("%s", format_file_size_message("BED", file))
  dt <- data.table::fread(file, header = FALSE, sep = "\t", data.table = TRUE, showProgress = FALSE)
  stage$tick()
  original_ncol <- ncol(dt)
  stop_if_not(original_ncol >= 3L, "A BED file must contain at least 3 columns.")
  progress_msg("[GeneTrackR] Loaded %s records with %s columns.", format(nrow(dt), big.mark = ","), ncol(dt))
  progress_msg("[GeneTrackR] Detected input format: BED")

  old_names <- names(dt)
  data.table::setnames(dt, old_names[seq_len(3L)], c("chrom", "bed_start", "bed_end"))
  dt[, "start" := if (coordinate == "bed") as.integer(dt[["bed_start"]]) + 1L else as.integer(dt[["bed_start"]])]
  dt[, "end" := as.integer(dt[["bed_end"]])]
  dt[, "feature_id" := paste0("BED_", seq_len(.N))]
  if (original_ncol >= name_col) {
    dt[, "name" := as.character(dt[[paste0("V", name_col)]])]
  } else {
    dt[, "name" := as.character(dt[["feature_id"]])]
  }
  dt[, "type" := as.character(type)]
  dt[, "level" := "feature"]
  dt[, "score" := if (original_ncol >= 5L) suppressWarnings(as.numeric(dt[["V5"]])) else NA_real_]
  dt[, "strand" := if (original_ncol >= 6L) as.character(dt[["V6"]]) else "*"]
  dt[, "source" := "BED"]
  dt[, "parent_id" := NA_character_]
  dt[, "gene_id" := NA_character_]
  dt[, "transcript_id" := NA_character_]
  stage$tick()

  out <- dt[, .(feature_id, name, chrom, start, end, type, level, score, strand, source, gene_id, transcript_id, parent_id)]
  progress_msg("[GeneTrackR] Building Feature object and derived gene/transcript/exon tables.")
  obj <- FeatureTrack(
    out,
    meta = list(
      source_file = input_file,
      format = "BED",
      coordinate_input = coordinate,
      coordinate_internal = "1-based closed"
    )
  )
  stage$tick()

  progress_msg(
    "[GeneTrackR] Finished. records: %s; genes: %s; transcripts: %s; exons: %s.",
    format(nrow(obj$data), big.mark = ","),
    format(if (!is.null(obj$genes)) nrow(obj$genes) else 0L, big.mark = ","),
    format(if (!is.null(obj$transcripts)) nrow(obj$transcripts) else 0L, big.mark = ","),
    format(if (!is.null(obj$exons)) nrow(obj$exons) else 0L, big.mark = ",")
  )
  obj
}
