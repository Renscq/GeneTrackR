# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.7
# Function: Write unified Feature annotation objects to GenePred, GFF, GTF, or BED files
# Input: Feature/GenePred-compatible annotation object
# Output: Annotation files

#' Write Feature annotation data
#'
#' @description
#' Writes a unified `Feature` annotation object to GenePred, GenePredExt, GFF3,
#' GTF, or BED format. This is the main output function for annotation tracks
#' produced by `read_genepred()`, `read_gff()`, `read_gtf()`, or `read_bed()`.
#'
#' @param object A `Feature`, `FeatureTrack`, or GenePred-compatible object.
#' @param file Output file path.
#' @param format Output format. Use `auto`, `genepred`, `genepredext`, `gff`, `gtf`, or `bed`.
#' @param coordinate Output coordinate system for GenePred/BED-like formats.
#' `ucsc` means 0-based half-open where applicable; `granges` means 1-based closed.
#' @param overwrite Whether to overwrite an existing file.
#' @param keep_attributes Whether to reuse existing `attribute` strings for GFF/GTF output when available.
#' @return Invisibly returns the output file path.
#' @examples
#' \dontrun{
#' write_feature(gp, "annotation.genePredExt", format = "genepredext")
#' write_feature(gtf, "annotation.gtf", format = "gtf")
#' write_feature(gff, "annotation.gff3", format = "gff")
#' write_feature(bed, "regions.bed", format = "bed")
#' }
#' @export
write_feature <- function(object,
                          file,
                          format = c("auto", "genepred", "genepredext", "gff", "gtf", "bed"),
                          coordinate = c("ucsc", "granges"),
                          overwrite = FALSE,
                          keep_attributes = TRUE) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"), "`object` must be a Feature/FeatureTrack or GenePred-compatible object.")
  format <- match.arg(format)
  coordinate <- match.arg(coordinate)
  if (format == "auto") format <- infer_write_feature_format(file)
  check_output_file(file, overwrite)

  if (format %in% c("genepred", "genepredext")) {
    gp <- as_genepred(object)
    write_feature_as_genepred(gp, file, format = format, coordinate = coordinate)
  } else if (format == "bed") {
    write_feature_as_bed(object, file, coordinate = coordinate)
  } else if (format == "gff") {
    write_feature_as_gff(object, file, keep_attributes = keep_attributes)
  } else if (format == "gtf") {
    write_feature_as_gtf(object, file, keep_attributes = keep_attributes)
  }
  invisible(file)
}

infer_write_feature_format <- function(file) {
  x <- tolower(basename(file))
  x <- sub("\\.(gz|bgz|bz2|xz)$", "", x)
  if (grepl("\\.genepredext$", x)) return("genepredext")
  if (grepl("\\.genepred$", x)) return("genepred")
  if (grepl("\\.(gff|gff3)$", x)) return("gff")
  if (grepl("\\.gtf$", x)) return("gtf")
  if (grepl("\\.bed$", x)) return("bed")
  stop("Cannot infer output format from file extension. Please set `format` explicitly.", call. = FALSE)
}

write_feature_as_genepred <- function(object, file, format = c("genepred", "genepredext"), coordinate = c("ucsc", "granges")) {
  format <- match.arg(format)
  coordinate <- match.arg(coordinate)
  tx <- data.table::copy(object$transcripts)
  ex <- data.table::copy(object$exons)
  stop_if_not(nrow(tx) > 0L && nrow(ex) > 0L, "GenePred output requires transcript and exon tables.")
  data.table::setorder(ex, transcript_id, exon_start, exon_end)

  if (!"exon_frame" %in% names(ex)) ex[, "exon_frame" := -1L]
  exon_agg <- ex[, .(
    exonStarts = paste_comma_integer(if (coordinate == "ucsc") as.integer(exon_start) - 1L else as.integer(exon_start)),
    exonEnds = paste_comma_integer(as.integer(exon_end)),
    exonFrames = paste_comma_integer(as.integer(exon_frame))
  ), by = transcript_id]

  out <- merge(tx, exon_agg, by = "transcript_id", all.x = TRUE)
  if (!"score" %in% names(out)) out[, "score" := 0]
  if (!"cds_start_stat" %in% names(out)) out[, "cds_start_stat" := "unk"]
  if (!"cds_end_stat" %in% names(out)) out[, "cds_end_stat" := "unk"]

  out[, `:=`(
    name = as.character(transcript_id),
    txStart = if (coordinate == "ucsc") as.integer(tx_start) - 1L else as.integer(tx_start),
    txEnd = as.integer(tx_end),
    cdsStart = if (coordinate == "ucsc") as.integer(cds_start) - 1L else as.integer(cds_start),
    cdsEnd = as.integer(cds_end),
    exonCount = as.integer(exon_count),
    name2 = as.character(gene_id),
    score = ifelse(is.na(score), 0, score),
    cdsStartStat = ifelse(is.na(cds_start_stat) | cds_start_stat == "", "unk", cds_start_stat),
    cdsEndStat = ifelse(is.na(cds_end_stat) | cds_end_stat == "", "unk", cds_end_stat)
  )]

  standard_cols <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds")
  if (format == "genepred") {
    data.table::fwrite(out[, ..standard_cols], file, sep = "\t", col.names = FALSE)
  } else {
    ext_cols <- c(standard_cols, "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
    data.table::fwrite(out[, ..ext_cols], file, sep = "\t", col.names = FALSE)
  }
}

write_feature_as_bed <- function(object, file, coordinate = c("ucsc", "granges")) {
  coordinate <- match.arg(coordinate)
  dt <- as_feature_table(object)
  out <- dt[, .(
    chrom = as.character(chrom),
    chromStart = if (coordinate == "ucsc") as.integer(start) - 1L else as.integer(start),
    chromEnd = as.integer(end),
    name = as.character(name),
    score = ifelse(is.na(score), 0, score),
    strand = as.character(strand)
  )]
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE)
}

write_feature_as_gff <- function(object, file, keep_attributes = TRUE) {
  dt <- as_feature_table(object)
  attr <- build_gff_attributes(dt, keep_attributes = keep_attributes)
  out <- dt[, .(
    seqid = as.character(chrom),
    source = ifelse(is.na(source) | source == "", ".", as.character(source)),
    type = as.character(type),
    start = as.integer(start),
    end = as.integer(end),
    score = ifelse(is.na(score), ".", as.character(score)),
    strand = ifelse(is.na(strand) | strand == "*", ".", as.character(strand)),
    phase = ifelse(is.na(phase) | phase == "", ".", as.character(phase)),
    attributes = attr
  )]
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE)
}

write_feature_as_gtf <- function(object, file, keep_attributes = TRUE) {
  dt <- as_feature_table(object)
  attr <- build_gtf_attributes(dt, keep_attributes = keep_attributes)
  out <- dt[, .(
    seqname = as.character(chrom),
    source = ifelse(is.na(source) | source == "", ".", as.character(source)),
    feature = normalize_gtf_feature_type(type),
    start = as.integer(start),
    end = as.integer(end),
    score = ifelse(is.na(score), ".", as.character(score)),
    strand = ifelse(is.na(strand) | strand == "*", ".", as.character(strand)),
    frame = ifelse(is.na(phase) | phase == "", ".", as.character(phase)),
    attribute = attr
  )]
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE)
}

build_gff_attributes <- function(dt, keep_attributes = TRUE) {
  if (isTRUE(keep_attributes) && "attribute" %in% names(dt)) {
    attr <- as.character(dt$attribute)
    has_attr <- !is.na(attr) & nzchar(attr)
  } else {
    attr <- rep(NA_character_, nrow(dt))
    has_attr <- rep(FALSE, nrow(dt))
  }
  out <- attr
  idx <- which(!has_attr)
  if (length(idx) > 0L) {
    out[idx] <- vapply(idx, function(i) {
      vals <- c(
        ID = dt$feature_id[i],
        Name = dt$name[i],
        Parent = dt$parent_id[i],
        gene_id = dt$gene_id[i],
        transcript_id = dt$transcript_id[i],
        gene_type = dt$gene_type[i]
      )
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0L) return(".")
      paste(paste0(names(vals), "=", vals), collapse = ";")
    }, character(1L))
  }
  out
}

build_gtf_attributes <- function(dt, keep_attributes = TRUE) {
  if (isTRUE(keep_attributes) && "attribute" %in% names(dt)) {
    attr <- as.character(dt$attribute)
    has_attr <- !is.na(attr) & nzchar(attr)
  } else {
    attr <- rep(NA_character_, nrow(dt))
    has_attr <- rep(FALSE, nrow(dt))
  }
  out <- attr
  idx <- which(!has_attr)
  if (length(idx) > 0L) {
    out[idx] <- vapply(idx, function(i) {
      vals <- c(
        gene_id = dt$gene_id[i],
        transcript_id = dt$transcript_id[i],
        gene_name = dt$name[i],
        gene_type = dt$gene_type[i],
        exon_number = if (!is.na(dt$exon_number[i])) as.character(dt$exon_number[i]) else NA_character_
      )
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0L) return(".")
      paste0(paste(paste0(names(vals), " \"", vals, "\""), collapse = "; "), ";")
    }, character(1L))
  }
  out
}

normalize_gtf_feature_type <- function(x) {
  x <- as.character(x)
  x[tolower(x) == "mrna"] <- "transcript"
  x
}

#' Write a GenePred-compatible annotation object
#'
#' @description Backward-compatible wrapper around `write_feature()`.
#' @inheritParams write_feature
#' @export
write_genepred <- function(object, file, format = c("genePred", "genePredExt"), coordinate = c("ucsc", "granges")) {
  format <- match.arg(format)
  fmt <- if (format == "genePred") "genepred" else "genepredext"
  write_feature(object = object, file = file, format = fmt, coordinate = coordinate)
}

#' Write a FeatureTrack object
#'
#' @description Backward-compatible wrapper around `write_feature()`.
#' @inheritParams write_feature
#' @export
write_feature_track <- function(object, file, format = c("bed", "gff", "gtf", "genepred", "genepredext")) {
  format <- match.arg(format)
  write_feature(object = object, file = file, format = format)
}
