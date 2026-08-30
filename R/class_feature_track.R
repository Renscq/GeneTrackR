# Author: Rensc
# Date: 2026-08-31
# Version: dev002
# Function: Define unified Feature annotation class and GenePred-compatible constructor
# Input: Standardized annotation or variant tables
# Output: Unified S3 track objects

#' Create a standardized Feature object
#'
#' @description
#' `Feature()` is the unified annotation class used by GeneTrackR. BED, GFF,
#' GTF, and GenePred annotations are all normalized to a `Feature`-compatible
#' object with a flat interval table in `object$data`. GenePred-compatible
#' inputs additionally store derived `genes`, `transcripts`, and `exons` tables,
#' allowing the same object to be used by gene model plotting functions.
#'
#' Required coordinate convention is 1-based closed intervals.
#'
#' @param data A data.frame or data.table containing standardized interval columns.
#' Required columns are `chrom`, `start`, and `end`. Recommended columns include
#' `feature_id`, `name`, `type`, `score`, `strand`, `source`, `gene_id`,
#' `transcript_id`, `parent_id`, `level`, and `gene_type`.
#' @param meta Metadata list such as source file, original format, and coordinate system.
#' @param genes Optional gene-level table.
#' @param transcripts Optional transcript-level table.
#' @param exons Optional exon-level table.
#' @param validation Optional validation list.
#' @return A Feature object, also inheriting from FeatureTrack.
#' @export
Feature <- function(data,
                    meta = list(),
                    genes = NULL,
                    transcripts = NULL,
                    exons = NULL,
                    validation = make_empty_validation()) {
  dt <- standardize_feature_table(data)
  hierarchy <- derive_feature_hierarchy(dt, genes = genes, transcripts = transcripts, exons = exons)
  structure(
    list(
      data = dt,
      genes = hierarchy$genes,
      transcripts = hierarchy$transcripts,
      exons = hierarchy$exons,
      meta = meta,
      validation = validation
    ),
    class = c("Feature", "FeatureTrack")
  )
}

#' Create a FeatureTrack object
#'
#' @description
#' `FeatureTrack()` is kept as a user-facing alias of `Feature()` for interval
#' tracks. It returns an object inheriting from both `Feature` and `FeatureTrack`.
#'
#' @inheritParams Feature
#' @return A Feature object inheriting from FeatureTrack.
#' @export
FeatureTrack <- function(data, meta = list(), genes = NULL, transcripts = NULL, exons = NULL, validation = make_empty_validation()) {
  Feature(data = data, meta = meta, genes = genes, transcripts = transcripts, exons = exons, validation = validation)
}

standardize_feature_table <- function(data) {
  dt <- data.table::as.data.table(data)
  required_cols <- c("chrom", "start", "end")
  stop_if_not(all(required_cols %in% names(dt)), "`data` must contain `chrom`, `start`, and `end` columns.")

  if (!"feature_id" %in% names(dt)) dt[, "feature_id" := paste0("feature_", seq_len(.N))]
  if (!"name" %in% names(dt)) dt[, "name" := as.character(dt[["feature_id"]])]
  if (!"type" %in% names(dt)) dt[, "type" := "feature"]
  if (!"level" %in% names(dt)) dt[, "level" := infer_feature_level(dt[["type"]])]
  if (!"score" %in% names(dt)) dt[, "score" := NA_real_]
  if (!"strand" %in% names(dt)) dt[, "strand" := "*"]
  if (!"source" %in% names(dt)) dt[, "source" := NA_character_]
  if (!"gene_id" %in% names(dt)) dt[, "gene_id" := NA_character_]
  if (!"transcript_id" %in% names(dt)) dt[, "transcript_id" := NA_character_]
  if (!"parent_id" %in% names(dt)) dt[, "parent_id" := NA_character_]
  if (!"gene_type" %in% names(dt)) dt[, "gene_type" := NA_character_]
  if (!"exon_number" %in% names(dt)) dt[, "exon_number" := NA_integer_]
  if (!"phase" %in% names(dt)) dt[, "phase" := NA_character_]
  if (!"attribute" %in% names(dt)) dt[, "attribute" := NA_character_]

  dt[, "chrom" := as.character(dt[["chrom"]])]
  dt[, "start" := as.integer(dt[["start"]])]
  dt[, "end" := as.integer(dt[["end"]])]
  dt[, "feature_id" := as.character(dt[["feature_id"]])]
  dt[, "name" := as.character(dt[["name"]])]
  dt[, "type" := as.character(dt[["type"]])]
  dt[, "level" := as.character(dt[["level"]])]
  dt[, "score" := suppressWarnings(as.numeric(dt[["score"]]))]
  dt[, "strand" := as.character(dt[["strand"]])]
  dt[is.na(strand) | !strand %in% c("+", "-", "*", "."), "strand" := "*"]
  dt[strand == ".", "strand" := "*"]
  dt[, "source" := as.character(dt[["source"]])]
  dt[, "gene_id" := as.character(dt[["gene_id"]])]
  dt[, "transcript_id" := as.character(dt[["transcript_id"]])]
  dt[, "parent_id" := as.character(dt[["parent_id"]])]
  dt[, "gene_type" := as.character(dt[["gene_type"]])]
  dt[, "exon_number" := suppressWarnings(as.integer(dt[["exon_number"]]))]
  dt[, "phase" := as.character(dt[["phase"]])]
  dt[, "attribute" := as.character(dt[["attribute"]])]

  dt[is.na(gene_id) | gene_id == "", "gene_id" := NA_character_]
  dt[is.na(transcript_id) | transcript_id == "", "transcript_id" := NA_character_]
  dt[is.na(parent_id) | parent_id == "", "parent_id" := NA_character_]
  dt[is.na(gene_type) | gene_type == "", "gene_type" := NA_character_]

  data.table::setorderv(dt, c("chrom", "start", "end", "feature_id"))
  dt[]
}

infer_feature_level <- function(type) {
  type <- tolower(as.character(type))
  out <- rep("feature", length(type))
  out[type %in% c("gene")] <- "gene"
  out[type %in% c("mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna", "mirna", "primary_transcript")] <- "transcript"
  out[type %in% c("exon")] <- "exon"
  out[type %in% c("cds", "utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr", "five_utr", "three_utr")] <- "subfeature"
  out
}

is_transcript_feature_type <- function(x) {
  tolower(as.character(x)) %in% c("mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna", "mirna", "primary_transcript")
}

is_gene_model_feature <- function(object) {
  inherits(object, "GenePred") || (
    inherits(object, "Feature") &&
      !is.null(object$transcripts) && nrow(object$transcripts) > 0L &&
      !is.null(object$exons) && nrow(object$exons) > 0L
  )
}

derive_feature_hierarchy <- function(dt, genes = NULL, transcripts = NULL, exons = NULL) {
  if (!is.null(genes) && !is.null(transcripts) && !is.null(exons)) {
    return(list(
      genes = data.table::as.data.table(genes),
      transcripts = data.table::as.data.table(transcripts),
      exons = data.table::as.data.table(exons)
    ))
  }

  x <- data.table::copy(data.table::as.data.table(dt))
  gene_rows <- x[tolower(x[["level"]]) == "gene" | tolower(x[["type"]]) == "gene"]
  tx_rows <- x[tolower(x[["level"]]) == "transcript" | is_transcript_feature_type(x[["type"]])]
  exon_rows <- x[tolower(x[["level"]]) == "exon" | tolower(x[["type"]]) == "exon"]
  cds_rows <- x[tolower(x[["type"]]) == "cds"]

  if (is.null(transcripts)) {
    if (nrow(tx_rows) > 0L) {
      transcripts <- data.table::copy(tx_rows)
      transcripts[is.na(transcript_id) | transcript_id == "", "transcript_id" := feature_id]
      transcripts[is.na(gene_id) | gene_id == "", "gene_id" := parent_id]
      transcripts[is.na(gene_id) | gene_id == "", "gene_id" := transcript_id]
      transcripts <- transcripts[, .(
        transcript_id = as.character(transcript_id),
        gene_id = as.character(gene_id),
        chrom = as.character(chrom),
        strand = as.character(strand),
        tx_start = as.integer(start),
        tx_end = as.integer(end),
        gene_type = as.character(gene_type)
      )]
    } else if (nrow(exon_rows) > 0L && any(!is.na(exon_rows[["transcript_id"]]))) {
      transcripts <- exon_rows[!is.na(transcript_id), .(
        chrom = chrom[1L],
        strand = strand[1L],
        tx_start = as.integer(min(start, na.rm = TRUE)),
        tx_end = as.integer(max(end, na.rm = TRUE)),
        gene_id = pick_first_nonempty(gene_id[1L], parent_id[1L], transcript_id[1L]),
        gene_type = NA_character_
      ), by = transcript_id]
    } else {
      transcripts <- data.table::data.table()
    }
  } else {
    transcripts <- data.table::as.data.table(transcripts)
  }

  if (nrow(transcripts) > 0L) {
    cds_summary <- data.table::data.table(transcript_id = character(), cds_start = integer(), cds_end = integer())
    if (nrow(cds_rows) > 0L) {
      cds_rows[is.na(transcript_id) | transcript_id == "", "transcript_id" := parent_id]
      cds_summary <- cds_rows[!is.na(transcript_id), .(
        cds_start = as.integer(min(start, na.rm = TRUE)),
        cds_end = as.integer(max(end, na.rm = TRUE))
      ), by = transcript_id]
    }
    transcripts <- merge(transcripts, cds_summary, by = "transcript_id", all.x = TRUE)
    transcripts[is.na(cds_start) | is.na(cds_end), `:=`(cds_start = as.integer(tx_start), cds_end = as.integer(tx_start - 1L))]
    transcripts[, "gene_type" := data.table::fifelse(as.integer(cds_start) <= as.integer(cds_end), "coding", "non-coding")]
    transcripts[, "exon_count" := NA_integer_]
    transcripts[, "score" := NA_real_]
    transcripts[, "cds_start_stat" := NA_character_]
    transcripts[, "cds_end_stat" := NA_character_]
    transcripts[, "exon_frames" := NA_character_]
  }

  if (is.null(exons)) {
    if (nrow(exon_rows) > 0L) {
      exons <- data.table::copy(exon_rows)
      exons[is.na(transcript_id) | transcript_id == "", "transcript_id" := parent_id]
      if (nrow(transcripts) > 0L) {
        tx_map <- transcripts[, .(transcript_id, gene_id)]
        exons <- merge(exons, tx_map, by = "transcript_id", all.x = TRUE, suffixes = c("", ".tx"))
        exons[is.na(gene_id) | gene_id == "", "gene_id" := gene_id.tx]
        exons[, "gene_id.tx" := NULL]
      }
      exons[is.na(gene_id) | gene_id == "", "gene_id" := parent_id]
      data.table::setorderv(exons, c("transcript_id", "start", "end"))
      exons[, "exon_number" := seq_len(.N), by = transcript_id]
      exons <- exons[, .(
        transcript_id = as.character(transcript_id),
        gene_id = as.character(gene_id),
        chrom = as.character(chrom),
        strand = as.character(strand),
        exon_number = as.integer(exon_number),
        exon_start = as.integer(start),
        exon_end = as.integer(end),
        exon_frame = NA_integer_
      )]
    } else {
      exons <- data.table::data.table()
    }
  } else {
    exons <- data.table::as.data.table(exons)
  }

  if (nrow(transcripts) > 0L && nrow(exons) > 0L) {
    exon_counts <- exons[, .(exon_count = as.integer(.N)), by = transcript_id]
    if ("exon_count" %in% names(transcripts)) {
      transcripts[, "exon_count" := NULL]
    }
    transcripts <- merge(transcripts, exon_counts, by = "transcript_id", all.x = TRUE)
    transcripts[is.na(exon_count), "exon_count" := 0L]
  }

  frame_info <- fill_cds_status_and_exon_frames(transcripts, exons, cds_rows)
  transcripts <- frame_info$transcripts
  exons <- frame_info$exons

  if (is.null(genes)) {
    if (nrow(transcripts) > 0L) {
      genes <- build_gene_table(transcripts)
    } else if (nrow(gene_rows) > 0L) {
      genes <- gene_rows[, .(
        gene_id = pick_first_nonempty(gene_id, feature_id, name),
        chrom = chrom,
        strand = strand,
        gene_start = as.integer(start),
        gene_end = as.integer(end),
        n_transcripts = 0L,
        gene_type = data.table::fifelse(!is.na(gene_type) & gene_type != "", gene_type, "unknown")
      )]
    } else {
      genes <- data.table::data.table()
    }
  } else {
    genes <- data.table::as.data.table(genes)
  }

  list(genes = genes, transcripts = transcripts, exons = exons)
}

#' Coerce an annotation object to a unified Feature object
#'
#' @param object A GenePred, FeatureTrack, or data.frame-like object.
#' @return A Feature object.
#' @export
as_feature <- function(object) {
  UseMethod("as_feature")
}

#' @export
as_feature.Feature <- function(object) object

#' @export
as_feature.FeatureTrack <- function(object) object

#' @export
as_feature.GenePred <- function(object) {
  Feature(
    data = object$data %||% genepred_to_feature_table(object$transcripts, object$exons, object$genes),
    genes = object$genes,
    transcripts = object$transcripts,
    exons = object$exons,
    meta = object$meta,
    validation = object$validation
  )
}

#' @export
as_feature.data.frame <- function(object) Feature(object)

#' Extract the standardized feature table
#'
#' @param object A Feature, FeatureTrack, or GenePred object.
#' @return A data.table with standardized feature columns.
#' @export
as_feature_table <- function(object) {
  if (inherits(object, "GenePred") && !is.null(object$data)) return(data.table::copy(object$data))
  if (inherits(object, "Feature") || inherits(object, "FeatureTrack")) return(data.table::copy(object$data))
  stop("`object` must be a Feature, FeatureTrack, or GenePred object.", call. = FALSE)
}

#' Convert a Feature object with gene-model information to GenePred
#'
#' @param object A Feature object containing transcript and exon features.
#' @return A GenePred object.
#' @export
as_genepred <- function(object) {
  if (inherits(object, "GenePred")) return(object)
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack"), "`object` must be a Feature/FeatureTrack or GenePred object.")
  stop_if_not(!is.null(object$transcripts) && nrow(object$transcripts) > 0L, "The Feature object does not contain transcript-level gene model information.")
  stop_if_not(!is.null(object$exons) && nrow(object$exons) > 0L, "The Feature object does not contain exon-level gene model information.")

  # Fast path: GTF/GFF/BED-derived Feature objects already store standardized
  # gene/transcript/exon tables. Rebuilding a GenePred object would regenerate
  # the full feature table and can be slow for large annotations. Adding the
  # GenePred-compatible class is sufficient for downstream accessors and plots.
  cls <- unique(c(class(object), "GenePred"))
  class(object) <- cls
  object
}

#' @export
print.Feature <- function(x, ...) {
  cat("<Feature>\n")
  cat("  records    : ", format(nrow(x$data), big.mark = ","), "\n", sep = "")
  cat("  genes      : ", if (!is.null(x$genes)) nrow(x$genes) else 0L, "\n", sep = "")
  cat("  transcripts: ", if (!is.null(x$transcripts)) nrow(x$transcripts) else 0L, "\n", sep = "")
  cat("  exons      : ", if (!is.null(x$exons)) nrow(x$exons) else 0L, "\n", sep = "")
  cat("  format     : ", x$meta$format %||% "unknown", "\n", sep = "")
  cat("  coordinate : ", x$meta$coordinate_internal %||% "1-based closed", "\n", sep = "")
  invisible(x)
}

#' @export
print.FeatureTrack <- function(x, ...) {
  print.Feature(x, ...)
}


####################################################################
# GenePred-compatible constructor and coercion helpers
####################################################################

#' Create a GenePred object
#'
#' @param transcripts Transcript-level annotation table.
#' @param exons Exon-level annotation table.
#' @param genes Gene-level annotation table.
#' @param meta Metadata list.
#' @param validation Validation result list.
#' @return A GenePred object.
#' @export
GenePred <- function(transcripts, exons, genes = NULL, meta = list(), validation = make_empty_validation()) {
  if (is.null(genes)) {
    genes <- build_gene_table(transcripts)
  }
  tx <- data.table::as.data.table(transcripts)
  ex <- data.table::as.data.table(exons)
  gn <- data.table::as.data.table(genes)
  ft <- genepred_to_feature_table(tx, ex, gn)
  structure(
    list(
      data = ft,
      transcripts = tx,
      exons = ex,
      genes = gn,
      meta = meta,
      validation = validation
    ),
    class = c("Feature", "FeatureTrack", "GenePred")
  )
}

#' @export
print.GenePred <- function(x, ...) {
  cat("<GenePred>\n")
  cat("  genes      : ", nrow(x$genes), "\n", sep = "")
  cat("  transcripts: ", nrow(x$transcripts), "\n", sep = "")
  cat("  exons      : ", nrow(x$exons), "\n", sep = "")
  cat("  coordinate : ", x$meta$coordinate_internal %||% "1-based closed", "\n", sep = "")
  invisible(x)
}

#' @export
summary.GenePred <- function(object, ...) {
  summary_feature(object, ...)
}

#' @export
summary.Feature <- function(object, ...) {
  summary_feature(object, ...)
}

#' @export
summary.FeatureTrack <- function(object, ...) {
  summary_feature(object, ...)
}

#' Extract transcript table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with transcript-level annotation.
#' @export
as_transcript_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$transcripts)
}

#' Extract exon table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with exon-level annotation.
#' @export
as_exon_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$exons)
}

#' Extract gene table from a GenePred object
#'
#' @param object A GenePred object.
#' @return A data.table with gene-level annotation.
#' @export
as_gene_table <- function(object) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  data.table::copy(object$genes)
}

genepred_to_feature_table <- function(transcripts, exons, genes = NULL) {
  tx <- data.table::copy(data.table::as.data.table(transcripts))
  ex <- data.table::copy(data.table::as.data.table(exons))
  if (is.null(genes)) genes <- build_gene_table(tx)
  gn <- data.table::copy(data.table::as.data.table(genes))

  # GenePred itself does not store the original GTF/GFF source column.
  # GeneTrackR uses "ribo" for generated GTF/GFF-compatible records so that
  # round-tripped files are closer to common ribosome-profiling annotations.
  default_source <- "ribo"

  out <- list()
  if (nrow(gn) > 0L) {
    out[[length(out) + 1L]] <- data.table::data.table(
      feature_id = as.character(gn[["gene_id"]]),
      name = as.character(gn[["gene_id"]]),
      gene_name = as.character(gn[["gene_id"]]),
      exon_id = NA_character_,
      chrom = as.character(gn[["chrom"]]),
      start = as.integer(gn[["gene_start"]]),
      end = as.integer(gn[["gene_end"]]),
      type = "gene",
      level = "gene",
      score = NA_real_,
      strand = as.character(gn[["strand"]]),
      source = default_source,
      gene_id = as.character(gn[["gene_id"]]),
      transcript_id = NA_character_,
      parent_id = NA_character_,
      gene_type = as.character(gn[["gene_type"]]),
      exon_number = NA_integer_,
      phase = NA_character_,
      attribute = NA_character_
    )
  }

  if (nrow(tx) > 0L) {
    out[[length(out) + 1L]] <- data.table::data.table(
      feature_id = as.character(tx[["transcript_id"]]),
      name = as.character(tx[["transcript_id"]]),
      gene_name = as.character(tx[["gene_id"]]),
      exon_id = NA_character_,
      chrom = as.character(tx[["chrom"]]),
      start = as.integer(tx[["tx_start"]]),
      end = as.integer(tx[["tx_end"]]),
      type = "transcript",
      level = "transcript",
      score = NA_real_,
      strand = as.character(tx[["strand"]]),
      source = default_source,
      gene_id = as.character(tx[["gene_id"]]),
      transcript_id = as.character(tx[["transcript_id"]]),
      parent_id = as.character(tx[["gene_id"]]),
      gene_type = as.character(tx[["gene_type"]]),
      exon_number = NA_integer_,
      phase = NA_character_,
      attribute = NA_character_
    )
  }

  if (nrow(ex) > 0L) {
    tx_small <- tx[, c("transcript_id", "cds_start", "cds_end", "gene_type"), with = FALSE]
    ex2 <- merge(ex, tx_small, by = "transcript_id", all.x = TRUE, sort = FALSE)
    ex2[, "exon_id" := paste0(as.character(transcript_id), ".", as.integer(exon_number))]

    out[[length(out) + 1L]] <- data.table::data.table(
      feature_id = as.character(ex2[["exon_id"]]),
      name = as.character(ex2[["exon_id"]]),
      gene_name = as.character(ex2[["gene_id"]]),
      exon_id = as.character(ex2[["exon_id"]]),
      chrom = as.character(ex2[["chrom"]]),
      start = as.integer(ex2[["exon_start"]]),
      end = as.integer(ex2[["exon_end"]]),
      type = "exon",
      level = "exon",
      score = NA_real_,
      strand = as.character(ex2[["strand"]]),
      source = default_source,
      gene_id = as.character(ex2[["gene_id"]]),
      transcript_id = as.character(ex2[["transcript_id"]]),
      parent_id = as.character(ex2[["transcript_id"]]),
      gene_type = as.character(ex2[["gene_type"]]),
      exon_number = as.integer(ex2[["exon_number"]]),
      phase = NA_character_,
      attribute = NA_character_
    )

    seg <- make_genepred_gtf_segments(tx, ex)
    if (nrow(seg) > 0L) {
      seg[, "exon_id" := paste0(as.character(transcript_id), ".", as.integer(exon_number))]
      out[[length(out) + 1L]] <- data.table::data.table(
        feature_id = paste(as.character(seg[["transcript_id"]]), as.character(seg[["feature"]]), as.integer(seg[["start"]]), as.integer(seg[["end"]]), sep = ":"),
        name = as.character(seg[["feature"]]),
        gene_name = as.character(seg[["gene_id"]]),
        exon_id = as.character(seg[["exon_id"]]),
        chrom = as.character(seg[["chrom"]]),
        start = as.integer(seg[["start"]]),
        end = as.integer(seg[["end"]]),
        type = as.character(seg[["feature"]]),
        level = "subfeature",
        score = NA_real_,
        strand = as.character(seg[["strand"]]),
        source = default_source,
        gene_id = as.character(seg[["gene_id"]]),
        transcript_id = as.character(seg[["transcript_id"]]),
        parent_id = as.character(seg[["transcript_id"]]),
        gene_type = NA_character_,
        exon_number = as.integer(seg[["exon_number"]]),
        phase = as.character(seg[["phase"]]),
        attribute = NA_character_
      )
    }
  }

  if (length(out) == 0L) {
    return(data.table::data.table(
      feature_id = character(), name = character(), gene_name = character(), exon_id = character(),
      chrom = character(), start = integer(), end = integer(), type = character(), level = character(),
      score = numeric(), strand = character(), source = character(), gene_id = character(),
      transcript_id = character(), parent_id = character(), gene_type = character(),
      exon_number = integer(), phase = character(), attribute = character()
    ))
  }
  standardize_feature_table(data.table::rbindlist(out, fill = TRUE))
}

make_genepred_gtf_segments <- function(tx, exons) {
  tx <- data.table::copy(data.table::as.data.table(tx))
  ex <- data.table::copy(data.table::as.data.table(exons))
  if (nrow(tx) == 0L || nrow(ex) == 0L) return(data.table::data.table())

  # Older objects may not have gene_type. Infer it from CDS coordinates when needed.
  if (!"gene_type" %in% names(tx)) {
    tx[, "gene_type" := data.table::fifelse(
      !is.na(as.integer(cds_start)) & !is.na(as.integer(cds_end)) &
        as.integer(cds_start) <= as.integer(cds_end),
      "coding",
      "non-coding"
    )]
  }

  tx_small <- tx[, .(
    transcript_id = as.character(.SD[["transcript_id"]]),
    cds_start = as.integer(.SD[["cds_start"]]),
    cds_end = as.integer(.SD[["cds_end"]]),
    gene_type = as.character(.SD[["gene_type"]])
  ), .SDcols = c("transcript_id", "cds_start", "cds_end", "gene_type")]

  # Avoid gene_type.x/gene_type.y after merge when exons already carry gene_type.
  if ("gene_type" %in% names(ex)) {
    ex[, "gene_type" := NULL]
  }
  ex <- merge(ex, tx_small, by = "transcript_id", all.x = TRUE, sort = FALSE)

  ex[, `:=`(
    exon_start = as.integer(.SD[["exon_start"]]),
    exon_end = as.integer(.SD[["exon_end"]]),
    cds_start = as.integer(.SD[["cds_start"]]),
    cds_end = as.integer(.SD[["cds_end"]]),
    strand = as.character(.SD[["strand"]]),
    gene_type = as.character(.SD[["gene_type"]])
  ), .SDcols = c("exon_start", "exon_end", "cds_start", "cds_end", "strand", "gene_type")]

  coding_key <- !is.na(ex[["gene_type"]]) & ex[["gene_type"]] == "coding" &
    !is.na(ex[["cds_start"]]) & !is.na(ex[["cds_end"]]) &
    as.integer(ex[["cds_start"]]) <= as.integer(ex[["cds_end"]])

  coding <- ex[coding_key]
  noncoding <- ex[!coding_key]
  base_cols <- c("transcript_id", "gene_id", "chrom", "strand", "exon_number")

  make_piece <- function(dt, feature, start, end, phase = NA_character_) {
    if (nrow(dt) == 0L) return(data.table::data.table())
    out <- dt[, base_cols, with = FALSE]
    out[, `:=`(
      feature = as.character(feature),
      start = as.integer(start),
      end = as.integer(end),
      phase = as.character(phase)
    )]
    out[!is.na(start) & !is.na(end) & start <= end]
  }

  out <- list()
  if (nrow(noncoding) > 0L) {
    out[[length(out) + 1L]] <- make_piece(noncoding, "exon", noncoding[["exon_start"]], noncoding[["exon_end"]])
  }

  if (nrow(coding) == 0L) {
    ans <- data.table::rbindlist(out, fill = TRUE)
    if (is.null(ans)) return(data.table::data.table())
    return(ans[])
  }

  # UTR intervals. The left genomic side is 5UTR on the plus strand but 3UTR on
  # the minus strand; the right genomic side is reversed.
  left_utr_type <- data.table::fifelse(coding[["strand"]] == "-", "3UTR", "5UTR")
  right_utr_type <- data.table::fifelse(coding[["strand"]] == "-", "5UTR", "3UTR")
  out[[length(out) + 1L]] <- make_piece(
    coding,
    left_utr_type,
    coding[["exon_start"]],
    pmin(coding[["exon_end"]], coding[["cds_start"]] - 1L)
  )
  out[[length(out) + 1L]] <- make_piece(
    coding,
    right_utr_type,
    pmax(coding[["exon_start"]], coding[["cds_end"]] + 1L),
    coding[["exon_end"]]
  )

  # GTF convention: stop codons are separate records, so CDS intervals exclude
  # the terminal stop codon. GenePred cds_start/cds_end cover the ORF span.
  tx_coding_key <- !is.na(tx[["gene_type"]]) & tx[["gene_type"]] == "coding" &
    !is.na(tx[["cds_start"]]) & !is.na(tx[["cds_end"]]) &
    as.integer(tx[["cds_start"]]) <= as.integer(tx[["cds_end"]])

  cds_tx <- unique(tx[tx_coding_key, .(
    transcript_id = as.character(.SD[["transcript_id"]]),
    cds_start = as.integer(.SD[["cds_start"]]),
    cds_end = as.integer(.SD[["cds_end"]]),
    strand = as.character(.SD[["strand"]])
  ), .SDcols = c("transcript_id", "cds_start", "cds_end", "strand")])

  if (nrow(cds_tx) > 0L) {
    cds_tx[, `:=`(
      cds_body_start = data.table::fifelse(strand == "-", cds_start + 3L, cds_start),
      cds_body_end = data.table::fifelse(strand == "-", cds_end, cds_end - 3L),
      start_codon_start = data.table::fifelse(strand == "-", cds_end - 2L, cds_start),
      start_codon_end = data.table::fifelse(strand == "-", cds_end, cds_start + 2L),
      stop_codon_start = data.table::fifelse(strand == "-", cds_start, cds_end - 2L),
      stop_codon_end = data.table::fifelse(strand == "-", cds_start + 2L, cds_end)
    )]
  }

  ex_base <- coding[, c(base_cols, "exon_start", "exon_end"), with = FALSE]
  add_tx_interval <- function(interval_dt, start_col, end_col, feature) {
    if (nrow(interval_dt) == 0L) return(data.table::data.table())
    iv <- interval_dt[, .(
      transcript_id = as.character(.SD[["transcript_id"]]),
      tx_feature_start = as.integer(.SD[[start_col]]),
      tx_feature_end = as.integer(.SD[[end_col]])
    ), .SDcols = c("transcript_id", start_col, end_col)]
    x <- merge(ex_base, iv, by = "transcript_id", all.x = FALSE, allow.cartesian = TRUE, sort = FALSE)
    x[, `:=`(
      start = pmax(as.integer(exon_start), as.integer(tx_feature_start)),
      end = pmin(as.integer(exon_end), as.integer(tx_feature_end)),
      feature = feature,
      phase = NA_character_
    )]
    x <- x[!is.na(start) & !is.na(end) & start <= end]
    x[, c(base_cols, "feature", "start", "end", "phase"), with = FALSE]
  }

  cds_seg <- add_tx_interval(cds_tx, "cds_body_start", "cds_body_end", "CDS")
  if (nrow(cds_seg) > 0L) {
    cds_seg[, "cds_len" := as.integer(end - start + 1L)]
    cds_seg[, "order_key" := data.table::fifelse(as.character(strand) == "-", -as.integer(end), as.integer(start))]
    data.table::setorderv(cds_seg, c("transcript_id", "order_key", "start", "end"))
    cds_seg[, "cum_before" := cumsum(data.table::shift(cds_len, fill = 0L)), by = transcript_id]
    cds_seg[, "phase" := as.character((3L - (as.integer(cum_before) %% 3L)) %% 3L)]
    cds_seg[, c("cds_len", "order_key", "cum_before") := NULL]
    out[[length(out) + 1L]] <- cds_seg
  }

  out[[length(out) + 1L]] <- add_tx_interval(cds_tx, "start_codon_start", "start_codon_end", "start_codon")
  out[[length(out) + 1L]] <- add_tx_interval(cds_tx, "stop_codon_start", "stop_codon_end", "stop_codon")

  ans <- data.table::rbindlist(out, fill = TRUE)
  if (nrow(ans) == 0L) return(ans)
  ans[]
}

build_gene_table <- function(transcripts) {
  dt <- data.table::as.data.table(transcripts)
  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      gene_id = character(), chrom = character(), strand = character(),
      gene_start = integer(), gene_end = integer(), n_transcripts = integer(),
      gene_type = character()
    ))
  }
  dt[, .(
    chrom = chrom[1],
    strand = strand[1],
    gene_start = min(tx_start),
    gene_end = max(tx_end),
    n_transcripts = data.table::uniqueN(transcript_id),
    gene_type = if (any(gene_type == "coding")) "coding" else "non-coding"
  ), by = gene_id][order(chrom, gene_start, gene_end, gene_id)]
}


fill_cds_status_and_exon_frames <- function(transcripts, exons, cds_rows) {
  tx <- data.table::copy(data.table::as.data.table(transcripts))
  ex <- data.table::copy(data.table::as.data.table(exons))
  cds <- data.table::copy(data.table::as.data.table(cds_rows))

  if (nrow(tx) == 0L) {
    return(list(transcripts = tx, exons = ex))
  }

  if (!"cds_start_stat" %in% names(tx)) tx[, "cds_start_stat" := NA_character_]
  if (!"cds_end_stat" %in% names(tx)) tx[, "cds_end_stat" := NA_character_]
  if (!"exon_frames" %in% names(tx)) tx[, "exon_frames" := NA_character_]
  if (!"exon_frame" %in% names(ex)) ex[, "exon_frame" := NA_integer_]

  tx[, "cds_start_stat" := data.table::fifelse(as.character(gene_type) == "coding", "cmpl", "none")]
  tx[, "cds_end_stat" := data.table::fifelse(as.character(gene_type) == "coding", "cmpl", "none")]

  if (nrow(ex) == 0L) {
    tx[, "exon_frames" := NA_character_]
    return(list(transcripts = tx, exons = ex))
  }

  ex[, "exon_frame" := -1L]

  if (nrow(cds) > 0L) {
    if (!"transcript_id" %in% names(cds)) cds[, "transcript_id" := NA_character_]
    if (!"parent_id" %in% names(cds)) cds[, "parent_id" := NA_character_]
    cds[is.na(transcript_id) | transcript_id == "", "transcript_id" := parent_id]
    cds <- cds[!is.na(transcript_id) & transcript_id != ""]
    if (nrow(cds) > 0L) {
      cds[, "phase_int" := suppressWarnings(as.integer(phase))]
      cds[is.na(phase_int), "phase_int" := 0L]
      cds_iv <- cds[, .(
        transcript_id = as.character(transcript_id),
        start = as.integer(start),
        end = as.integer(end),
        cds_start = as.integer(start),
        cds_end = as.integer(end),
        phase_int = as.integer(phase_int)
      )]
      ex_iv <- ex[, .(
        exon_row = .I,
        transcript_id = as.character(transcript_id),
        start = as.integer(exon_start),
        end = as.integer(exon_end),
        strand = as.character(strand)
      )]
      cds_iv <- cds_iv[!is.na(start) & !is.na(end) & start <= end]
      ex_iv <- ex_iv[!is.na(start) & !is.na(end) & start <= end & !is.na(transcript_id) & transcript_id != ""]
      if (nrow(cds_iv) > 0L && nrow(ex_iv) > 0L) {
        data.table::setkey(cds_iv, transcript_id, start, end)
        ov <- data.table::foverlaps(
          x = ex_iv,
          y = cds_iv,
          by.x = c("transcript_id", "start", "end"),
          by.y = c("transcript_id", "start", "end"),
          type = "any",
          nomatch = 0L
        )
        if (nrow(ov) > 0L) {
          exon_col <- if ("i.exon_row" %in% names(ov)) "i.exon_row" else "exon_row"
          strand_col <- if ("i.strand" %in% names(ov)) "i.strand" else "strand"
          ov[, "exon_row_key" := as.integer(.SD[[exon_col]]), .SDcols = exon_col]
          ov[, "strand_key" := as.character(.SD[[strand_col]]), .SDcols = strand_col]
          ov[, "rank_pos" := data.table::fifelse(as.character(strand_key) == "-", -as.integer(cds_end), as.integer(cds_start))]
          data.table::setorderv(ov, c("exon_row_key", "rank_pos"))
          chosen <- ov[, .SD[1L], by = exon_row_key]
          ex[chosen[["exon_row_key"]], "exon_frame" := as.integer(chosen[["phase_int"]])]
        }
      }
    }
  }

  data.table::setorderv(ex, c("transcript_id", "exon_number", "exon_start", "exon_end"))
  frame_summary <- ex[, .(
    exon_frames = paste(as.integer(exon_frame), collapse = ",")
  ), by = transcript_id]
  if ("exon_frames" %in% names(tx)) tx[, "exon_frames" := NULL]
  tx <- merge(tx, frame_summary, by = "transcript_id", all.x = TRUE)
  tx[is.na(exon_frames), "exon_frames" := NA_character_]

  list(transcripts = tx, exons = ex)
}
