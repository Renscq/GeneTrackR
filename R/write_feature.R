# Author: Rensc
# Date: 2026-08-31
# Version: dev003
# Function: Write unified Feature annotation objects to GenePred, GFF, GTF, BED6, or BED12 files
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
#' @param format Output format. Use `auto`, `genepred`, `genepredext`, `gff`, `gtf`, `bed6`, or `bed12`.
#' `bed6` writes six-column BED intervals. For gene-model objects, BED6 is transcript-level. For generic interval objects, BED6 is feature-level.
#' `bed12` writes transcript-level BED12 gene models with exon blocks.
#' @param coordinate Output coordinate system for GenePred/BED-like formats.
#' `ucsc` means 0-based half-open where applicable; `granges` means 1-based closed.
#' @param overwrite Whether to overwrite an existing file.
#' @param keep_attributes Whether to reuse existing `attribute` strings for GFF/GTF output when available.
#' @param sort_output Whether to sort records before writing. Default TRUE. GenePred output is sorted by chromosome and transcript start. GFF/GTF output is sorted by chromosome, start, and feature hierarchy.
#' @param chrom_order Optional chromosome order used for sorting. It can be NULL, a character vector of chromosome names, a data frame whose first column contains chromosome names, or a FASTA index `.fai` file path. If NULL, a natural chromosome order is used.
#' @return Invisibly returns the output file path.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' outfile <- tempfile(fileext = ".bed")
#' write_feature(gp, outfile, format = "bed12", overwrite = TRUE)
#' @export
write_feature <- function(object,
                          file,
                          format = c("auto", "genepred", "genepredext", "gff", "gtf", "bed6", "bed12"),
                          coordinate = c("ucsc", "granges"),
                          overwrite = FALSE,
                          keep_attributes = TRUE,
                          sort_output = TRUE,
                          chrom_order = NULL) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"), "`object` must be a Feature/FeatureTrack or GenePred-compatible object.")
  format <- match.arg(format)
  coordinate <- match.arg(coordinate)
  if (format == "auto") format <- infer_write_feature_format(file)
  check_output_file(file, overwrite)

  if (format %in% c("genepred", "genepredext")) {
    gp <- as_genepred(object)
    write_feature_as_genepred(gp, file, format = format, coordinate = coordinate, sort_output = sort_output, chrom_order = chrom_order)
  } else if (format == "bed6") {
    write_feature_as_bed6(
      object,
      file,
      coordinate = coordinate,
      sort_output = sort_output,
      chrom_order = chrom_order
    )
  } else if (format == "bed12") {
    write_feature_as_bed12(
      object,
      file,
      coordinate = coordinate,
      sort_output = sort_output,
      chrom_order = chrom_order
    )
  } else if (format == "gff") {
    write_feature_as_gff(object, file, keep_attributes = keep_attributes, sort_output = sort_output, chrom_order = chrom_order)
  } else if (format == "gtf") {
    write_feature_as_gtf(object, file, keep_attributes = keep_attributes, sort_output = sort_output, chrom_order = chrom_order)
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
  if (grepl("\\.bed6$", x)) return("bed6")
  if (grepl("\\.bed12$", x)) return("bed12")
  if (grepl("\\.bed$", x)) return("bed12")
  stop("Cannot infer output format from file extension. Please set `format` explicitly.", call. = FALSE)
}

write_feature_as_genepred <- function(object, file, format = c("genepred", "genepredext"), coordinate = c("ucsc", "granges"), sort_output = TRUE, chrom_order = NULL) {
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

  if (isTRUE(sort_output)) {
    out <- sort_genepred_output_table(out, chrom_order = chrom_order)
  }

  standard_cols <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds")
  if (format == "genepred") {
    data.table::fwrite(out[, ..standard_cols], file, sep = "\t", col.names = FALSE, quote = FALSE)
  } else {
    ext_cols <- c(standard_cols, "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
    data.table::fwrite(out[, ..ext_cols], file, sep = "\t", col.names = FALSE, quote = FALSE)
  }
}

write_feature_as_bed6 <- function(object,
                                  file,
                                  coordinate = c("ucsc", "granges"),
                                  sort_output = TRUE,
                                  chrom_order = NULL) {
  out <- make_bed6_table(
    object,
    coordinate = coordinate,
    sort_output = sort_output,
    chrom_order = chrom_order
  )
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE, quote = FALSE)
}

write_feature_as_bed12 <- function(object,
                                   file,
                                   coordinate = c("ucsc", "granges"),
                                   sort_output = TRUE,
                                   chrom_order = NULL) {
  out <- make_bed12_table(
    object,
    coordinate = coordinate,
    sort_output = sort_output,
    chrom_order = chrom_order
  )
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE, quote = FALSE)
}

make_bed6_table <- function(object,
                            coordinate = c("ucsc", "granges"),
                            sort_output = TRUE,
                            chrom_order = NULL) {
  coordinate <- match.arg(coordinate)

  has_gene_model <- inherits(object, "GenePred") || (
    (inherits(object, "Feature") || inherits(object, "FeatureTrack")) &&
      !is.null(object$transcripts) && nrow(object$transcripts) > 0L
  )

  if (has_gene_model) {
    gp <- as_genepred(object)
    tx <- data.table::copy(data.table::as.data.table(gp$transcripts))
    if (!"score" %in% names(tx)) tx[, "score" := 0]
    out <- tx[, .(
      chrom = as.character(.SD[["chrom"]]),
      chromStart = if (coordinate == "ucsc") as.integer(.SD[["tx_start"]]) - 1L else as.integer(.SD[["tx_start"]]),
      chromEnd = as.integer(.SD[["tx_end"]]),
      name = as.character(.SD[["transcript_id"]]),
      score = normalize_bed_score(.SD[["score"]]),
      strand = normalize_bed_strand(.SD[["strand"]])
    ), .SDcols = intersect(names(tx), c("chrom", "tx_start", "tx_end", "transcript_id", "score", "strand"))]
  } else {
    dt <- as_feature_table(object)
    if (isTRUE(sort_output)) {
      dt <- sort_feature_output_table(dt, format = "bed", chrom_order = chrom_order)
    }
    name_vec <- pick_bed_feature_name(dt)
    out <- dt[, .(
      chrom = as.character(.SD[["chrom"]]),
      chromStart = if (coordinate == "ucsc") as.integer(.SD[["start"]]) - 1L else as.integer(.SD[["start"]]),
      chromEnd = as.integer(.SD[["end"]]),
      name = name_vec,
      score = normalize_bed_score(.SD[["score"]]),
      strand = normalize_bed_strand(.SD[["strand"]])
    ), .SDcols = intersect(names(dt), c("chrom", "start", "end", "score", "strand"))]
  }

  out <- out[!is.na(chrom) & !is.na(chromStart) & !is.na(chromEnd) & chromStart <= chromEnd]
  if (isTRUE(sort_output) && has_gene_model) {
    out <- sort_bed_output_table(out, chrom_order = chrom_order)
  }
  out[]
}

make_bed12_table <- function(object,
                             coordinate = c("ucsc", "granges"),
                             sort_output = TRUE,
                             chrom_order = NULL) {
  coordinate <- match.arg(coordinate)
  gp <- as_genepred(object)
  tx <- data.table::copy(data.table::as.data.table(gp$transcripts))
  ex <- data.table::copy(data.table::as.data.table(gp$exons))
  stop_if_not(nrow(tx) > 0L && nrow(ex) > 0L, "BED12 output requires non-empty transcript and exon tables.")

  if (!"score" %in% names(tx)) tx[, "score" := 0]
  if (!"cds_start" %in% names(tx)) tx[, "cds_start" := as.integer(tx_start)]
  if (!"cds_end" %in% names(tx)) tx[, "cds_end" := as.integer(tx_start - 1L)]
  if (!"gene_type" %in% names(tx)) {
    tx[, "gene_type" := data.table::fifelse(
      !is.na(as.integer(cds_start)) & !is.na(as.integer(cds_end)) & as.integer(cds_start) <= as.integer(cds_end),
      "coding",
      "non-coding"
    )]
  }

  ex[, `:=`(
    exon_start = as.integer(.SD[["exon_start"]]),
    exon_end = as.integer(.SD[["exon_end"]])
  ), .SDcols = c("exon_start", "exon_end")]
  data.table::setorder(ex, transcript_id, exon_start, exon_end)

  block_dt <- ex[!is.na(transcript_id) & !is.na(exon_start) & !is.na(exon_end) & exon_start <= exon_end,
    .(
      blockCount = .N,
      blockSizes = paste0(as.integer(exon_end - exon_start + 1L), collapse = ","),
      exonStartsInternal = paste0(as.integer(exon_start), collapse = ",")
    ),
    by = transcript_id
  ]

  out <- merge(tx, block_dt, by = "transcript_id", all.x = FALSE, sort = FALSE)
  tx_start0 <- if (coordinate == "ucsc") as.integer(out[["tx_start"]]) - 1L else as.integer(out[["tx_start"]])
  tx_start_internal <- as.integer(out[["tx_start"]])
  block_starts <- character(nrow(out))
  for (i in seq_len(nrow(out))) {
    starts_i <- as.integer(strsplit(as.character(out[["exonStartsInternal"]][i]), ",", fixed = TRUE)[[1L]])
    starts_i <- starts_i[!is.na(starts_i)]
    block_starts[i] <- paste0(starts_i - tx_start_internal[i], collapse = ",")
  }

  is_coding <- !is.na(out[["gene_type"]]) & as.character(out[["gene_type"]]) == "coding" &
    !is.na(as.integer(out[["cds_start"]])) & !is.na(as.integer(out[["cds_end"]])) &
    as.integer(out[["cds_start"]]) <= as.integer(out[["cds_end"]])

  thick_start <- tx_start0
  thick_end <- tx_start0
  thick_start[is_coding] <- if (coordinate == "ucsc") as.integer(out[["cds_start"]][is_coding]) - 1L else as.integer(out[["cds_start"]][is_coding])
  thick_end[is_coding] <- as.integer(out[["cds_end"]][is_coding])

  ans <- data.table::data.table(
    chrom = as.character(out[["chrom"]]),
    chromStart = tx_start0,
    chromEnd = as.integer(out[["tx_end"]]),
    name = as.character(out[["transcript_id"]]),
    score = normalize_bed_score(out[["score"]]),
    strand = normalize_bed_strand(out[["strand"]]),
    thickStart = as.integer(thick_start),
    thickEnd = as.integer(thick_end),
    itemRgb = "0",
    blockCount = as.integer(out[["blockCount"]]),
    blockSizes = as.character(out[["blockSizes"]]),
    blockStarts = as.character(block_starts)
  )

  ans <- ans[!is.na(chrom) & !is.na(chromStart) & !is.na(chromEnd) & chromStart <= chromEnd]
  if (isTRUE(sort_output)) {
    ans <- sort_bed_output_table(ans, chrom_order = chrom_order)
  }
  ans[]
}

normalize_bed_score <- function(x) {
  score <- suppressWarnings(as.numeric(x))
  score[is.na(score)] <- 0
  score <- pmin(pmax(round(score), 0), 1000)
  as.integer(score)
}

normalize_bed_strand <- function(x) {
  strand <- as.character(x)
  strand[is.na(strand) | !strand %in% c("+", "-")] <- "."
  strand
}

pick_bed_feature_name <- function(dt) {
  get_col <- function(nm) {
    if (nm %in% names(dt)) as.character(dt[[nm]]) else rep(NA_character_, nrow(dt))
  }
  name <- get_col("name")
  feature_id <- get_col("feature_id")
  transcript_id <- get_col("transcript_id")
  gene_id <- get_col("gene_id")
  type <- get_col("type")
  idx <- is.na(name) | name == "" | name %in% c("CDS", "UTR", "5UTR", "3UTR", "start_codon", "stop_codon", "exon", "gene", "transcript")
  name[idx & !is.na(feature_id) & feature_id != ""] <- feature_id[idx & !is.na(feature_id) & feature_id != ""]
  idx <- is.na(name) | name == ""
  name[idx & !is.na(transcript_id) & transcript_id != ""] <- transcript_id[idx & !is.na(transcript_id) & transcript_id != ""]
  idx <- is.na(name) | name == ""
  name[idx & !is.na(gene_id) & gene_id != ""] <- gene_id[idx & !is.na(gene_id) & gene_id != ""]
  idx <- is.na(name) | name == ""
  name[idx] <- paste0(type[idx], "_", seq_len(sum(idx)))
  name
}

sort_bed_output_table <- function(dt, chrom_order = NULL) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  dt[, "chrom_rank" := chrom_sort_rank(.SD[["chrom"]], chrom_order = chrom_order)]
  data.table::setorder(dt, chrom_rank, chrom, chromStart, chromEnd, name)
  dt[, "chrom_rank" := NULL]
  dt
}

write_feature_as_gff <- function(object, file, keep_attributes = TRUE, sort_output = TRUE, chrom_order = NULL) {
  dt <- as_feature_table(object)
  if (isTRUE(sort_output)) {
    dt <- sort_feature_output_table(dt, format = "gff", chrom_order = chrom_order)
  }
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
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE, quote = FALSE)
}

write_feature_as_gtf <- function(object, file, keep_attributes = TRUE, sort_output = TRUE, chrom_order = NULL) {
  dt <- as_feature_table(object)
  if (isTRUE(sort_output)) {
    dt <- sort_feature_output_table(dt, format = "gtf", chrom_order = chrom_order)
  }
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
  data.table::fwrite(out, file, sep = "\t", col.names = FALSE, quote = FALSE)
}


sort_genepred_output_table <- function(dt, chrom_order = NULL) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  dt[, "chrom_rank" := chrom_sort_rank(.SD[["chrom"]], chrom_order = chrom_order)]
  dt[, "strand_rank" := data.table::fifelse(.SD[["strand"]] == "+", 1L, data.table::fifelse(.SD[["strand"]] == "-", 2L, 3L))]
  data.table::setorder(dt, chrom_rank, chrom, txStart, txEnd, strand_rank, name)
  dt[, c("chrom_rank", "strand_rank") := NULL]
  dt
}

sort_feature_output_table <- function(dt, format = c("gff", "gtf", "bed"), chrom_order = NULL) {
  format <- match.arg(format)
  dt <- data.table::copy(data.table::as.data.table(dt))
  ensure_col <- function(x, nm, value) {
    if (!nm %in% names(x)) x[, (nm) := value]
    x
  }
  dt <- ensure_col(dt, "gene_id", NA_character_)
  dt <- ensure_col(dt, "transcript_id", NA_character_)
  dt <- ensure_col(dt, "feature_id", NA_character_)
  dt <- ensure_col(dt, "name", NA_character_)
  dt <- ensure_col(dt, "type", "feature")
  dt <- ensure_col(dt, "level", NA_character_)
  dt <- ensure_col(dt, "exon_number", NA_integer_)
  dt <- ensure_col(dt, "strand", "*")

  dt[, "chrom_rank" := chrom_sort_rank(.SD[["chrom"]], chrom_order = chrom_order)]
  dt[, "sort_group" := data.table::fifelse(
    !is.na(.SD[["gene_id"]]) & nzchar(.SD[["gene_id"]]),
    .SD[["gene_id"]],
    data.table::fifelse(
      !is.na(.SD[["transcript_id"]]) & nzchar(.SD[["transcript_id"]]),
      .SD[["transcript_id"]],
      data.table::fifelse(!is.na(.SD[["feature_id"]]) & nzchar(.SD[["feature_id"]]), .SD[["feature_id"]], .SD[["name"]])
    )
  )]

  group_ranges <- dt[, .(
    group_start = suppressWarnings(min(as.integer(.SD[["start"]]), na.rm = TRUE)),
    group_end = suppressWarnings(max(as.integer(.SD[["end"]]), na.rm = TRUE))
  ), by = sort_group]
  group_ranges[!is.finite(group_start), "group_start" := NA_integer_]
  group_ranges[!is.finite(group_end), "group_end" := NA_integer_]
  dt <- merge(dt, group_ranges, by = "sort_group", all.x = TRUE, sort = FALSE)
  dt[, "feature_rank" := feature_sort_rank(.SD[["type"]], .SD[["level"]])]
  dt[, "strand_rank" := data.table::fifelse(.SD[["strand"]] == "+", 1L, data.table::fifelse(.SD[["strand"]] == "-", 2L, 3L))]
  dt[, "exon_rank" := suppressWarnings(as.integer(.SD[["exon_number"]]))]
  dt[, "sort_start" := data.table::fifelse(
    tolower(as.character(.SD[["type"]])) %in% c("start_codon", "stop_codon"),
    as.integer(.SD[["group_end"]]) + as.integer(.SD[["feature_rank"]]),
    as.integer(.SD[["start"]])
  )]
  dt[is.na(exon_rank), "exon_rank" := 2147483647L]

  if (format == "bed") {
    data.table::setorder(dt, chrom_rank, chrom, start, end, name, type)
  } else {
    # GFF/GTF standard sorting: seqname -> start -> feature hierarchy -> end.
    # This keeps gene/transcript before exon/CDS/UTR when records share the same locus/start.
    data.table::setorder(
      dt,
      chrom_rank,
      chrom,
      sort_start,
      feature_rank,
      transcript_id,
      exon_rank,
      end,
      type
    )
  }
  drop_cols <- intersect(c("chrom_rank", "sort_group", "group_start", "group_end", "feature_rank", "strand_rank", "exon_rank", "sort_start"), names(dt))
  dt[, (drop_cols) := NULL]
  dt
}

feature_sort_rank <- function(type, level = NULL) {
  x <- tolower(as.character(type))
  lvl <- if (is.null(level)) rep(NA_character_, length(x)) else tolower(as.character(level))
  rank <- rep(50L, length(x))
  rank[x %in% c("gene") | lvl %in% c("gene")] <- 1L
  rank[x %in% c("mrna", "transcript") | lvl %in% c("transcript")] <- 2L
  rank[x %in% c("exon")] <- 3L
  rank[x %in% c("utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr", "5'utr", "3'utr")] <- 4L
  rank[x %in% c("cds")] <- 5L
  rank[x %in% c("start_codon")] <- 98L
  rank[x %in% c("stop_codon")] <- 99L
  rank
}

chrom_sort_rank <- function(chrom, chrom_order = NULL) {
  x <- as.character(chrom)
  order_vec <- normalize_chrom_order(chrom_order)
  if (!is.null(order_vec)) {
    ord <- stats::setNames(seq_along(order_vec), order_vec)
    rank <- unname(ord[x])
    missing <- is.na(rank)
    if (any(missing)) {
      rank[missing] <- 100000L + chrom_sort_rank(x[missing], chrom_order = NULL)
    }
    return(as.integer(rank))
  }

  clean <- toupper(x)
  clean <- sub("^CHR", "", clean)
  clean <- sub("^CHROMOSOME", "", clean)
  roman_map <- c(I = 1L, II = 2L, III = 3L, IV = 4L, V = 5L, VI = 6L, VII = 7L, VIII = 8L, IX = 9L, X = 10L,
                 XI = 11L, XII = 12L, XIII = 13L, XIV = 14L, XV = 15L, XVI = 16L, XVII = 17L, XVIII = 18L, XIX = 19L, XX = 20L)
  numeric_part <- suppressWarnings(as.integer(clean))
  roman_part <- unname(roman_map[clean])
  rank <- rep(100000L, length(clean))
  rank[!is.na(numeric_part)] <- numeric_part[!is.na(numeric_part)]
  rank[is.na(numeric_part) & !is.na(roman_part)] <- roman_part[is.na(numeric_part) & !is.na(roman_part)]
  rank[clean %in% c("M", "MT", "MITO", "MITOCHONDRIA")] <- 90000L
  rank[clean %in% c("C", "CP", "PT", "CHLOROPLAST", "PLASTID")] <- 90001L
  as.integer(rank)
}

normalize_chrom_order <- function(chrom_order = NULL) {
  if (is.null(chrom_order)) return(NULL)
  if (is.character(chrom_order) && length(chrom_order) == 1L && file.exists(chrom_order)) {
    fai <- data.table::fread(chrom_order, header = FALSE, sep = "\t", data.table = FALSE, showProgress = FALSE)
    return(as.character(fai[[1L]]))
  }
  if (is.data.frame(chrom_order)) {
    return(as.character(chrom_order[[1L]]))
  }
  if (is.character(chrom_order)) {
    return(as.character(chrom_order))
  }
  stop("`chrom_order` must be NULL, a character vector, a data frame, or a .fai file path.", call. = FALSE)
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
    has_attr <- !is.na(attr) & nzchar(attr) & attr != "."
  } else {
    attr <- rep(NA_character_, nrow(dt))
    has_attr <- rep(FALSE, nrow(dt))
  }

  out <- attr
  idx <- which(!has_attr)
  if (length(idx) > 0L) {
    get_col <- function(nm) {
      if (nm %in% names(dt)) as.character(dt[[nm]]) else rep(NA_character_, nrow(dt))
    }
    gene_id <- get_col("gene_id")
    transcript_id <- get_col("transcript_id")
    gene_name <- get_col("gene_name")
    exon_id <- get_col("exon_id")
    type <- tolower(as.character(dt$type))
    exon_number <- if ("exon_number" %in% names(dt)) suppressWarnings(as.integer(dt$exon_number)) else rep(NA_integer_, nrow(dt))

    gene_name[is.na(gene_name) | gene_name == ""] <- gene_id[is.na(gene_name) | gene_name == ""]
    exon_id[is.na(exon_id) | exon_id == ""] <- paste0(transcript_id[is.na(exon_id) | exon_id == ""], ".", exon_number[is.na(exon_id) | exon_id == ""])

    out[idx] <- vapply(idx, function(i) {
      vals <- list(gene_id = gene_id[i])
      if (!type[i] %in% c("gene")) {
        vals$transcript_id <- transcript_id[i]
      }
      if (type[i] %in% c("exon", "cds", "utr", "5utr", "3utr", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon")) {
        vals$exon_number <- if (!is.na(exon_number[i])) as.character(exon_number[i]) else NA_character_
        vals$exon_id <- exon_id[i]
      }
      vals$gene_name <- gene_name[i]
      vals <- vals[!is.na(unlist(vals, use.names = FALSE)) & nzchar(unlist(vals, use.names = FALSE))]
      if (length(vals) == 0L) return(".")
      paste0(paste(sprintf('%s "%s"', names(vals), unlist(vals, use.names = FALSE)), collapse = "; "), ";")
    }, character(1L))
  }
  out
}

normalize_gtf_feature_type <- function(x) {
  x <- as.character(x)
  lx <- tolower(x)
  x[lx == "mrna"] <- "transcript"
  x[lx %in% c("five_prime_utr", "five_utr", "5utr", "5'utr")] <- "5UTR"
  x[lx %in% c("three_prime_utr", "three_utr", "3utr", "3'utr")] <- "3UTR"
  x
}


#' Write a GenePred-compatible annotation object
#'
#' @description Backward-compatible wrapper around `write_feature()`.
#' @inheritParams write_feature
#' @return Invisibly returns the output file path.
#' @examples
#' gp_file <- system.file(
#'   "extdata", "gtr_demo.genePredExt", package = "GeneTrackR"
#' )
#' gp <- read_genepred(
#'   gp_file, format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' outfile <- tempfile(fileext = ".genePred")
#' write_genepred(gp, outfile)
#' @export
write_genepred <- function(object, file, format = c("genePred", "genePredExt"), coordinate = c("ucsc", "granges"), sort_output = TRUE, chrom_order = NULL) {
  format <- match.arg(format)
  fmt <- if (format == "genePred") "genepred" else "genepredext"
  write_feature(object = object, file = file, format = fmt, coordinate = coordinate, sort_output = sort_output, chrom_order = chrom_order)
}

#' Write a FeatureTrack object
#'
#' @description Backward-compatible wrapper around `write_feature()`.
#' @inheritParams write_feature
#' @return Invisibly returns the output file path.
#' @examples
#' bed_file <- system.file(
#'   "extdata", "gtr_demo_features.bed", package = "GeneTrackR"
#' )
#' features <- read_bed(bed_file, verbose = FALSE, progress = FALSE)
#' outfile <- tempfile(fileext = ".bed")
#' write_feature_track(features, outfile)
#' @export
write_feature_track <- function(object,
                                file,
                                format = c("bed12", "bed6", "gff", "gtf", "genepred", "genepredext"),
                                sort_output = TRUE,
                                chrom_order = NULL) {
  format <- match.arg(format)
  write_feature(
    object = object,
    file = file,
    format = format,
    sort_output = sort_output,
    chrom_order = chrom_order
  )
}
