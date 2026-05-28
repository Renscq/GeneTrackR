# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.0
# Function: Summarize, slice, merge, and write GenePred objects
# Input: GenePred objects and genomic regions
# Output: Tables, subset objects, merged objects, and GenePred files

#' Summarize a GenePred object
#'
#' @param object A GenePred object.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start.
#' @param end Optional region end.
#' @param level Summary level.
#' @return A data.table summary.
#' @export
summary_genepred <- function(object, chrom = NULL, start = NULL, end = NULL, level = c("gene", "transcript", "exon")) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  level <- match.arg(level)
  check_region(chrom, start, end)

  obj <- if (!is.null(chrom)) {
    slice_genepred(
      object = object,
      chrom = chrom,
      start = start %||% 1L,
      end = end %||% .Machine$integer.max,
      mode = "overlap"
    )
  } else {
    object
  }

  if (level == "gene") {
    dt <- as_gene_table(obj)

    if (nrow(dt) == 0L) {
      return(data.table::data.table(
        chrom = character(),
        strand = character(),
        n_genes = integer(),
        n_coding_genes = integer(),
        n_noncoding_genes = integer(),
        median_gene_length = numeric()
      ))
    }

    return(dt[, .(
      n_genes = as.integer(.N),
      n_coding_genes = as.integer(sum(gene_type == "coding", na.rm = TRUE)),
      n_noncoding_genes = as.integer(sum(gene_type != "coding" | is.na(gene_type), na.rm = TRUE)),
      median_gene_length = as.numeric(stats::median(as.numeric(gene_end - gene_start + 1L), na.rm = TRUE))
    ), by = .(chrom, strand)][order(chrom, strand)])
  }

  if (level == "transcript") {
    dt <- as_transcript_table(obj)

    if (nrow(dt) == 0L) {
      return(data.table::data.table(
        chrom = character(),
        strand = character(),
        n_transcripts = integer(),
        n_coding_transcripts = integer(),
        n_noncoding_transcripts = integer(),
        median_transcript_length = numeric(),
        median_exon_count = numeric()
      ))
    }

    return(dt[, .(
      n_transcripts = as.integer(.N),
      n_coding_transcripts = as.integer(sum(gene_type == "coding", na.rm = TRUE)),
      n_noncoding_transcripts = as.integer(sum(gene_type != "coding" | is.na(gene_type), na.rm = TRUE)),
      median_transcript_length = as.numeric(stats::median(as.numeric(tx_end - tx_start + 1L), na.rm = TRUE)),
      median_exon_count = as.numeric(stats::median(as.numeric(exon_count), na.rm = TRUE))
    ), by = .(chrom, strand)][order(chrom, strand)])
  }

  dt <- as_exon_table(obj)

  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      chrom = character(),
      strand = character(),
      n_exons = integer(),
      median_exon_length = numeric()
    ))
  }

  dt[, .(
    n_exons = as.integer(.N),
    median_exon_length = as.numeric(stats::median(as.numeric(exon_end - exon_start + 1L), na.rm = TRUE))
  ), by = .(chrom, strand)][order(chrom, strand)]
}

#' Slice a GenePred object by genomic region
#'
#' @param object A GenePred object.
#' @param chrom Chromosome name.
#' @param start Region start in internal 1-based closed coordinates.
#' @param end Region end in internal 1-based closed coordinates.
#' @param mode Selection mode. `within` keeps only fully contained transcripts, `overlap` keeps any overlapping transcript without trimming, and `trim` clips transcript/exon/CDS coordinates to the requested interval.
#' @return A sliced GenePred object.
#' @examples
#' \dontrun{
#' slice_genepred(gp, chrom = "chr1", start = 1, end = 1000, mode = "overlap")
#' slice_genepred(gp, chrom = "chr1", start = 1, end = 1000, mode = "trim")
#' }
#' @export
slice_genepred <- function(object, chrom, start, end, mode = c("within", "overlap", "trim")) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  mode <- match.arg(mode)
  check_region(chrom, start, end)

  tx <- data.table::copy(data.table::as.data.table(object$transcripts))
  ex <- data.table::copy(data.table::as.data.table(object$exons))
  chrom_value <- as.character(chrom)
  start_value <- as.integer(start)
  end_value <- as.integer(end)

  if (mode == "within") {
    keep_idx <- tx[["chrom"]] == chrom_value &
      tx[["tx_start"]] >= start_value &
      tx[["tx_end"]] <= end_value
  } else {
    keep_idx <- tx[["chrom"]] == chrom_value &
      tx[["tx_start"]] <= end_value &
      tx[["tx_end"]] >= start_value
  }

  keep_tx <- tx[["transcript_id"]][keep_idx]
  tx <- tx[tx[["transcript_id"]] %in% keep_tx]
  ex <- ex[ex[["transcript_id"]] %in% keep_tx]

  if (mode == "trim" && nrow(tx) > 0L) {
    tx[, "tx_start" := as.integer(pmax(tx[["tx_start"]], start_value))]
    tx[, "tx_end" := as.integer(pmin(tx[["tx_end"]], end_value))]
    tx[, "cds_start" := as.integer(pmax(tx[["cds_start"]], start_value))]
    tx[, "cds_end" := as.integer(pmin(tx[["cds_end"]], end_value))]

    ex[, "exon_start" := as.integer(pmax(ex[["exon_start"]], start_value))]
    ex[, "exon_end" := as.integer(pmin(ex[["exon_end"]], end_value))]
    ex <- ex[ex[["exon_start"]] <= ex[["exon_end"]]]
    tx <- tx[tx[["transcript_id"]] %in% ex[["transcript_id"]]]

    exon_counts <- ex[, .(exon_count_new = as.integer(.N)), by = "transcript_id"]
    tx <- merge(tx[, setdiff(names(tx), "exon_count"), with = FALSE], exon_counts, by = "transcript_id", all.x = TRUE)
    data.table::setnames(tx, "exon_count_new", "exon_count")
  }

  GenePred(
    transcripts = tx,
    exons = ex,
    genes = build_gene_table(tx),
    meta = object$meta,
    validation = object$validation
  )
}

#' Merge GenePred objects
#'
#' @param ... GenePred objects.
#' @param conflict Conflict strategy for duplicated transcript IDs. Use `error` to stop, `keep_first` to keep the first annotation, `keep_all` to allow duplicates, or `rename` to add source prefixes.
#' @param rename_prefix Optional prefixes used when conflict is rename.
#' @param sort Sort merged records by genomic coordinate.
#' @return A merged GenePred object.
#' @examples
#' \dontrun{
#' merge_genepred(gp1, gp2, conflict = "error")
#' merge_genepred(gp1, gp2, conflict = "rename", rename_prefix = c("ref", "new"))
#' }
#' @export
merge_genepred <- function(..., conflict = c("error", "keep_first", "keep_all", "rename"), rename_prefix = NULL, sort = TRUE) {
  objs <- list(...)
  stop_if_not(length(objs) >= 1L, "At least one GenePred object is required.")
  stop_if_not(all(vapply(objs, inherits, logical(1), "GenePred")), "All inputs must be GenePred objects.")
  conflict <- match.arg(conflict)

  tx_list <- lapply(seq_along(objs), function(i) {
    x <- data.table::copy(objs[[i]]$transcripts)
    x[, source_index := i]
    x
  })
  ex_list <- lapply(seq_along(objs), function(i) {
    x <- data.table::copy(objs[[i]]$exons)
    x[, source_index := i]
    x
  })

  tx <- data.table::rbindlist(tx_list, fill = TRUE)
  ex <- data.table::rbindlist(ex_list, fill = TRUE)

  dup_ids <- tx[duplicated(transcript_id) | duplicated(transcript_id, fromLast = TRUE), unique(transcript_id)]
  if (length(dup_ids) > 0L) {
    if (conflict == "error") {
      stop("Duplicated transcript IDs found. Use `conflict = 'rename'`, 'keep_first', or 'keep_all'.", call. = FALSE)
    }
    if (conflict == "keep_first") {
      keep <- tx[, .I[1L], by = transcript_id]$V1
      tx <- tx[keep]
      ex <- ex[transcript_id %in% tx$transcript_id]
    }
    if (conflict == "rename") {
      if (is.null(rename_prefix)) {
        rename_prefix <- paste0("set", seq_along(objs), "_")
      }
      stop_if_not(length(rename_prefix) == length(objs), "`rename_prefix` must match the number of input objects.")
      tx[, old_transcript_id := transcript_id]
      tx[transcript_id %in% dup_ids, transcript_id := paste0(rename_prefix[source_index], transcript_id)]
      tx[transcript_id %in% dup_ids, gene_id := paste0(rename_prefix[source_index], gene_id)]
      map <- tx[, .(source_index, old_transcript_id, transcript_id, gene_id)]
      ex <- merge(ex, map, by.x = c("source_index", "transcript_id"), by.y = c("source_index", "old_transcript_id"), all.x = TRUE)
      ex[!is.na(i.transcript_id), transcript_id := i.transcript_id]
      ex[!is.na(i.gene_id), gene_id := i.gene_id]
      ex[, c("i.transcript_id", "i.gene_id") := NULL]
      tx[, old_transcript_id := NULL]
    }
  }

  tx[, source_index := NULL]
  ex[, source_index := NULL]
  if (sort) {
    data.table::setorder(tx, chrom, tx_start, tx_end, gene_id, transcript_id)
    data.table::setorder(ex, chrom, exon_start, exon_end, gene_id, transcript_id)
  }

  GenePred(
    transcripts = tx,
    exons = ex,
    genes = build_gene_table(tx),
    meta = list(format = "merged", coordinate_internal = "1-based closed"),
    validation = make_empty_validation()
  )
}

#' Write a GenePred object
#'
#' @param object A GenePred object.
#' @param file Output file path.
#' @param format Output format.
#' @param coordinate Output coordinate system.
#' @return Invisibly returns the output file path.
#' @export
write_genepred <- function(object, file, format = c("genePred", "genePredExt"), coordinate = c("ucsc", "granges")) {
  stop_if_not(inherits(object, "GenePred"), "`object` must be a GenePred object.")
  format <- match.arg(format)
  coordinate <- match.arg(coordinate)

  tx <- data.table::copy(object$transcripts)
  ex <- data.table::copy(object$exons)
  data.table::setorder(ex, transcript_id, exon_start, exon_end)

  exon_agg <- ex[, .(
    exonStarts = paste_comma_integer(if (coordinate == "ucsc") exon_start - 1L else exon_start),
    exonEnds = paste_comma_integer(exon_end),
    exonFrames = paste_comma_integer(exon_frame)
  ), by = transcript_id]

  out <- merge(tx, exon_agg, by = "transcript_id", all.x = TRUE)
  out[, `:=`(
    name = transcript_id,
    txStart = if (coordinate == "ucsc") tx_start - 1L else tx_start,
    txEnd = tx_end,
    cdsStart = if (coordinate == "ucsc") cds_start - 1L else cds_start,
    cdsEnd = cds_end,
    exonCount = exon_count,
    name2 = gene_id
  )]

  standard_cols <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds")
  if (format == "genePred") {
    data.table::fwrite(out[, ..standard_cols], file, sep = "\t", col.names = FALSE)
  } else {
    out[, `:=`(
      score = ifelse(is.na(score), 0, score),
      cdsStartStat = cds_start_stat %||% "unk",
      cdsEndStat = cds_end_stat %||% "unk"
    )]
    ext_cols <- c(standard_cols, "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
    data.table::fwrite(out[, ..ext_cols], file, sep = "\t", col.names = FALSE)
  }
  invisible(file)
}
