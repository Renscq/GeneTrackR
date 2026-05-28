# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.11
# Function: Unified retrieve and merge APIs for feature, signal, and variant tracks
# Input: Feature, BwgTrack, and VariantTrack objects
# Output: Retrieved tables/objects and merged track objects

filter_by_region_internal <- function(dt, chrom = NULL, start = NULL, end = NULL,
                                      start_col = "start", end_col = "end",
                                      mode = c("overlap", "within")) {
  mode <- match.arg(mode)
  dt <- data.table::as.data.table(dt)
  if (!is.null(chrom) && "chrom" %in% names(dt)) {
    dt <- dt[dt[["chrom"]] %in% as.character(chrom)]
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

#' Retrieve annotation features from a Feature-compatible object
#'
#' @description
#' `retrieve_feature()` is the unified retrieval API for annotation objects.
#' It is the single retrieval API for annotation features. It can retrieve standardized feature rows,
#' gene rows, transcript rows, or exon rows from GenePred/GTF/GFF/BED-derived
#' Feature objects.
#'
#' @param object A Feature-compatible annotation object.
#' @param pattern Optional character pattern to search IDs/names.
#' @param level Retrieval level. Use `feature`, `gene`, `transcript`, or `exon`.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start.
#' @param end Optional region end.
#' @param mode Region selection mode. `overlap` keeps overlapping records and
#' `within` keeps records fully contained in the region.
#' @param gene_id Optional gene ID filter.
#' @param transcript_id Optional transcript ID filter.
#' @param type Optional feature type filter. Only used for `level = "feature"`.
#' @param ignore_case Logical. Whether pattern matching ignores case.
#' @param fixed Logical. Whether `pattern` should be treated as a fixed string.
#'
#' @return A data.table containing retrieved records.
#'
#' @examples
#' \dontrun{
#' retrieve_feature(gp, pattern = "GeneA")
#' retrieve_feature(gp, level = "gene", pattern = "GeneA")
#' retrieve_feature(gp, level = "transcript", gene_id = "GeneA")
#' retrieve_feature(gp, chrom = "chr1", start = 1, end = 1000)
#' }
#' @export
retrieve_feature <- function(object,
                             pattern = NULL,
                             level = c("feature", "gene", "transcript", "exon"),
                             chrom = NULL,
                             start = NULL,
                             end = NULL,
                             mode = c("overlap", "within", "trim"),
                             gene_id = NULL,
                             transcript_id = NULL,
                             type = NULL,
                             ignore_case = TRUE,
                             fixed = FALSE,
                             as = c("data.table", "Feature")) {
  stop_if_not(inherits(object, "Feature") || inherits(object, "FeatureTrack") || inherits(object, "GenePred"),
              "`object` must be a Feature-compatible annotation object.")
  level <- match.arg(level)
  mode <- match.arg(mode)
  as <- match.arg(as)

  if (as == "Feature") {
    stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end),
                "`chrom`, `start`, and `end` are required when `as = 'Feature'`.")
    stop_if_not(is_gene_model_feature(object) || level == "feature",
                "`as = 'Feature'` requires a Feature-compatible annotation object.")
    if (is_gene_model_feature(object) && level %in% c("feature", "gene", "transcript", "exon")) {
      gp <- as_genepred(object)
      tx <- data.table::copy(gp$transcripts)
      ex <- data.table::copy(gp$exons)
      chrom_value <- as.character(chrom)[1L]
      start_value <- as.integer(start)[1L]
      end_value <- as.integer(end)[1L]
      if (mode == "within") {
        keep_idx <- tx[["chrom"]] == chrom_value & tx[["tx_start"]] >= start_value & tx[["tx_end"]] <= end_value
      } else {
        keep_idx <- tx[["chrom"]] == chrom_value & tx[["tx_start"]] <= end_value & tx[["tx_end"]] >= start_value
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
      out <- Feature(
        data = genepred_to_feature_table(tx, ex, build_gene_table(tx)),
        genes = build_gene_table(tx),
        transcripts = tx,
        exons = ex,
        meta = modifyList(gp$meta, list(format = "retrieved", coordinate_internal = "1-based closed")),
        validation = make_empty_validation()
      )
      class(out) <- unique(c(class(out), "GenePred"))
      return(out)
    }
    dt_feature <- retrieve_feature(object, pattern = pattern, level = "feature", chrom = chrom, start = start, end = end,
                                   mode = if (mode == "trim") "overlap" else mode, gene_id = gene_id,
                                   transcript_id = transcript_id, type = type, ignore_case = ignore_case,
                                   fixed = fixed, as = "data.table")
    return(Feature(dt_feature, meta = modifyList(object$meta %||% list(), list(format = "retrieved", coordinate_internal = "1-based closed"))))
  }

  if (level == "feature") {
    dt <- as_feature_table(object)
    region_mode <- if (mode == "trim") "overlap" else mode
    dt <- filter_by_region_internal(dt, chrom = chrom, start = start, end = end, mode = region_mode)
    if (!is.null(gene_id) && "gene_id" %in% names(dt)) {
      dt <- dt[dt[["gene_id"]] %in% as.character(gene_id)]
    }
    if (!is.null(transcript_id) && "transcript_id" %in% names(dt)) {
      dt <- dt[dt[["transcript_id"]] %in% as.character(transcript_id)]
    }
    if (!is.null(type) && "type" %in% names(dt)) {
      dt <- dt[dt[["type"]] %in% as.character(type)]
    }
    dt <- match_pattern_internal(
      dt,
      pattern = pattern,
      fields = c("feature_id", "name", "gene_id", "transcript_id", "parent_id", "type"),
      ignore_case = ignore_case,
      fixed = fixed
    )
    if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "start", "end", "feature_id"))
    return(dt[])
  }

  gp <- as_genepred(object)
  if (level == "gene") {
    dt <- data.table::copy(gp$genes)
    region_mode <- if (mode == "trim") "overlap" else mode
    dt <- filter_by_region_internal(dt, chrom = chrom, start = start, end = end,
                                    start_col = "gene_start", end_col = "gene_end", mode = region_mode)
    if (!is.null(gene_id)) dt <- dt[dt[["gene_id"]] %in% as.character(gene_id)]
    dt <- match_pattern_internal(dt, pattern = pattern, fields = c("gene_id", "gene_type"),
                                 ignore_case = ignore_case, fixed = fixed)
    if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "gene_start", "gene_end", "gene_id"))
    return(dt[])
  }

  if (level == "transcript") {
    dt <- data.table::copy(gp$transcripts)
    region_mode <- if (mode == "trim") "overlap" else mode
    dt <- filter_by_region_internal(dt, chrom = chrom, start = start, end = end,
                                    start_col = "tx_start", end_col = "tx_end", mode = region_mode)
    if (!is.null(gene_id)) dt <- dt[dt[["gene_id"]] %in% as.character(gene_id)]
    if (!is.null(transcript_id)) dt <- dt[dt[["transcript_id"]] %in% as.character(transcript_id)]
    dt <- match_pattern_internal(dt, pattern = pattern, fields = c("transcript_id", "gene_id", "gene_type"),
                                 ignore_case = ignore_case, fixed = fixed)
    if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "tx_start", "tx_end", "transcript_id"))
    return(dt[])
  }

  dt <- data.table::copy(gp$exons)
  region_mode <- if (mode == "trim") "overlap" else mode
  dt <- filter_by_region_internal(dt, chrom = chrom, start = start, end = end,
                                  start_col = "exon_start", end_col = "exon_end", mode = region_mode)
  if (!is.null(gene_id)) dt <- dt[dt[["gene_id"]] %in% as.character(gene_id)]
  if (!is.null(transcript_id)) dt <- dt[dt[["transcript_id"]] %in% as.character(transcript_id)]
  dt <- match_pattern_internal(dt, pattern = pattern, fields = c("transcript_id", "gene_id"),
                               ignore_case = ignore_case, fixed = fixed)
  if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "exon_start", "exon_end", "transcript_id", "exon_number"))
  dt[]
}

#' Retrieve signal records from a BwgTrack object
#'
#' @description
#' A user-facing wrapper around `query_bwg()` that matches the unified
#' `retrieve_*` naming scheme.
#'
#' @inheritParams query_bwg
#' @param as Output type. Use `data.table`, `BwgTrack`, or `GRanges`.
#' @return A data.table, BwgTrack object, or GRanges object.
#' @export
retrieve_bwg <- function(object,
                         chrom,
                         start,
                         end,
                         samples = NULL,
                         strand = c("ignore", "+", "-", "both", "auto"),
                         strand_policy = c("ignore_unstranded", "strict"),
                         as = c("data.table", "BwgTrack", "GRanges"),
                         verbose = FALSE,
                         progress = interactive() && isTRUE(verbose),
                         keep_empty_samples = FALSE,
                         tabix_empty_fallback = NULL) {
  as <- match.arg(as)
  strand <- match.arg(strand)
  strand_policy <- match.arg(strand_policy)
  dt <- query_bwg(
    object = object,
    chrom = chrom,
    start = start,
    end = end,
    samples = samples,
    strand = strand,
    strand_policy = strand_policy,
    verbose = verbose,
    progress = progress,
    keep_empty_samples = keep_empty_samples,
    tabix_empty_fallback = tabix_empty_fallback
  )
  if (as == "data.table") return(dt)
  if (as == "GRanges") {
    return(GenomicRanges::GRanges(
      seqnames = dt$chrom,
      ranges = IRanges::IRanges(dt$start, dt$end),
      sample_id = dt$sample_id,
      score = dt$value,
      strand = dt$strand
    ))
  }
  sample_tbl <- object$samples
  if (!is.null(samples)) {
    sample_tbl <- sample_tbl[sample_tbl[["sample_id"]] %in% as.character(samples)]
  }
  BwgTrack(samples = sample_tbl, data = dt, meta = modifyList(object$meta, list(mode = "memory")), validation = make_empty_validation())
}

#' Retrieve variants from a VariantTrack object
#'
#' @param object A VariantTrack object.
#' @param pattern Optional pattern to search variant ID, alleles, INFO, or type.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start.
#' @param end Optional region end.
#' @param variant_id Optional variant ID filter.
#' @param variant_type Optional variant type filter.
#' @param ignore_case Logical. Whether pattern matching ignores case.
#' @param fixed Logical. Whether `pattern` should be treated as fixed text.
#' @param as Output type. Use `data.table` or `VariantTrack`.
#'
#' @return A data.table or VariantTrack object.
#' @export
retrieve_vcf <- function(object,
                         pattern = NULL,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         variant_id = NULL,
                         variant_type = NULL,
                         ignore_case = TRUE,
                         fixed = FALSE,
                         as = c("data.table", "VariantTrack")) {
  stop_if_not(inherits(object, "VariantTrack"), "`object` must be a VariantTrack object.")
  as <- match.arg(as)
  dt <- data.table::copy(object$data)
  if (!is.null(chrom)) dt <- dt[dt[["chrom"]] %in% as.character(chrom)]
  if (!is.null(start) && !is.null(end)) {
    s <- as.integer(start)[1L]
    e <- as.integer(end)[1L]
    dt <- dt[dt[["pos"]] >= s & dt[["pos"]] <= e]
  }
  if (!is.null(variant_id)) dt <- dt[dt[["variant_id"]] %in% as.character(variant_id)]
  if (!is.null(variant_type)) dt <- dt[dt[["variant_type"]] %in% as.character(variant_type)]
  dt <- match_pattern_internal(dt, pattern = pattern,
                               fields = c("variant_id", "ref", "alt", "info", "variant_type"),
                               ignore_case = ignore_case, fixed = fixed)
  if (nrow(dt) > 0L) data.table::setorderv(dt, c("chrom", "pos", "variant_id"))
  if (as == "data.table") return(dt[])
  VariantTrack(dt, meta = object$meta)
}

prefix_ids_for_merge <- function(dt, prefix, id_cols = c("feature_id", "name", "gene_id", "transcript_id", "parent_id")) {
  dt <- data.table::copy(dt)
  for (col in intersect(id_cols, names(dt))) {
    idx <- !is.na(dt[[col]]) & nzchar(as.character(dt[[col]]))
    dt[idx, (col) := paste0(prefix, as.character(dt[[col]]))]
  }
  dt
}

#' Merge Feature-compatible annotation objects
#'
#' @description
#' `merge_feature()` is the unified merge API for GenePred/GTF/GFF/BED-derived
#' annotation objects.
#'
#' @param ... Feature-compatible annotation objects.
#' @param source_names Optional source labels stored in `track_source`.
#' @param conflict Conflict strategy for duplicated feature/gene/transcript IDs:
#' `keep_all`, `error`, `rename`, or `keep_first`.
#' @param rename_prefix Optional prefixes used when `conflict = "rename"`.
#' @param sort Logical. Whether to sort merged tables.
#'
#' @return A merged Feature object. If gene model tables are present, the object
#' also inherits from GenePred for compatibility with gene model plotting.
#' @export
merge_feature <- function(...,
                          source_names = NULL,
                          conflict = c("keep_all", "error", "rename", "keep_first"),
                          rename_prefix = NULL,
                          sort = TRUE) {
  objs <- list(...)
  stop_if_not(length(objs) > 0L, "At least one Feature-compatible object is required.")
  stop_if_not(all(vapply(objs, function(x) inherits(x, "Feature") || inherits(x, "FeatureTrack") || inherits(x, "GenePred"), logical(1L))),
              "All inputs must be Feature-compatible objects.")
  conflict <- match.arg(conflict)
  if (is.null(source_names)) source_names <- paste0("track", seq_along(objs))
  stop_if_not(length(source_names) == length(objs), "`source_names` must match the number of input objects.")

  if (conflict == "rename" && is.null(rename_prefix)) {
    rename_prefix <- paste0(source_names, "_")
  }
  if (!is.null(rename_prefix)) {
    stop_if_not(length(rename_prefix) == length(objs), "`rename_prefix` must match the number of input objects.")
  }

  data_list <- vector("list", length(objs))
  gene_list <- list()
  tx_list <- list()
  exon_list <- list()

  for (i in seq_along(objs)) {
    x <- objs[[i]]
    prefix <- if (!is.null(rename_prefix)) as.character(rename_prefix[i]) else ""
    dt <- as_feature_table(x)
    if (conflict == "rename") dt <- prefix_ids_for_merge(dt, prefix)
    dt[, "track_source" := as.character(source_names[i])]
    data_list[[i]] <- dt

    if (!is.null(x$genes) && nrow(x$genes) > 0L) {
      g <- data.table::copy(x$genes)
      if (conflict == "rename") {
        idx <- !is.na(g[["gene_id"]]) & nzchar(g[["gene_id"]])
        g[idx, "gene_id" := paste0(prefix, g[["gene_id"]])]
      }
      g[, "track_source" := as.character(source_names[i])]
      gene_list[[length(gene_list) + 1L]] <- g
    }
    if (!is.null(x$transcripts) && nrow(x$transcripts) > 0L) {
      tx <- data.table::copy(x$transcripts)
      if (conflict == "rename") {
        for (col in intersect(c("gene_id", "transcript_id"), names(tx))) {
          idx <- !is.na(tx[[col]]) & nzchar(tx[[col]])
          tx[idx, (col) := paste0(prefix, tx[[col]])]
        }
      }
      tx[, "track_source" := as.character(source_names[i])]
      tx_list[[length(tx_list) + 1L]] <- tx
    }
    if (!is.null(x$exons) && nrow(x$exons) > 0L) {
      ex <- data.table::copy(x$exons)
      if (conflict == "rename") {
        for (col in intersect(c("gene_id", "transcript_id"), names(ex))) {
          idx <- !is.na(ex[[col]]) & nzchar(ex[[col]])
          ex[idx, (col) := paste0(prefix, ex[[col]])]
        }
      }
      ex[, "track_source" := as.character(source_names[i])]
      exon_list[[length(exon_list) + 1L]] <- ex
    }
  }

  data <- data.table::rbindlist(data_list, fill = TRUE)
  genes <- if (length(gene_list) > 0L) data.table::rbindlist(gene_list, fill = TRUE) else NULL
  transcripts <- if (length(tx_list) > 0L) data.table::rbindlist(tx_list, fill = TRUE) else NULL
  exons <- if (length(exon_list) > 0L) data.table::rbindlist(exon_list, fill = TRUE) else NULL

  if (conflict == "error") {
    if ("feature_id" %in% names(data)) {
      dup_feature <- data[duplicated(feature_id) | duplicated(feature_id, fromLast = TRUE), unique(feature_id)]
      stop_if_not(length(dup_feature) == 0L, "Duplicated feature IDs found. Use `conflict = 'rename'`, `keep_all`, or `keep_first`.")
    }
    if (!is.null(transcripts) && "transcript_id" %in% names(transcripts)) {
      dup_tx <- transcripts[duplicated(transcript_id) | duplicated(transcript_id, fromLast = TRUE), unique(transcript_id)]
      stop_if_not(length(dup_tx) == 0L, "Duplicated transcript IDs found. Use `conflict = 'rename'`, `keep_all`, or `keep_first`.")
    }
  }

  if (conflict == "keep_first") {
    if ("feature_id" %in% names(data)) data <- data[!duplicated(feature_id)]
    if (!is.null(genes) && "gene_id" %in% names(genes)) genes <- genes[!duplicated(gene_id)]
    if (!is.null(transcripts) && "transcript_id" %in% names(transcripts)) transcripts <- transcripts[!duplicated(transcript_id)]
    if (!is.null(exons) && "transcript_id" %in% names(exons)) {
      exons <- exons[transcript_id %in% transcripts$transcript_id]
    }
  }

  if (isTRUE(sort)) {
    if (nrow(data) > 0L) data.table::setorderv(data, intersect(c("chrom", "start", "end", "gene_id", "transcript_id", "feature_id"), names(data)))
    if (!is.null(genes) && nrow(genes) > 0L) data.table::setorderv(genes, intersect(c("chrom", "gene_start", "gene_end", "gene_id"), names(genes)))
    if (!is.null(transcripts) && nrow(transcripts) > 0L) data.table::setorderv(transcripts, intersect(c("chrom", "tx_start", "tx_end", "gene_id", "transcript_id"), names(transcripts)))
    if (!is.null(exons) && nrow(exons) > 0L) data.table::setorderv(exons, intersect(c("chrom", "exon_start", "exon_end", "gene_id", "transcript_id", "exon_number"), names(exons)))
  }

  out <- Feature(
    data = data,
    genes = genes,
    transcripts = transcripts,
    exons = exons,
    meta = list(format = "merged", coordinate_internal = "1-based closed"),
    validation = make_empty_validation()
  )
  if (!is.null(transcripts) && nrow(transcripts) > 0L && !is.null(exons) && nrow(exons) > 0L) {
    class(out) <- unique(c(class(out), "GenePred"))
  }
  out
}

#' Merge VCF-like VariantTrack objects
#'
#' @param ... VariantTrack objects.
#' @param source_names Optional source labels stored in `track_source`.
#' @param conflict Conflict strategy for duplicated variant IDs.
#' @param rename_prefix Optional prefixes when `conflict = "rename"`.
#' @param sort Logical. Whether to sort the merged variants.
#'
#' @return A merged VariantTrack object.
#' @export
merge_vcf <- function(...,
                      source_names = NULL,
                      conflict = c("keep_all", "error", "rename", "keep_first"),
                      rename_prefix = NULL,
                      sort = TRUE) {
  tracks <- list(...)
  stop_if_not(length(tracks) > 0L, "At least one VariantTrack object is required.")
  stop_if_not(all(vapply(tracks, inherits, logical(1L), "VariantTrack")), "All inputs must be VariantTrack objects.")
  conflict <- match.arg(conflict)
  if (is.null(source_names)) source_names <- paste0("track", seq_along(tracks))
  stop_if_not(length(source_names) == length(tracks), "`source_names` must match the number of input objects.")
  if (conflict == "rename" && is.null(rename_prefix)) rename_prefix <- paste0(source_names, "_")
  if (!is.null(rename_prefix)) stop_if_not(length(rename_prefix) == length(tracks), "`rename_prefix` must match the number of input objects.")

  out <- lapply(seq_along(tracks), function(i) {
    dt <- data.table::copy(tracks[[i]]$data)
    if (conflict == "rename" && "variant_id" %in% names(dt)) {
      idx <- !is.na(dt[["variant_id"]]) & nzchar(dt[["variant_id"]])
      dt[idx, "variant_id" := paste0(rename_prefix[i], dt[["variant_id"]])]
    }
    dt[, "track_source" := as.character(source_names[i])]
    dt
  })
  dt <- data.table::rbindlist(out, fill = TRUE)

  if (conflict == "error" && "variant_id" %in% names(dt)) {
    dup <- dt[duplicated(variant_id) | duplicated(variant_id, fromLast = TRUE), unique(variant_id)]
    stop_if_not(length(dup) == 0L, "Duplicated variant IDs found. Use `conflict = 'rename'`, `keep_all`, or `keep_first`.")
  }
  if (conflict == "keep_first" && "variant_id" %in% names(dt)) {
    dt <- dt[!duplicated(variant_id)]
  }
  if (isTRUE(sort) && nrow(dt) > 0L) {
    data.table::setorderv(dt, intersect(c("chrom", "pos", "variant_id"), names(dt)))
  }
  VariantTrack(dt, meta = list(format = "merged", coordinate_internal = "1-based position"))
}

# Redundant search/slice/merge wrapper APIs were removed in 0.2.14.
