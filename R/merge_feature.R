# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.19
# Function: Merge Feature/GenePred-compatible annotation objects
# Input: Feature-compatible annotation objects
# Output: Merged Feature/GenePred-compatible object

prefix_ids_for_merge <- function(dt, prefix, id_cols = c("feature_id", "name", "gene_id", "transcript_id", "parent_id")) {
  dt <- data.table::copy(dt)
  for (col in intersect(id_cols, names(dt))) {
    idx <- !is.na(dt[[col]]) & nzchar(as.character(dt[[col]]))
    dt[idx, (col) := paste0(prefix, as.character(dt[[col]]))]
  }
  dt
}

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

