# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.28
# Function: Internal plotting utilities for gene models and signal tracks
# Input: Annotation and signal tables
# Output: ggplot objects and transformed plotting tables

make_feature_segments <- function(tx, exons) {
  if (nrow(tx) == 0L || nrow(exons) == 0L) {
    return(data.table::data.table())
  }

  tx_small <- data.table::as.data.table(tx)[
    , c("transcript_id", "cds_start", "cds_end", "gene_type"), with = FALSE
  ]
  ex <- merge(
    data.table::copy(data.table::as.data.table(exons)),
    tx_small,
    by = "transcript_id",
    all.x = TRUE,
    sort = FALSE
  )

  ex[, `:=`(
    exon_start = as.integer(exon_start),
    exon_end = as.integer(exon_end),
    cds_start = as.integer(cds_start),
    cds_end = as.integer(cds_end),
    gene_type = as.character(gene_type)
  )]

  base_cols <- c("transcript_id", "gene_id", "chrom", "strand", "exon_number")
  make_segment <- function(dt, feature, start, end) {
    if (nrow(dt) == 0L) {
      return(data.table::data.table())
    }
    out <- dt[, base_cols, with = FALSE]
    out[, `:=`(
      feature = as.character(feature),
      start = as.integer(start),
      end = as.integer(end)
    )]
    out[!is.na(start) & !is.na(end) & start <= end]
  }

  noncoding <- ex[
    is.na(gene_type) | gene_type != "coding" |
      is.na(cds_start) | is.na(cds_end) | cds_start > cds_end
  ]
  coding <- ex[
    !(is.na(gene_type) | gene_type != "coding" |
        is.na(cds_start) | is.na(cds_end) | cds_start > cds_end)
  ]

  out <- list(
    make_segment(noncoding, "exon", noncoding$exon_start, noncoding$exon_end)
  )

  if (nrow(coding) > 0L) {
    out[[length(out) + 1L]] <- make_segment(
      coding,
      "UTR",
      coding$exon_start,
      pmin(coding$exon_end, coding$cds_start - 1L)
    )
    out[[length(out) + 1L]] <- make_segment(
      coding,
      "CDS",
      pmax(coding$exon_start, coding$cds_start),
      pmin(coding$exon_end, coding$cds_end)
    )
    out[[length(out) + 1L]] <- make_segment(
      coding,
      "UTR",
      pmax(coding$exon_start, coding$cds_end + 1L),
      coding$exon_end
    )
  }

  out <- data.table::rbindlist(out, fill = TRUE)
  if (nrow(out) == 0L) {
    return(data.table::data.table())
  }
  data.table::setorderv(out, c("transcript_id", "exon_number", "start", "end", "feature"))
  out[]
}


map_to_transcript_coordinate <- function(exons, segments = NULL) {
  ex <- data.table::copy(data.table::as.data.table(exons))
  data.table::setorderv(ex, c("transcript_id", "exon_start", "exon_end"))

  ex[, "exon_width" := as.integer(ex[["exon_end"]] - ex[["exon_start"]] + 1L)]
  ex[, "exon_offset" := as.integer(cumsum(data.table::shift(ex[["exon_width"]], fill = 0L))), by = "transcript_id"]
  ex[, "plot_start" := as.integer(ex[["exon_offset"]] + 1L)]
  ex[, "plot_end" := as.integer(ex[["exon_offset"]] + ex[["exon_width"]])]

  if (is.null(segments)) {
    return(list(exons = ex, segments = NULL))
  }

  seg <- merge(
    data.table::copy(data.table::as.data.table(segments)),
    ex[, c("transcript_id", "exon_number", "exon_start", "exon_offset"), with = FALSE],
    by = c("transcript_id", "exon_number"),
    all.x = TRUE
  )
  seg[, "plot_start" := as.integer(seg[["exon_offset"]] + (seg[["start"]] - seg[["exon_start"]]) + 1L)]
  seg[, "plot_end" := as.integer(seg[["exon_offset"]] + (seg[["end"]] - seg[["exon_start"]]) + 1L)]

  list(exons = ex, segments = seg)
}

add_highlight_layer <- function(p, highlight, ymin = -Inf, ymax = Inf) {
  if (is.null(highlight) || nrow(highlight) == 0L) {
    return(p)
  }
  h <- data.table::as.data.table(highlight)
  if (!all(c("start", "end") %in% names(h))) {
    stop("`highlight` must contain `start` and `end` columns.", call. = FALSE)
  }
  p + ggplot2::geom_rect(
    data = h,
    ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    alpha = 0.15
  )
}

collapse_transcripts <- function(tx, exons, collapse = c("none", "union_exon", "longest")) {
  collapse <- match.arg(collapse)
  tx <- data.table::as.data.table(tx)
  exons <- data.table::as.data.table(exons)

  if (collapse == "none" || nrow(tx) == 0L) {
    return(list(transcripts = tx, exons = exons))
  }

  if (collapse == "longest") {
    tx_copy <- data.table::copy(tx)
    tx_copy[, "tx_width" := as.numeric(tx_copy[["tx_end"]] - tx_copy[["tx_start"]] + 1L)]
    data.table::setorderv(tx_copy, c("gene_id", "tx_width"), order = c(1L, -1L))
    keep <- tx_copy[, .SD[1L], by = "gene_id"]$transcript_id
    return(list(
      transcripts = tx[tx[["transcript_id"]] %in% keep],
      exons = exons[exons[["transcript_id"]] %in% keep]
    ))
  }

  genes <- build_gene_table(tx)

  union_ex <- exons[
    , .(
      exon_start = as.integer(min(.SD[["exon_start"]], na.rm = TRUE)),
      exon_end = as.integer(max(.SD[["exon_end"]], na.rm = TRUE))
    ),
    by = c("gene_id", "chrom", "strand", "exon_number"),
    .SDcols = c("exon_start", "exon_end")
  ]
  union_ex[, "transcript_id" := paste0(union_ex[["gene_id"]], "_collapsed")]
  union_ex[, "exon_frame" := NA_integer_]
  union_ex <- union_ex[, c("transcript_id", "gene_id", "chrom", "strand", "exon_number", "exon_start", "exon_end", "exon_frame"), with = FALSE]

  exon_counts <- union_ex[, .(exon_count = as.integer(.N)), by = "gene_id"]
  tx_new <- merge(genes, exon_counts, by = "gene_id", all.x = TRUE)
  tx_new[, "transcript_id" := paste0(tx_new[["gene_id"]], "_collapsed")]
  tx_new[, "tx_start" := as.integer(tx_new[["gene_start"]])]
  tx_new[, "tx_end" := as.integer(tx_new[["gene_end"]])]
  tx_new[, "cds_start" := as.integer(tx_new[["gene_start"]])]
  tx_new[, "cds_end" := as.integer(tx_new[["gene_start"]] - 1L)]
  tx_new <- tx_new[, c(
    "transcript_id", "gene_id", "chrom", "strand",
    "tx_start", "tx_end", "cds_start", "cds_end",
    "exon_count", "gene_type"
  ), with = FALSE]

  list(transcripts = tx_new, exons = union_ex)
}

make_extended_brewer_palette <- function(n, color_palette = "Paired") {
  n <- max(1L, as.integer(n))
  color_palette <- as.character(color_palette)[1L]
  if (is.na(color_palette) || !nzchar(color_palette)) {
    color_palette <- "Paired"
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_info <- RColorBrewer::brewer.pal.info[color_palette, , drop = FALSE]
    max_colors <- as.integer(pal_info[["maxcolors"]][1L])
    min_colors <- max(3L, min(max_colors, n))
    base_n <- max(3L, min(max_colors, max(3L, min_colors)))
    base <- RColorBrewer::brewer.pal(base_n, color_palette)
    return(grDevices::colorRampPalette(base)(n))
  }

  warning(
    sprintf("Unknown `color_palette`: %s. Falling back to 'Paired'.", color_palette),
    call. = FALSE
  )
  base <- RColorBrewer::brewer.pal(12L, "Paired")
  grDevices::colorRampPalette(base)(n)
}

normalize_discrete_fill_colors <- function(n, color_palette = "Paired", fill_colors = NULL) {
  n <- max(1L, as.integer(n))

  if (is.null(fill_colors)) {
    return(make_extended_brewer_palette(n, color_palette = color_palette))
  }

  fill_colors <- as.character(fill_colors)
  fill_colors <- fill_colors[!is.na(fill_colors) & nzchar(fill_colors)]
  if (length(fill_colors) == 0L) {
    return(make_extended_brewer_palette(n, color_palette = color_palette))
  }

  if (length(fill_colors) >= n) {
    return(fill_colors[seq_len(n)])
  }

  if (length(fill_colors) == 1L) {
    return(rep(fill_colors, n))
  }

  grDevices::colorRampPalette(fill_colors)(n)
}

apply_discrete_fill_scale <- function(p, color_palette = "Paired", fill_colors = NULL) {
  if (!is.null(fill_colors)) {
    fill_colors <- as.character(fill_colors)
    color_names <- names(fill_colors)
    if (!is.null(color_names) && any(nzchar(color_names))) {
      return(p + ggplot2::scale_fill_manual(values = fill_colors, na.value = "grey80"))
    }
  }

  p + ggplot2::discrete_scale(
    aesthetics = "fill",
    scale_name = "extended_brewer_fill",
    palette = function(n) normalize_discrete_fill_colors(
      n = n,
      color_palette = color_palette,
      fill_colors = fill_colors
    )
  )
}

normalize_border_color <- function(border_color = NA) {
  if (is.null(border_color) || length(border_color) == 0L) {
    return(NA_character_)
  }
  as.character(border_color)[1L]
}
