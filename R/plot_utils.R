# Author: Rensc
# Date: 2026-05-27
# Version: dev003
# Function: Internal plotting utilities for gene models and signal tracks
# Input: Annotation and signal tables
# Output: ggplot objects and transformed plotting tables

normalize_plot_theme <- function(plot_theme = c("bw", "classic", "light", "minimal")) {
  match.arg(plot_theme)
}

normalize_show_panel_border <- function(show_panel_border = NULL) {
  if (is.null(show_panel_border)) {
    return(NULL)
  }
  if (!is.logical(show_panel_border) ||
      length(show_panel_border) != 1L ||
      is.na(show_panel_border)) {
    stop("`show_panel_border` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  show_panel_border
}

make_track_theme <- function(
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
) {
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)

  out <- switch(
    plot_theme,
    bw = ggplot2::theme_bw(),
    classic = ggplot2::theme_classic(),
    light = ggplot2::theme_light(),
    minimal = ggplot2::theme_minimal()
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black")
    )

  if (isTRUE(show_panel_border)) {
    out <- out + ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  } else if (identical(show_panel_border, FALSE)) {
    out <- out + ggplot2::theme(panel.border = ggplot2::element_blank())
  }

  out
}

normalize_signal_alpha <- function(signal_alpha = 0.85) {
  signal_alpha <- suppressWarnings(as.numeric(signal_alpha)[1L])
  if (!is.finite(signal_alpha) || signal_alpha < 0 || signal_alpha > 1) {
    stop("`signal_alpha` must be a numeric value between 0 and 1.", call. = FALSE)
  }
  signal_alpha
}

normalize_track_height <- function(value, name, default) {
  value <- suppressWarnings(as.numeric(value)[1L])
  if (!is.finite(value) || value <= 0) {
    warning(
      sprintf("`%s` must be a positive numeric value. Falling back to %s.", name, default),
      call. = FALSE
    )
    return(as.numeric(default))
  }
  value
}

normalize_signal_bar_width <- function(signal_bar_width = 1) {
  signal_bar_width <- suppressWarnings(as.numeric(signal_bar_width)[1L])
  if (!is.finite(signal_bar_width) ||
      signal_bar_width <= 0 ||
      signal_bar_width > 1) {
    stop(
      "`signal_bar_width` must be a numeric value greater than 0 and no greater than 1.",
      call. = FALSE
    )
  }
  signal_bar_width
}

normalize_signal_y_limits <- function(signal_y_limits = NULL) {
  if (is.null(signal_y_limits)) {
    return(NULL)
  }
  signal_y_limits <- suppressWarnings(as.numeric(signal_y_limits))
  if (length(signal_y_limits) != 2L ||
      any(!is.finite(signal_y_limits)) ||
      signal_y_limits[1L] >= signal_y_limits[2L]) {
    stop(
      "`signal_y_limits` must contain two finite increasing numeric values.",
      call. = FALSE
    )
  }
  signal_y_limits
}

resolve_signal_y_scale <- function(signal_y_scale, signal_y_limits = NULL) {
  signal_y_scale <- match.arg(signal_y_scale, c("free", "fixed"))
  if (!is.null(signal_y_limits) && signal_y_scale == "free") {
    warning(
      "`signal_y_limits` requires a shared y-axis; `signal_y_scale` was changed to 'fixed'.",
      call. = FALSE
    )
    signal_y_scale <- "fixed"
  }
  signal_y_scale
}

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


# Stable gene-model feature levels used by plot_gene(), plot_transcript(),
# plot_region(), and haplotype/LD gene tracks. Keeping these levels fixed avoids
# color shifts when a plot contains only CDS, only UTR, or only exon segments.
gene_model_feature_levels <- function() {
  c("UTR", "CDS", "exon")
}

normalize_gene_model_feature <- function(feature) {
  x <- as.character(feature)
  xl <- tolower(x)
  out <- x
  out[xl %in% c("cds")] <- "CDS"
  out[xl %in% c(
    "utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr",
    "five_utr", "three_utr", "5'utr", "3'utr", "five prime utr",
    "three prime utr", "five-prime-utr", "three-prime-utr"
  )] <- "UTR"
  out[xl %in% c("exon")] <- "exon"
  out
}

make_gene_model_fill_colors <- function(color_palette = "Paired", gene_colors = NULL, features = NULL) {
  canonical <- gene_model_feature_levels()
  feature_values <- normalize_gene_model_feature(features %||% canonical)
  feature_values <- unique(as.character(feature_values))
  feature_values <- feature_values[!is.na(feature_values) & nzchar(feature_values)]
  all_levels <- unique(c(canonical, feature_values))

  default_cols <- normalize_discrete_fill_colors(
    n = length(canonical),
    color_palette = color_palette,
    fill_colors = NULL
  )
  names(default_cols) <- canonical

  if (length(all_levels) > length(canonical)) {
    extra_levels <- setdiff(all_levels, canonical)
    full_cols <- normalize_discrete_fill_colors(
      n = length(all_levels),
      color_palette = color_palette,
      fill_colors = NULL
    )
    extra_cols <- full_cols[(length(canonical) + 1L):length(all_levels)]
    names(extra_cols) <- extra_levels
    default_cols <- c(default_cols, extra_cols)
  }

  if (is.null(gene_colors)) {
    return(default_cols[all_levels])
  }

  gene_colors <- as.character(gene_colors)
  gene_colors <- gene_colors[!is.na(gene_colors) & nzchar(gene_colors)]
  if (length(gene_colors) == 0L) {
    return(default_cols[all_levels])
  }

  color_names <- names(gene_colors)
  if (!is.null(color_names) && any(nzchar(color_names))) {
    color_names <- normalize_gene_model_feature(color_names)
    names(gene_colors) <- color_names
    named_colors <- gene_colors[nzchar(names(gene_colors))]
    # If aliases such as five_prime_UTR and three_prime_UTR are both supplied,
    # keep the first explicit color after normalization to the shared UTR level.
    named_colors <- named_colors[!duplicated(names(named_colors))]
    out <- default_cols
    common <- intersect(names(named_colors), names(out))
    out[common] <- named_colors[common]
    extra_named <- named_colors[setdiff(names(named_colors), names(out))]
    if (length(extra_named) > 0L) {
      out <- c(out, extra_named)
    }
    return(out[unique(c(all_levels, intersect(names(out), names(named_colors))))])
  }

  unnamed_cols <- normalize_discrete_fill_colors(
    n = length(canonical),
    color_palette = color_palette,
    fill_colors = gene_colors
  )
  names(unnamed_cols) <- canonical
  out <- default_cols
  out[canonical] <- unnamed_cols[canonical]
  out[all_levels]
}

apply_gene_model_fill_scale <- function(p,
                                        color_palette = "Paired",
                                        gene_colors = NULL,
                                        features = NULL,
                                        drop = TRUE) {
  colors <- make_gene_model_fill_colors(
    color_palette = color_palette,
    gene_colors = gene_colors,
    features = features
  )
  observed <- normalize_gene_model_feature(features %||% names(colors))
  observed <- unique(as.character(observed))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  breaks <- intersect(names(colors), observed)
  if (length(breaks) == 0L) {
    breaks <- intersect(names(colors), gene_model_feature_levels())
  }

  p + ggplot2::scale_fill_manual(
    values = colors,
    limits = names(colors),
    breaks = breaks,
    drop = drop,
    na.value = "grey80"
  )
}


# Stable variant marker levels used by plot_hap_variant() and plot_ld_block().
# They keep marker colors fixed across figures even when a region contains only
# SNPs, only indels, or additional/unknown variant classes.
variant_marker_levels <- function() {
  c("SNP", "Ind", "...")
}

normalize_variant_marker_type <- function(variant_type) {
  x <- as.character(variant_type)
  xl <- tolower(x)
  out <- rep("...", length(x))
  out[xl %in% c("snp", "snv")] <- "SNP"
  out[xl %in% c(
    "ind", "indel", "ins", "del", "mnv", "insertion", "deletion",
    "insert", "delete"
  )] <- "Ind"
  out[is.na(x) | !nzchar(x)] <- "..."
  out
}

make_variant_marker_fill_colors <- function(variant_palette = "Paired",
                                            variant_colors = NULL,
                                            variant_types = NULL,
                                            alpha = 0.6) {
  canonical <- variant_marker_levels()
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  if (!is.finite(alpha)) alpha <- 0.6
  alpha <- max(0, min(1, alpha))

  observed <- normalize_variant_marker_type(variant_types %||% canonical)
  observed <- unique(as.character(observed))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  all_levels <- unique(c(canonical, observed))

  default_cols <- normalize_discrete_fill_colors(
    n = length(canonical),
    color_palette = variant_palette,
    fill_colors = NULL
  )
  names(default_cols) <- canonical

  if (length(all_levels) > length(canonical)) {
    extra_levels <- setdiff(all_levels, canonical)
    full_cols <- normalize_discrete_fill_colors(
      n = length(all_levels),
      color_palette = variant_palette,
      fill_colors = NULL
    )
    extra_cols <- full_cols[(length(canonical) + 1L):length(all_levels)]
    names(extra_cols) <- extra_levels
    default_cols <- c(default_cols, extra_cols)
  }

  if (!is.null(variant_colors)) {
    variant_colors <- as.character(variant_colors)
    variant_colors <- variant_colors[!is.na(variant_colors) & nzchar(variant_colors)]
    if (length(variant_colors) > 0L) {
      color_names <- names(variant_colors)
      if (!is.null(color_names) && any(nzchar(color_names))) {
        color_names <- normalize_variant_marker_type(color_names)
        names(variant_colors) <- color_names
        named_colors <- variant_colors[nzchar(names(variant_colors))]
        named_colors <- named_colors[!duplicated(names(named_colors))]
        common <- intersect(names(named_colors), names(default_cols))
        default_cols[common] <- named_colors[common]
        extra_named <- named_colors[setdiff(names(named_colors), names(default_cols))]
        if (length(extra_named) > 0L) {
          default_cols <- c(default_cols, extra_named)
        }
      } else {
        unnamed_cols <- normalize_discrete_fill_colors(
          n = length(canonical),
          color_palette = variant_palette,
          fill_colors = variant_colors
        )
        names(unnamed_cols) <- canonical
        default_cols[canonical] <- unnamed_cols[canonical]
      }
    }
  }

  out <- default_cols[unique(c(all_levels, intersect(names(default_cols), canonical)))]
  grDevices::adjustcolor(out, alpha.f = alpha)
}


resolve_variant_marker_colors <- function(variant_types,
                                          variant_palette = "Paired",
                                          variant_colors = NULL,
                                          alpha = 0.6) {
  labels <- normalize_variant_marker_type(variant_types)
  colors <- make_variant_marker_fill_colors(
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    variant_types = labels,
    alpha = alpha
  )
  out <- colors[labels]
  missing <- is.na(out) | !nzchar(out)
  if (any(missing)) {
    out[missing] <- grDevices::adjustcolor("grey85", alpha.f = alpha)
  }
  unname(out)
}

apply_variant_marker_color_scale <- function(p,
                                             variant_types,
                                             variant_palette = "Paired",
                                             variant_colors = NULL,
                                             alpha = 0.6,
                                             name = "Variant type",
                                             drop = FALSE) {
  colors <- make_variant_marker_fill_colors(
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    variant_types = variant_types,
    alpha = alpha
  )
  observed <- normalize_variant_marker_type(variant_types %||% names(colors))
  observed <- unique(as.character(observed))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  breaks <- intersect(variant_marker_levels(), observed)
  if (length(breaks) == 0L) {
    breaks <- variant_marker_levels()
  }
  legend_cols <- colors[breaks]
  legend_cols <- legend_cols[!is.na(legend_cols)]

  p +
    ggplot2::scale_color_manual(
      values = colors,
      limits = names(colors),
      breaks = breaks,
      name = name,
      drop = drop,
      na.value = grDevices::adjustcolor("grey85", alpha.f = alpha)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          shape = 17,
          color = unname(legend_cols),
          alpha = 1
        )
      )
    )
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
