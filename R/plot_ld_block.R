# Author: Rensc
# Date: 2026-07-02
# Version: dev002
# Function: Plot LD block triangular heatmaps
# Input: LDTrack objects, VariantTrack objects, or VCF paths
# Output: LDTrack objects with stored figures

#' Plot a triangular LD heatmap
#'
#' @description
#' Draws an inverted triangular LD heatmap using ggplot2. For exactly two
#' variants, the single pairwise LD value is drawn as one rectangular heatmap
#' cell. Interior grid lines are suppressed; only the outside frame is drawn. The function accepts
#' either an `LDTrack` object or a `VariantTrack`/VCF path, in which case LD is
#' computed first by `compute_ld_block()`. The default return value is an
#' updated `LDTrack` object: LD calculation results remain in the object and the
#' generated plot is stored in `LDTrack$figure`. When `show_region = TRUE`, the region
#' track follows the compact gene-model style used by `plot_hap_variant()` in
#' GeneTrackR 0.3.15 and connects genomic variant markers to LD heatmap columns
#' with a shared x scale, so line endpoints remain aligned after resizing.
#'
#' @param object An `LDTrack` object, a `VariantTrack` object, or a VCF path.
#' @param chrom,start,end Optional region used when `object` is not an `LDTrack`.
#' @param variant_type One of `both`, `snp`, or `ind` when computing from VCF.
#' @param method LD method used when computing from VCF.
#' @param color_palette RColorBrewer sequential palette name. Default `Reds`.
#' @param font Base font size. Font color is always black.
#' @param title Plot title. Default is `chrom:start-end LD`.
#' @param label_by Variant label column, either `pos`, `variant_id`, or another
#' column in `LDTrack$variants`.
#' @param show_variant_labels Logical. Whether to show labels above the heatmap.
#' @param show_region Logical. Whether to add a compact region/variant track
#' above the LD heatmap.
#' @param show_variant_marker Logical. Whether to draw natural-variant triangle markers in the region track.
#' @param variant_marker_size Size of natural-variant triangle markers. Set to 0 to hide markers.
#' @param annotation Optional gene annotation used by the compact region track.
#' @param connect_region Logical. Whether to connect region-track variant
#' positions to heatmap columns.
#' @param region_height Relative height of the compact region track.
#' @param connector_height Relative height of connector lines.
#' @param heatmap_height Relative height of the heatmap. If NULL, it is chosen
#' from the number of variants.
#' @param samples,min_pair_samples,ploidy,verbose Passed to `compute_ld_block()`
#' when `object` is not already an `LDTrack`.
#' @param return_object Logical. Default TRUE. If TRUE, return an updated
#' `LDTrack` containing the figure in `$figure`; if FALSE, return the figure
#' directly for compatibility with older scripts.
#' @return By default, an updated `LDTrack` object with the generated figure
#' stored in `$figure`. If `return_object = FALSE`, the ggplot/patchwork figure
#' is returned directly.
#' @export
plot_ld_block <- function(object,
                          chrom = NULL,
                          start = NULL,
                          end = NULL,
                          variant_type = c("both", "snp", "ind"),
                          method = c("r2", "Dprime"),
                          color_palette = "Reds",
                          font = 14,
                          title = NULL,
                          label_by = c("pos", "variant_id"),
                          show_variant_labels = TRUE,
                          show_region = FALSE,
                          show_variant_marker = TRUE,
                          variant_marker_size = 2.8,
                          annotation = NULL,
                          connect_region = TRUE,
                          region_height = 1.25,
                          connector_height = 0.35,
                          heatmap_height = NULL,
                          samples = NULL,
                          min_pair_samples = 3L,
                          ploidy = 2L,
                          verbose = TRUE,
                          return_object = TRUE) {
  method <- match.arg(method)
  variant_type <- match.arg(variant_type)
  show_variant_marker <- isTRUE(show_variant_marker)
  variant_marker_size <- as.numeric(variant_marker_size)[1L]
  if (is.na(variant_marker_size)) variant_marker_size <- 2.8
  if (variant_marker_size <= 0) {
    show_variant_marker <- FALSE
    variant_marker_size <- 0
  }

  if (inherits(object, "LDTrack")) {
    ld <- object
  } else {
    ld <- compute_ld_block(
      object,
      chrom = chrom,
      start = start,
      end = end,
      variant_type = variant_type,
      method = method,
      samples = samples,
      min_pair_samples = min_pair_samples,
      ploidy = ploidy,
      verbose = verbose
    )
  }

  variants <- data.table::copy(data.table::as.data.table(ld$variants))
  pair_dt <- data.table::copy(data.table::as.data.table(ld$data))
  stop_if_not(nrow(variants) >= 2L, "At least two variants are required for LD plotting.")
  stop_if_not(nrow(pair_dt) > 0L, "No pairwise LD values are available for plotting.")

  if (!"variant_index" %in% names(variants)) {
    variants[, "variant_index" := seq_len(.N)]
  }
  data.table::setorderv(variants, "variant_index")

  label_by <- as.character(label_by)[1L]
  if (is.na(label_by) || !label_by %in% names(variants)) {
    label_by <- "pos"
  }
  variants[, "variant_label" := as.character(.SD[[1L]]), .SDcols = label_by]
  variants[is.na(variant_label) | variant_label == "", "variant_label" := as.character(variant_id)]

  plot_method <- ld$meta$method %||% method
  if (is.null(title)) {
    region <- ld$region %||% list()
    if (!is.null(region$chrom) && !is.null(region$start) && !is.null(region$end)) {
      title <- paste0(region$chrom, ":", region$start, "-", region$end, " LD (", plot_method, ")")
    } else {
      title <- paste0("LD (", plot_method, ")")
    }
  }

  x_limits <- c(0.5, nrow(variants) + 0.5)
  p_heat <- draw_ld_triangle_heatmap(
    ld = ld,
    variants = variants,
    pair_dt = pair_dt,
    color_palette = color_palette,
    font = font,
    title = title,
    show_variant_labels = show_variant_labels,
    x_limits = x_limits
  )

  if (is.null(heatmap_height)) {
    heatmap_height <- max(2.0, min(8.0, nrow(variants) * 0.06))
  }

  gene_data <- NULL
  if (isTRUE(show_region)) {
    region_mapper <- make_ld_region_mapper(variants, ld$region %||% list())
    variants[, "region_x" := region_mapper(pos)]
    variants[, "gene_x" := region_x]
    variants[, "ld_x" := as.numeric(variant_index)]

    p_region <- draw_ld_region_track(
      ld = ld,
      variants = variants,
      annotation = annotation,
      x_limits = x_limits,
      font = font,
      show_variant_marker = show_variant_marker,
      variant_marker_size = variant_marker_size
    )
    gene_data <- attr(p_region, "gene_data")

    if (isTRUE(connect_region)) {
      p_connector <- draw_ld_connector_track(variants = variants, x_limits = x_limits)
      p <- patchwork::wrap_plots(
        p_region,
        p_connector,
        p_heat,
        ncol = 1,
        heights = c(region_height, connector_height, heatmap_height)
      )
    } else {
      p <- patchwork::wrap_plots(
        p_region,
        p_heat,
        ncol = 1,
        heights = c(region_height, heatmap_height)
      )
    }
  } else {
    p <- p_heat
  }

  attr(p, "ld_data") <- pair_dt[]
  attr(p, "variant_data") <- variants[]
  attr(p, "gene_data") <- gene_data
  ld$figure <- p
  if (isTRUE(return_object)) return(ld)
  p
}

draw_ld_triangle_heatmap <- function(ld,
                                     variants,
                                     pair_dt,
                                     color_palette = "Reds",
                                     font = 14,
                                     title = NULL,
                                     show_variant_labels = TRUE,
                                     x_limits = NULL) {
  n_var <- nrow(variants)
  poly_dt <- make_ld_triangle_polygons(pair_dt, n_var = n_var)
  frame_dt <- make_ld_triangle_frame(n_var)
  colors <- make_ld_continuous_palette(color_palette)
  legend_title <- if (identical(ld$meta$method %||% "r2", "r2")) expression(r^2) else "D'"
  if (is.null(x_limits)) x_limits <- c(0.5, n_var + 0.5)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_dt,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$pair_id, fill = .data$ld),
      color = NA
    ) +
    ggplot2::geom_path(
      data = frame_dt,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_gradientn(
      colors = colors,
      limits = c(0, 1),
      na.value = "grey95",
      name = legend_title
    ) +
    ggplot2::scale_x_continuous(
      breaks = variants[["variant_index"]],
      labels = variants[["variant_label"]],
      limits = x_limits,
      expand = c(0, 0),
      position = "top"
    ) +
    ggplot2::scale_y_reverse(
      limits = c(max(frame_dt[["y"]], na.rm = TRUE) + 0.05, -0.05),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = font, color = "black"),
      plot.title = ggplot2::element_text(size = font, color = "black", hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = font * 0.75, color = "black", angle = 90, hjust = 0, vjust = 0.5),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = font * 0.85, color = "black"),
      legend.text = ggplot2::element_text(size = font * 0.75, color = "black"),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  if (!isTRUE(show_variant_labels)) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  p
}

make_ld_triangle_polygons <- function(pair_dt, n_var = NULL) {
  x <- data.table::copy(data.table::as.data.table(pair_dt))
  x[, "center_x" := (as.numeric(index_i) + as.numeric(index_j)) / 2]
  x[, "center_y" := (as.numeric(index_j) - as.numeric(index_i)) / 2]
  x[, "pair_id" := paste0(index_i, "_", index_j)]

  if (!is.null(n_var) && as.integer(n_var)[1L] == 2L) {
    stop_if_not(nrow(x) == 1L, "Exactly one LD pair is required when plotting two variants.")
    return(data.table::rbindlist(list(
      x[, .(pair_id, ld, x = center_x - 0.5, y = 0)],
      x[, .(pair_id, ld, x = center_x + 0.5, y = 0)],
      x[, .(pair_id, ld, x = center_x + 0.5, y = 1)],
      x[, .(pair_id, ld, x = center_x - 0.5, y = 1)]
    ), fill = TRUE)[])
  }

  data.table::rbindlist(list(
    x[, .(pair_id, ld, x = center_x, y = center_y - 0.5)],
    x[, .(pair_id, ld, x = center_x + 0.5, y = center_y)],
    x[, .(pair_id, ld, x = center_x, y = center_y + 0.5)],
    x[, .(pair_id, ld, x = center_x - 0.5, y = center_y)]
  ), fill = TRUE)[order(pair_id)][]
}

make_ld_triangle_frame <- function(n_var) {
  if (n_var == 2L) {
    return(data.table::data.table(
      x = c(1, 2, 2, 1, 1),
      y = c(0, 0, 1, 1, 0)
    ))
  }
  data.table::data.table(
    x = c(1, 1.5, n_var - 0.5, n_var, (n_var + 1) / 2, 1),
    y = c(0.5, 0, 0, 0.5, n_var / 2, 0.5)
  )
}

make_ld_continuous_palette <- function(color_palette = "Reds") {
  color_palette <- as.character(color_palette)[1L]
  if (is.na(color_palette) || !nzchar(color_palette)) color_palette <- "Reds"
  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_info <- RColorBrewer::brewer.pal.info[color_palette, , drop = FALSE]
    max_colors <- as.integer(pal_info[["maxcolors"]][1L])
    base <- RColorBrewer::brewer.pal(max(3L, max_colors), color_palette)
    return(grDevices::colorRampPalette(base)(100L))
  }
  warning(sprintf("Unknown `color_palette`: %s. Falling back to 'Reds'.", color_palette), call. = FALSE)
  base <- RColorBrewer::brewer.pal(9L, "Reds")
  grDevices::colorRampPalette(base)(100L)
}

draw_ld_region_track <- function(ld,
                                 variants,
                                 annotation = NULL,
                                 x_limits = NULL,
                                 font = 14,
                                 show_variant_marker = TRUE,
                                 variant_marker_size = 2.8) {
  vars <- data.table::copy(data.table::as.data.table(variants))
  if (is.null(x_limits)) x_limits <- c(0.5, nrow(vars) + 0.5)
  show_variant_marker <- isTRUE(show_variant_marker)
  variant_marker_size <- as.numeric(variant_marker_size)[1L]
  if (is.na(variant_marker_size)) variant_marker_size <- 2.8
  if (variant_marker_size <= 0) {
    show_variant_marker <- FALSE
    variant_marker_size <- 0
  }
  if (!"region_x" %in% names(vars)) vars[, "region_x" := as.numeric(variant_index)]
  if (!"gene_x" %in% names(vars)) vars[, "gene_x" := region_x]

  gene_data <- NULL
  if (!is.null(annotation)) {
    gene_data <- tryCatch(
      prepare_hap_gene_track(
        annotation = annotation,
        hap = list(region = ld$region %||% list()),
        vars = vars,
        mapper = function(x) make_ld_region_mapper(vars, ld$region %||% list())(x)
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(gene_data) && nrow(gene_data$transcripts) > 0L && nrow(gene_data$segments) > 0L) {
    p <- draw_hap_gene_track(
      gene_data = gene_data,
      vars = vars,
      x_limits = x_limits,
      x_breaks = vars$variant_index,
      x_labels = vars$variant_label,
      text_size = font,
      gene_text_size = font,
      exon_height = 0.22,
      cds_height = 0.44,
      gene_palette = "Paired",
      gene_colors = NULL,
      gene_border_color = NA,
      variant_palette = "Set1",
      variant_colors = NULL,
      variant_alpha = 0.6,
      show_variant_marker = show_variant_marker,
      variant_marker_size = variant_marker_size,
      gene_track_legend_position = "right",
      show_gene_pos_axis = TRUE,
      gene_pos_axis_n = 5L,
      gene_pos_axis_label = NULL,
      gene_pos_x_angle = 0
    )
    attr(p, "gene_data") <- gene_data
    return(p)
  }

  if (!"variant_type" %in% names(vars)) {
    vars[, "variant_type_label" := "..."]
  } else {
    vars[, "variant_type_label" := normalize_variant_marker_type(variant_type)]
  }

  p <- ggplot2::ggplot(vars, ggplot2::aes(x = .data$gene_x, y = 1))
  if (isTRUE(show_variant_marker)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$variant_type_label),
        size = variant_marker_size,
        shape = 17
      )
    p <- apply_variant_marker_color_scale(
      p,
      variant_types = vars$variant_type_label,
      variant_palette = "Set1",
      variant_colors = NULL,
      alpha = 0.6,
      name = "Variant type",
      drop = FALSE
    )
  }
  p <- p +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0), position = "top") +
    ggplot2::scale_y_continuous(limits = c(0.7, 1.35), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = font, color = "black"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 12, 0, 8)
    )
  attr(p, "gene_data") <- NULL
  p
}

draw_ld_connector_track <- function(variants, x_limits = NULL) {
  vars <- data.table::copy(data.table::as.data.table(variants))
  if (!"gene_x" %in% names(vars)) vars[, "gene_x" := as.numeric(variant_index)]
  if (!"ld_x" %in% names(vars)) vars[, "ld_x" := as.numeric(variant_index)]
  if (is.null(x_limits)) x_limits <- c(0.5, nrow(vars) + 0.5)

  ggplot2::ggplot(vars) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$gene_x, xend = .data$ld_x, y = 1, yend = 0),
      linewidth = 0.35,
      color = "grey35",
      lineend = "butt"
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 12, 0, 8))
}

make_ld_region_mapper <- function(variants, region) {
  vars <- data.table::as.data.table(variants)
  pos <- sort(unique(as.numeric(vars[["pos"]])))
  n <- nrow(vars)
  region_start <- suppressWarnings(as.numeric(region$start %||% min(pos, na.rm = TRUE)))
  region_end <- suppressWarnings(as.numeric(region$end %||% max(pos, na.rm = TRUE)))
  if (!is.finite(region_start)) region_start <- min(pos, na.rm = TRUE)
  if (!is.finite(region_end)) region_end <- max(pos, na.rm = TRUE)
  if (region_start > min(pos, na.rm = TRUE)) region_start <- min(pos, na.rm = TRUE)
  if (region_end < max(pos, na.rm = TRUE)) region_end <- max(pos, na.rm = TRUE)
  if (!is.finite(region_start) || !is.finite(region_end) || region_end <= region_start) {
    region_start <- min(pos, na.rm = TRUE) - 1
    region_end <- max(pos, na.rm = TRUE) + 1
  }
  if (region_end <= region_start) region_end <- region_start + 1

  function(pos) {
    pos <- as.numeric(pos)
    0.5 + (pos - region_start) / (region_end - region_start) * max(1, n)
  }
}
