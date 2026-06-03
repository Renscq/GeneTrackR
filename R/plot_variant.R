# Author: Rensc
# Date: 2026-05-29
# Version: 0.3.0
# Function: Plot variant records from VariantTrack objects
# Input: VariantTrack objects and genomic regions
# Output: ggplot variant figures

#' Plot variants in a genomic region
#'
#' @description
#' Draws VCF-derived variant records in a genomic region. This function is
#' dedicated to variant visualization and can be used directly or through
#' `plot_tracks()`.
#'
#' @param variant A VariantTrack object.
#' @param chrom Chromosome name. If NULL, all chromosomes in `variant` are used.
#' @param start Region start in 1-based coordinates. If NULL, the minimum variant position is used.
#' @param end Region end in 1-based coordinates. If NULL, the maximum variant position is used.
#' @param variant_type Optional variant type filter, such as `SNP`, `INS`, `DEL`, or `MNV`.
#' @param color_by Column used for variant colors. Use `none` to disable grouping.
#' @param label_by Column used for labels. Use `none` to hide labels.
#' @param variant_shape Plot geometry. `lollipop` draws vertical stems and points,
#'   `point` draws points only, and `rug` draws bottom ticks.
#' @param variant_palette RColorBrewer palette name for variant colors.
#' @param variant_colors Optional named or unnamed variant color vector.
#' @param point_size Point size.
#' @param line_width Line width for lollipop/rug stems.
#' @param text_size Text size.
#' @return A ggplot object.
#' @export
plot_variant <- function(variant,
                         chrom = NULL,
                         start = NULL,
                         end = NULL,
                         variant_type = NULL,
                         color_by = c("variant_type", "filter", "none"),
                         label_by = c("none", "variant_id", "variant_type", "ref", "alt"),
                         variant_shape = c("lollipop", "point", "rug"),
                         variant_palette = "Set1",
                         variant_colors = NULL,
                         point_size = 2,
                         line_width = 0.35,
                         text_size = 14) {
  stop_if_not(inherits(variant, "VariantTrack"), "`variant` must be a VariantTrack object.")

  color_by <- match.arg(color_by)
  label_by <- match.arg(label_by)
  variant_shape <- match.arg(variant_shape)
  text_color <- "black"
  track_name <- NULL

  vt <- retrieve_vcf(
    variant,
    chrom = chrom,
    start = start,
    end = end,
    variant_type = variant_type,
    as = "data.table"
  )

  if (nrow(vt) > 0L) {
    if (is.null(chrom)) chrom <- unique(vt[["chrom"]])
    if (is.null(start)) start <- min(vt[["pos"]], na.rm = TRUE)
    if (is.null(end)) end <- max(vt[["pos"]], na.rm = TRUE)
  } else {
    if (is.null(start)) start <- 0L
    if (is.null(end)) end <- start + 1L
  }

  start <- as.integer(start)[1L]
  end <- as.integer(end)[1L]
  stop_if_not(is.finite(start) && is.finite(end) && start <= end, "`start` and `end` must define a valid region.")

  chrom_label <- if (length(chrom) == 1L && !is.na(chrom)) as.character(chrom) else "selected region"
  x_label <- if (length(chrom) == 1L && !is.na(chrom)) {
    paste0("Chromosome ", chrom_label, " position (bp)")
  } else {
    "Genomic position (bp)"
  }

  if (nrow(vt) == 0L) {
    empty_dt <- data.table::data.table(pos = c(start, end), y = c(1, 1))
    return(
      ggplot2::ggplot(empty_dt, ggplot2::aes(x = .data$pos, y = .data$y)) +
        ggplot2::geom_blank() +
        ggplot2::coord_cartesian(xlim = c(start, end)) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::labs(x = x_label, y = track_name %||% "Variant") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
          axis.title.x = ggplot2::element_text(color = text_color, size = text_size),
          axis.title.y = ggplot2::element_text(color = text_color, size = text_size),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank()
        )
    )
  }

  vt[, "y" := 1]
  vt[, "label" := if (label_by == "none") "" else as.character(.SD[[label_by]])]
  if (color_by == "none") {
    vt[, "color_group" := "variant"]
  } else {
    stop_if_not(color_by %in% names(vt), paste0("Column not found for `color_by`: ", color_by))
    vt[, "color_group" := as.character(.SD[[color_by]])]
    vt[is.na(color_group) | !nzchar(color_group), "color_group" := "NA"]
  }

  p <- ggplot2::ggplot(vt, ggplot2::aes(x = .data$pos, y = .data$y, color = .data$color_group))

  if (variant_shape == "lollipop") {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(ymin = 0.72, ymax = 1.18),
      linewidth = line_width
    ) +
      ggplot2::geom_point(size = point_size)
  } else if (variant_shape == "point") {
    p <- p + ggplot2::geom_point(size = point_size)
  } else {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(ymin = 0.72, ymax = 1.08),
      linewidth = line_width
    )
  }

  p <- p +
    ggplot2::coord_cartesian(xlim = c(start, end)) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      expand = ggplot2::expansion(mult = c(0.20, 0.20), add = c(0.20, 0.20))
    ) +
    ggplot2::labs(x = x_label, y = track_name %||% "Variant") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.y = ggplot2::element_text(color = text_color, size = text_size),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = text_color, size = text_size),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  color_values <- make_variant_palette(unique(vt[["color_group"]]), palette = variant_palette, colors = variant_colors)
  p <- p + ggplot2::scale_color_manual(values = color_values)

  if (color_by == "none") {
    p <- p + ggplot2::guides(color = "none")
  }

  if (label_by != "none") {
    p <- p + ggplot2::geom_text(
      data = vt[nzchar(label)],
      ggplot2::aes(label = .data$label),
      inherit.aes = TRUE,
      size = text_size / 2.845276,
      color = text_color,
      vjust = -0.8,
      show.legend = FALSE
    )
  }

  p
}

#' Plot a VariantTrack object
#'
#' @description Backward-compatible alias of `plot_variant()`.
#' @inheritParams plot_variant
#' @return A ggplot object.
#' @export
plot_variant_track <- function(track,
                               chrom,
                               start,
                               end,
                               color_by = c("variant_type", "filter", "none"),
                               label_by = c("none", "variant_id", "variant_type", "ref", "alt"),
                               variant_palette = "Set1",
                               variant_colors = NULL,
                               text_size = 14) {
  plot_variant(
    variant = track,
    chrom = chrom,
    start = start,
    end = end,
    color_by = color_by,
    label_by = label_by,
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    text_size = text_size
  )
}

make_variant_palette <- function(keys, palette = "Set2", colors = NULL) {
  keys <- as.character(keys)
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0L) return(character())

  if (!is.null(colors)) {
    colors <- as.character(colors)
    if (!is.null(names(colors)) && all(keys %in% names(colors))) {
      return(colors[keys])
    }
    if (length(colors) >= length(keys)) {
      out <- colors[seq_along(keys)]
      names(out) <- keys
      return(out)
    }
    out <- grDevices::colorRampPalette(colors)(length(keys))
    names(out) <- keys
    return(out)
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE) && palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_col <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    min_col <- max(3L, min(max_col, length(keys)))
    base <- RColorBrewer::brewer.pal(min_col, palette)
    out <- if (length(keys) > length(base)) grDevices::colorRampPalette(base)(length(keys)) else base[seq_along(keys)]
  } else {
    out <- grDevices::hcl.colors(length(keys), palette = "Dynamic")
  }

  names(out) <- keys
  out
}
