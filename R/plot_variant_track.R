# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.2
# Function: Plot genome-level VariantTrack objects
# Input: VariantTrack objects and genomic regions
# Output: ggplot variant track figures

#' Plot a VariantTrack object
#'
#' @description
#' Draws VCF-derived variant positions in a genomic region. This function is
#' designed to be compatible with `plot_tracks()`.
#'
#' @param track A VariantTrack object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param color_by Column used for variant colors. Usually `variant_type` or `filter`.
#' @param label_by Column used for labels. Use `none` to hide labels.
#' @param color_palette RColorBrewer palette name.
#' @param fill_colors Optional named or unnamed color vector.
#' @param text_size Text size.
#' @param text_color Text color.
#' @param track_name Optional y-axis label.
#' @return A ggplot object.
#' @export
plot_variant_track <- function(track,
                               chrom,
                               start,
                               end,
                               color_by = c("variant_type", "filter"),
                               label_by = c("none", "variant_id", "variant_type"),
                               color_palette = "Set1",
                               fill_colors = NULL,
                               text_size = 14,
                               text_color = "black",
                               track_name = NULL) {
  stop_if_not(inherits(track, "VariantTrack"), "`track` must be a VariantTrack object.")
  color_by <- match.arg(color_by)
  label_by <- match.arg(label_by)
  vt <- slice_variant_track(track, chrom = chrom, start = start, end = end)$data
  x_label <- paste0("Chromosome ", as.character(chrom)[1L], " position (bp)")

  if (nrow(vt) == 0L) {
    vt <- data.table::data.table(pos = c(start, end), y = c(1, 1), fill_group = "empty", label = "")
    p <- ggplot2::ggplot(vt, ggplot2::aes(x = .data$pos, y = .data$y)) +
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
    return(p)
  }

  vt[, "y" := 1]
  vt[, "fill_group" := as.character(vt[[color_by]])]
  vt[, "label" := if (label_by == "none") "" else as.character(vt[[label_by]])]

  p <- ggplot2::ggplot(vt, ggplot2::aes(x = .data$pos, y = .data$y, color = .data$fill_group)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = 0.75, ymax = 1.25), linewidth = 0.35) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_cartesian(xlim = c(start, end)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = c(0.2, 0.2), add = c(0.2, 0.2))) +
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
  color_values <- make_extended_palette(unique(vt[["fill_group"]]), palette = color_palette, colors = fill_colors)
  p <- p + ggplot2::scale_color_manual(values = color_values)

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

make_extended_palette <- function(keys, palette = "Set2", colors = NULL) {
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
