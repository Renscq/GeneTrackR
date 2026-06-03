# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.0
# Function: Plot generic feature and variant tracks
# Input: FeatureTrack or VariantTrack objects
# Output: ggplot track panels

#' Plot a FeatureTrack object
#'
#' @description
#' Draws BED/GFF/GTF-derived interval features in a genomic region. This function
#' is designed to be compatible with `plot_tracks()`.
#'
#' @param track A FeatureTrack object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param mode `overlap`, `within`, or `trim`.
#' @param color_by Column used for fill colors. Usually `type`, `source`, or `name`.
#' @param label_by Column used for labels. Use `none` to hide labels.
#' @param feature_palette RColorBrewer palette name for feature fills.
#' @param feature_colors Optional named or unnamed feature fill color vector.
#' @param feature_border_color Rectangle border color. Use NA to hide borders.
#' @param text_size Text size.
#' @return A ggplot object.
#' @export
plot_feature_track <- function(track,
                               chrom,
                               start,
                               end,
                               mode = c("overlap", "within", "trim"),
                               color_by = c("type", "source", "name"),
                               label_by = c("none", "name", "feature_id", "type"),
                               feature_palette = "Set2",
                               feature_colors = NULL,
                               feature_border_color = NA,
                               text_size = 14) {
  stop_if_not(inherits(track, "FeatureTrack"), "`track` must be a FeatureTrack object.")
  mode <- match.arg(mode)
  color_by <- match.arg(color_by)
  label_by <- match.arg(label_by)
  text_color <- "black"
  track_name <- NULL
  feature_border_color <- normalize_border_color(feature_border_color)

  ft <- retrieve_feature(
    track,
    chrom = chrom,
    start = start,
    end = end,
    mode = mode,
    level = "feature",
    as = "data.table"
  )
  if (is.null(ft)) {
    ft <- data.table::data.table()
  } else {
    ft <- data.table::as.data.table(ft)
  }
  x_label <- paste0("Chromosome ", as.character(chrom)[1L], " position (bp)")

  if (NROW(ft) == 0L) {
    ft <- data.table::data.table(
      feature_id = "empty_feature_track",
      name = "",
      chrom = as.character(chrom)[1L],
      start = as.integer(start)[1L],
      end = as.integer(end)[1L],
      type = "empty",
      source = "empty",
      y = 1L,
      ymin = 0.75,
      ymax = 1.25,
      fill_group = "empty",
      label = ""
    )
    p <- ggplot2::ggplot(ft) +
      ggplot2::geom_blank(ggplot2::aes(x = .data$start, y = .data$y)) +
      ggplot2::coord_cartesian(xlim = c(start, end)) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::labs(x = x_label, y = track_name %||% "Feature") +
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

  data.table::setorderv(ft, c("start", "end", "feature_id"))
  ft[, "y" := seq_len(.N)]
  ft[, "ymin" := as.numeric(ft[["y"]]) - 0.30]
  ft[, "ymax" := as.numeric(ft[["y"]]) + 0.30]
  ft[, "fill_group" := as.character(ft[[color_by]])]
  ft[, "label" := if (label_by == "none") "" else as.character(ft[[label_by]])]

  p <- ggplot2::ggplot(ft) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$fill_group
      ),
      color = feature_border_color,
      linewidth = 0.2
    ) +
    ggplot2::coord_cartesian(xlim = c(start, end)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = c(0.2, 0.2), add = c(0.2, 0.2))) +
    ggplot2::labs(x = x_label, y = track_name %||% "Feature") +
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
  p <- apply_discrete_fill_scale(p, color_palette = feature_palette, fill_colors = feature_colors)

  if (label_by != "none") {
    ft[, "label_x" := (as.numeric(ft[["start"]]) + as.numeric(ft[["end"]])) / 2]
    p <- p + ggplot2::geom_text(
      data = ft[nzchar(label)],
      ggplot2::aes(x = .data$label_x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      size = text_size / 2.845276,
      color = text_color,
      vjust = -0.6
    )
  }
  p
}

