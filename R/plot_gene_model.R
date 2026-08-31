# Author: Rensc
# Date: 2026-08-31
# Version: dev002
# Function: Plot transcript, gene, and genomic region structures
# Input: GenePred objects and feature identifiers or regions
# Output: ggplot gene model figures

#' Plot a transcript structure
#'
#' @param object A GenePred object.
#' @param transcript_id Transcript ID.
#' @param coordinate transcript removes introns, genomic keeps genomic scale.
#' @param show_cds Whether to distinguish CDS and UTR segments.
#' @param cds_height Vertical thickness of CDS rectangles.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles.
#' @param highlight Optional data frame used to shade genomic or transcript intervals. It must contain `start` and `end` columns. Optional columns such as `label` or `group` are allowed but are not required. In `coordinate = "genomic"`, `start` and `end` are genomic positions; in `coordinate = "transcript"`, they are spliced transcript coordinates.
#' @param gene_palette RColorBrewer palette name used for discrete fills. If the number of discrete groups exceeds the palette maximum, colors are automatically interpolated.
#' @param gene_colors Optional custom fill colors for gene model features. Use a named vector such as `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")` to map colors explicitly. Unnamed colors are matched to the fixed gene-model levels `UTR`, `CDS`, and `exon`, so colors remain stable even when only one feature type is present.
#' @param gene_border_color Optional rectangle border color. Use `NA` to hide borders, or a color such as `"black"` or `"grey30"` to draw feature outlines.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw the panel border. `NULL` preserves the selected theme default.
#' @param label_position Where to draw feature labels. `axis` draws labels on the y axis and `feature` draws labels at the center of each gene/transcript structure.
#' @param label_by Which identifier to use for feature labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_size Text size in points for axis text, axis titles, and legends.
#' @param direction_mode Direction-arrow style. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @details
#' Common feature names used by `gene_colors` are `UTR`, `CDS`, and `exon`.
#' `highlight` must contain `start` and `end` columns. In genomic coordinate
#' plots these are genomic positions; in transcript coordinate plots they are
#' spliced transcript positions. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot object.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' plot_transcript(gp, transcript_id = "TxA1", coordinate = "transcript")
#' @export
plot_transcript <- function(object, transcript_id, coordinate = c("transcript", "genomic"), show_cds = TRUE, cds_height = 0.50, utr_height = 0.25, direction_mode = c("transcript", "gene", "end", "none"), highlight = NULL, gene_palette = "Paired", gene_colors = NULL, gene_border_color = NA, label_position = c("axis", "feature"), label_by = c("gene", "transcript"), text_size = 14, plot_theme = c("bw", "classic", "light", "minimal"), show_panel_border = NULL) {
  stop_if_not(is_gene_model_feature(object), "`object` must be a GenePred object or a Feature object with transcript/exon records.")
  object <- as_genepred(object)
  coordinate <- match.arg(coordinate)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  transcript_id_value <- as.character(transcript_id)

  tx_all <- data.table::as.data.table(object$transcripts)
  ex_all <- data.table::as.data.table(object$exons)
  tx <- tx_all[tx_all[["transcript_id"]] == transcript_id_value]
  stop_if_not(nrow(tx) > 0L, "Transcript ID was not found.")
  ex <- ex_all[ex_all[["transcript_id"]] == transcript_id_value]

  plot_gene_model_core(
    tx = tx,
    exons = ex,
    coordinate = coordinate,
    collapse = "none",
    show_cds = show_cds,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    highlight = highlight,
    gene_palette = gene_palette,
    gene_colors = gene_colors,
    gene_border_color = gene_border_color,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  )
}

#' Plot a gene structure
#'
#' @param object A GenePred object.
#' @param gene_id Gene ID.
#' @param collapse Transcript collapse mode.
#' @param coordinate genomic or transcript.
#' @param show_cds Whether to distinguish CDS and UTR segments.
#' @param cds_height Vertical thickness of CDS rectangles.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles.
#' @param highlight Optional data frame used to shade genomic or transcript intervals. It must contain `start` and `end` columns. Optional columns such as `label` or `group` are allowed but are not required. In `coordinate = "genomic"`, `start` and `end` are genomic positions; in `coordinate = "transcript"`, they are spliced transcript coordinates.
#' @param gene_palette RColorBrewer palette name used for discrete fills. If the number of discrete groups exceeds the palette maximum, colors are automatically interpolated.
#' @param gene_colors Optional custom fill colors for gene model features. Use a named vector such as `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")` to map colors explicitly. Unnamed colors are matched to the fixed gene-model levels `UTR`, `CDS`, and `exon`, so colors remain stable even when only one feature type is present.
#' @param gene_border_color Optional rectangle border color. Use `NA` to hide borders, or a color such as `"black"` or `"grey30"` to draw feature outlines.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw the panel border. `NULL` preserves the selected theme default.
#' @param label_position Where to draw feature labels. `axis` draws labels on the y axis and `feature` draws labels at the center of each gene/transcript structure.
#' @param label_by Which identifier to use for feature labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_size Text size in points for axis text, axis titles, and legends.
#' @param direction_mode Direction-arrow style. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @details
#' Common feature names used by `gene_colors` are `UTR`, `CDS`, and `exon`.
#' `highlight` must contain `start` and `end` columns. In genomic coordinate
#' plots these are genomic positions; in transcript coordinate plots they are
#' spliced transcript positions. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot object.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' plot_gene(gp, gene_id = "GeneA")
#' @export
plot_gene <- function(object, gene_id, collapse = c("none", "union_exon", "longest"), coordinate = c("genomic", "transcript"), show_cds = TRUE, cds_height = 0.50, utr_height = 0.25, direction_mode = c("transcript", "gene", "end", "none"), highlight = NULL, gene_palette = "Paired", gene_colors = NULL, gene_border_color = NA, label_position = c("axis", "feature"), label_by = c("gene", "transcript"), text_size = 14, plot_theme = c("bw", "classic", "light", "minimal"), show_panel_border = NULL) {
  stop_if_not(is_gene_model_feature(object), "`object` must be a GenePred object or a Feature object with transcript/exon records.")
  object <- as_genepred(object)
  collapse <- match.arg(collapse)
  coordinate <- match.arg(coordinate)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  direction_mode <- match.arg(direction_mode)
  gene_id_value <- as.character(gene_id)

  tx_all <- data.table::as.data.table(object$transcripts)
  ex_all <- data.table::as.data.table(object$exons)
  tx <- tx_all[tx_all[["gene_id"]] == gene_id_value]
  stop_if_not(nrow(tx) > 0L, "Gene ID was not found.")
  ex <- ex_all[ex_all[["gene_id"]] == gene_id_value]

  collapsed <- collapse_transcripts(tx, ex, collapse)
  plot_gene_model_core(
    tx = collapsed$transcripts,
    exons = collapsed$exons,
    coordinate = coordinate,
    collapse = collapse,
    show_cds = show_cds,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    highlight = highlight,
    gene_palette = gene_palette,
    gene_colors = gene_colors,
    gene_border_color = gene_border_color,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  )
}

#' Plot gene structures in a genomic region
#'
#' @param object A GenePred object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param mode within, overlap, or trim.
#' @param collapse Transcript collapse mode.
#' @param show_cds Whether to distinguish CDS and UTR segments.
#' @param cds_height Vertical thickness of CDS rectangles.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles.
#' @param highlight Optional data frame used to shade genomic or transcript intervals. It must contain `start` and `end` columns. Optional columns such as `label` or `group` are allowed but are not required. In `coordinate = "genomic"`, `start` and `end` are genomic positions; in `coordinate = "transcript"`, they are spliced transcript coordinates.
#' @param gene_palette RColorBrewer palette name used for discrete fills. If the number of discrete groups exceeds the palette maximum, colors are automatically interpolated.
#' @param gene_colors Optional custom fill colors for gene model features. Use a named vector such as `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")` to map colors explicitly. Unnamed colors are matched to the fixed gene-model levels `UTR`, `CDS`, and `exon`, so colors remain stable even when only one feature type is present.
#' @param gene_border_color Optional rectangle border color. Use `NA` to hide borders, or a color such as `"black"` or `"grey30"` to draw feature outlines.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw the panel border. `NULL` preserves the selected theme default.
#' @param label_position Where to draw feature labels. `axis` draws labels on the y axis and `feature` draws labels at the center of each gene/transcript structure.
#' @param label_by Which identifier to use for feature labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_size Text size in points for axis text, axis titles, and legends.
#' @param direction_mode Direction-arrow style. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @details
#' Common feature names used by `gene_colors` are `UTR`, `CDS`, and `exon`.
#' `highlight` must contain `start` and `end` columns. In genomic coordinate
#' plots these are genomic positions; in transcript coordinate plots they are
#' spliced transcript positions. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot object.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' plot_region(gp, chrom = "chr1", start = 12339001, end = 12352000)
#' @export
plot_region <- function(object, chrom, start, end, mode = c("within", "overlap", "trim"), collapse = c("none", "union_exon", "longest"), show_cds = TRUE, cds_height = 0.50, utr_height = 0.25, direction_mode = c("transcript", "gene", "end", "none"), highlight = NULL, gene_palette = "Paired", gene_colors = NULL, gene_border_color = NA, label_position = c("axis", "feature"), label_by = c("gene", "transcript"), text_size = 14, plot_theme = c("bw", "classic", "light", "minimal"), show_panel_border = NULL) {
  stop_if_not(is_gene_model_feature(object), "`object` must be a GenePred object or a Feature object with transcript/exon records.")
  object <- as_genepred(object)
  mode <- match.arg(mode)
  collapse <- match.arg(collapse)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  direction_mode <- match.arg(direction_mode)
  chrom_value <- as.character(chrom)
  start_value <- as.integer(start)
  end_value <- as.integer(end)

  gp <- retrieve_feature(
    object = object,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
    mode = mode,
    as = "Feature"
  )
  tx <- gp$transcripts
  ex <- gp$exons
  stop_if_not(nrow(tx) > 0L, "No transcripts were found in the specified region.")

  collapsed <- collapse_transcripts(tx, ex, collapse)
  plot_gene_model_core(
    tx = collapsed$transcripts,
    exons = collapsed$exons,
    coordinate = "genomic",
    collapse = collapse,
    show_cds = show_cds,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    highlight = highlight,
    gene_palette = gene_palette,
    gene_colors = gene_colors,
    gene_border_color = gene_border_color,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  ) +
    ggplot2::coord_cartesian(xlim = c(start_value, end_value))
}

plot_gene_model_core <- function(tx, exons, coordinate = c("genomic", "transcript"), collapse = "none", show_cds = TRUE, cds_height = 0.50, utr_height = 0.25, direction_mode = c("transcript", "gene", "end", "none"), highlight = NULL, gene_palette = "Paired", gene_colors = NULL, gene_border_color = NA, plot_theme = c("bw", "classic", "light", "minimal"), show_panel_border = NULL, label_position = c("axis", "feature"), label_by = c("gene", "transcript"), text_size = 14) {
  coordinate <- match.arg(coordinate)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  label_side <- "center"
  label_offset <- 0.45
  text_color <- "black"
  text_size <- as.numeric(text_size)[1L]
  if (!is.finite(text_size) || text_size <= 0) {
    stop("`text_size` must be a positive numeric value.", call. = FALSE)
  }
  label_size <- text_size * 0.85
  show_direction <- !identical(direction_mode, "none")
  label_offset <- as.numeric(label_offset)[1L]
  if (!is.finite(label_offset) || label_offset < 0) {
    stop("`label_offset` must be a non-negative numeric value.", call. = FALSE)
  }
  cds_height <- as.numeric(cds_height)[1L]
  utr_height <- as.numeric(utr_height)[1L]
  if (!is.finite(cds_height) || cds_height <= 0) {
    stop("`cds_height` must be a positive numeric value.", call. = FALSE)
  }
  if (!is.finite(utr_height) || utr_height <= 0) {
    stop("`utr_height` must be a positive numeric value.", call. = FALSE)
  }
  gene_border_color <- normalize_border_color(gene_border_color)

  tx <- data.table::copy(data.table::as.data.table(tx))
  ex <- data.table::copy(data.table::as.data.table(exons))

  if (nrow(tx) == 0L || nrow(ex) == 0L) {
    stop("No transcript or exon records are available for plotting.", call. = FALSE)
  }

  tx[, ".tx_span" := pmax(1L, as.integer(tx[["tx_end"]]) - as.integer(tx[["tx_start"]]) + 1L)]
  data.table::setorderv(tx, c("gene_id", ".tx_span", "tx_start", "tx_end", "transcript_id"))
  tx[, "y" := rev(seq_len(.N))]
  tx[, "plot_label" := if (label_by == "gene") as.character(tx[["gene_id"]]) else as.character(tx[["transcript_id"]])]

  y_map <- tx[, c("transcript_id", "y"), with = FALSE]
  ex <- merge(ex, y_map, by = "transcript_id", all.x = TRUE)

  if (show_cds) {
    seg <- make_feature_segments(tx, ex)
    seg <- merge(seg, y_map, by = "transcript_id", all.x = TRUE)
  } else {
    seg <- data.table::data.table(
      transcript_id = ex[["transcript_id"]],
      gene_id = ex[["gene_id"]],
      chrom = ex[["chrom"]],
      strand = ex[["strand"]],
      exon_number = ex[["exon_number"]],
      feature = "exon",
      start = as.integer(ex[["exon_start"]]),
      end = as.integer(ex[["exon_end"]]),
      y = as.integer(ex[["y"]])
    )
  }

  if (coordinate == "transcript") {
    mapped <- map_to_transcript_coordinate(ex, seg)
    ex <- mapped$exons
    seg <- mapped$segments
    line_dt <- ex[
      , .(
        line_start = as.integer(min(.SD[["plot_start"]], na.rm = TRUE)),
        line_end = as.integer(max(.SD[["plot_end"]], na.rm = TRUE)),
        y = as.integer(.SD[["y"]][1L])
      ),
      by = "transcript_id",
      .SDcols = c("plot_start", "plot_end", "y")
    ]
    label_lookup <- tx[, c("transcript_id", "gene_id", "strand", "plot_label"), with = FALSE]
    line_dt <- merge(line_dt, label_lookup, by = "transcript_id", all.x = TRUE)
  } else {
    ex[, "plot_start" := as.integer(ex[["exon_start"]])]
    ex[, "plot_end" := as.integer(ex[["exon_end"]])]
    seg[, "plot_start" := as.integer(seg[["start"]])]
    seg[, "plot_end" := as.integer(seg[["end"]])]
    line_dt <- data.table::data.table(
      transcript_id = tx[["transcript_id"]],
      line_start = as.integer(tx[["tx_start"]]),
      line_end = as.integer(tx[["tx_end"]]),
      y = as.integer(tx[["y"]]),
      gene_id = as.character(tx[["gene_id"]]),
      strand = as.character(tx[["strand"]]),
      plot_label = as.character(tx[["plot_label"]])
    )
  }

  seg[, "feature" := normalize_gene_model_feature(seg[["feature"]])]
  gene_feature_levels <- names(make_gene_model_fill_colors(
    color_palette = gene_palette,
    gene_colors = gene_colors,
    features = seg[["feature"]]
  ))
  seg[, "feature" := factor(as.character(seg[["feature"]]), levels = gene_feature_levels)]

  seg[, "feature_width" := ifelse(as.character(seg[["feature"]]) == "CDS", cds_height, utr_height)]
  if (!show_cds) {
    seg[, "feature_width" := cds_height]
  }
  seg[, "ymin" := as.numeric(seg[["y"]]) - as.numeric(seg[["feature_width"]]) / 2]
  seg[, "ymax" := as.numeric(seg[["y"]]) + as.numeric(seg[["feature_width"]]) / 2]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = line_dt,
      ggplot2::aes(x = .data$line_start, xend = .data$line_end, y = .data$y, yend = .data$y),
      linewidth = 0.35
    ) +
    ggplot2::geom_rect(
      data = seg,
      ggplot2::aes(
        xmin = .data$plot_start,
        xmax = .data$plot_end,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$feature
      ),
      color = gene_border_color,
      linewidth = 0.2
    )

  y_breaks <- tx[["y"]]
  if (label_position == "axis") {
    if (label_by == "gene") {
      y_labels <- rep("", length(y_breaks))
      gene_ids_for_axis <- as.character(tx[["gene_id"]])
      first_idx <- !duplicated(gene_ids_for_axis)
      y_labels[first_idx] <- as.character(tx[["plot_label"]][first_idx])
    } else {
      y_labels <- as.character(tx[["plot_label"]])
    }
  } else {
    y_labels <- rep("", length(y_breaks))
  }
  chrom_values <- unique(as.character(tx[["chrom"]]))
  chrom_values <- chrom_values[!is.na(chrom_values) & nzchar(chrom_values)]
  x_label <- if (coordinate == "transcript") {
    "Transcript coordinate (bp)"
  } else if (length(chrom_values) == 1L) {
    paste0("Chromosome ", chrom_values[1L], " position (bp)")
  } else {
    "Genomic coordinate (bp)"
  }

  p <- apply_gene_model_fill_scale(
    p,
    color_palette = gene_palette,
    gene_colors = gene_colors,
    features = seg[["feature"]]
  ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = ggplot2::expansion(mult = c(0.15, 0.15), add = c(0.6, 0.6))
    ) +
    ggplot2::labs(
      x = x_label,
      y = NULL
    ) +
    make_track_theme(
      plot_theme = plot_theme,
      show_panel_border = show_panel_border
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.text.y = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.x = ggplot2::element_text(color = text_color, size = text_size),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = text_color, size = text_size)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  if (label_position == "feature") {
    p <- p + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }

  if (label_position == "feature") {
    if (label_by == "gene") {
      label_dt <- line_dt[
        , .(
          line_start = as.integer(min(.SD[["line_start"]], na.rm = TRUE)),
          line_end = as.integer(max(.SD[["line_end"]], na.rm = TRUE)),
          y = as.numeric(mean(.SD[["y"]], na.rm = TRUE)),
          plot_label = as.character(.SD[["plot_label"]][1L])
        ),
        by = "gene_id",
        .SDcols = c("line_start", "line_end", "y", "plot_label")
      ]
    } else {
      label_dt <- data.table::copy(line_dt)
    }
    label_dt[, "label_x" := (as.numeric(label_dt[["line_start"]]) + as.numeric(label_dt[["line_end"]])) / 2]
    if (label_side == "above") {
      label_dt[, "label_y" := as.numeric(label_dt[["y"]]) + label_offset]
      label_vjust <- 0
    } else if (label_side == "below") {
      label_dt[, "label_y" := as.numeric(label_dt[["y"]]) - label_offset]
      label_vjust <- 1
    } else {
      label_dt[, "label_y" := as.numeric(label_dt[["y"]])]
      label_vjust <- 0.5
    }
    p <- p + ggplot2::geom_text(
      data = label_dt,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$plot_label),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = label_vjust,
      color = text_color,
      size = label_size / 2.845276
    )
  }

  if (show_direction && coordinate == "genomic") {
    if (direction_mode == "transcript") {
      arrow_dt <- data.table::data.table(
        x = as.integer(tx[["tx_start"]]),
        xend = as.integer(tx[["tx_end"]]),
        y = as.numeric(tx[["y"]]) + 0.38,
        strand = as.character(tx[["strand"]])
      )
    } else {
      arrow_dt <- line_dt[
        , .(
          x = as.integer(min(.SD[["line_start"]], na.rm = TRUE)),
          xend = as.integer(max(.SD[["line_end"]], na.rm = TRUE)),
          y = as.numeric(mean(.SD[["y"]], na.rm = TRUE)) + 0.38,
          strand = as.character(.SD[["strand"]][1L])
        ),
        by = "gene_id",
        .SDcols = c("line_start", "line_end", "y", "strand")
      ]
      if (direction_mode == "end") {
        span <- pmax(1, abs(as.numeric(arrow_dt[["xend"]]) - as.numeric(arrow_dt[["x"]])))
        short_len <- pmax(20, pmin(span * 0.15, span))
        plus_idx <- which(arrow_dt[["strand"]] != "-")
        minus_idx <- which(arrow_dt[["strand"]] == "-")
        if (length(plus_idx) > 0L) {
          arrow_dt[plus_idx, "x" := as.integer(as.numeric(arrow_dt[["xend"]][plus_idx]) - short_len[plus_idx])]
        }
        if (length(minus_idx) > 0L) {
          arrow_dt[minus_idx, "xend" := as.integer(as.numeric(arrow_dt[["x"]][minus_idx]) + short_len[minus_idx])]
        }
      }
    }
    minus_idx <- which(arrow_dt[["strand"]] == "-")
    if (length(minus_idx) > 0L) {
      old_x <- arrow_dt[["x"]][minus_idx]
      arrow_dt[minus_idx, "x" := as.integer(arrow_dt[["xend"]][minus_idx])]
      arrow_dt[minus_idx, "xend" := as.integer(old_x)]
    }
    p <- p + ggplot2::geom_segment(
      data = arrow_dt,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      arrow = grid::arrow(length = grid::unit(0.08, "inches")),
      linewidth = 0.25
    )
  }

  add_highlight_layer(p, highlight)
}
