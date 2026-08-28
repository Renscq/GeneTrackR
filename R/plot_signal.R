# Author: Rensc
# Date: 2026-05-27
# Version: dev006
# Function: Plot signal tracks for transcripts, genes, and genomic regions
# Input: BwgTrack and optional GenePred objects
# Output: ggplot signal figures

#' Plot signal over a transcript
#'
#' @param signal A BwgTrack object.
#' @param annotation A GenePred object.
#' @param transcript_id Transcript ID.
#' @param samples Optional sample IDs to plot. If NULL, all samples are used.
#' @param sample_groups Optional sample group mapping for group-level coloring or replicate summaries. Use a named character vector, a data frame with `sample_id` and `group`, or an unnamed vector with one group per selected sample.
#' @param signal_color_by Color signal tracks by `sample` or `group`.
#' @param signal_summary Replicate summary mode. Use `none` to plot individual samples, or `mean`, `median`, or `sum` to summarize samples within each group. Summary is performed on the current signal intervals, so using `bin_size` is recommended when raw interval boundaries differ among samples.
#' @param coordinate Coordinate mode. `transcript` removes introns and displays spliced transcript coordinates, while `genomic` keeps genomic coordinates.
#' @param plot_type Signal plot type: `bar`, `line`, `heatmap`, or `frame`. `frame` is designed for Ribo-seq style transcript plots and colors CDS positions by frame0/frame1/frame2.
#' @param strand Strand selector. Use `auto` to use the transcript strand.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns. For `coordinate = "genomic"`, these are genomic coordinates; for `coordinate = "transcript"`, these are spliced transcript coordinates.
#' @param show_gene_model Whether to append the transcript gene model track. Default TRUE.
#' @param signal_track_height Relative height of the signal panel when the gene model is shown. Default 3.
#' @param gene_track_height Relative height of the gene model panel when it is shown. Default 1.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. Discrete sample/group colors are assigned in the standard RColorBrewer class order; heatmaps use the corresponding continuous gradient.
#' @param signal_palette_direction Direction for generated signal colors. Use `1` for the standard palette order and `-1` for the reversed palette order. Discrete sample/group colors preserve level-to-color order; heatmap gradients reverse continuously.
#' @param frame_palette RColorBrewer palette for CDS frame colors when `plot_type = "frame"`. Default is `Paired`. Colors are assigned in order to `frame0`, `frame1`, and `frame2`.
#' @param frame_colors Optional named colors for `frame0`, `frame1`, and `frame2`.
#' @param signal_colors Optional named or unnamed vector of explicit colors for samples. If supplied, it overrides `signal_palette`; explicit colors are not modified by `signal_palette_direction`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param signal_y_limits Optional two-element numeric vector giving the plotted y-axis limits after `signal_transform`. Supplying limits changes `signal_y_scale` to `fixed`.
#' @param signal_alpha Signal geometry transparency from 0 to 1.
#' @param signal_bar_width Relative width of bar intervals from greater than 0 to 1. A value below 1 creates proportional gaps without changing interval centers.
#' @param cds_height Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis and `feature` draws labels on the model.
#' @param label_by Which identifier to use for gene model labels.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw panel borders. `NULL` preserves the selected theme default.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @details
#' `samples` selects the samples to draw. `sample_groups` can be a named vector
#' or a data frame with `sample_id` and `group`. Use `signal_color_by = "group"`
#' to color by group, and `signal_summary = "mean"`, `"median"`, or `"sum"`
#' to collapse replicates within groups. `signal_colors` overrides
#' `signal_palette` and may be named by sample ID or group name. See also
#' [GeneTrackR-advanced-parameters].
#' @return A ggplot or patchwork object.
#' @examples
#' \dontrun{
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt",
#'   verbose = FALSE
#' )
#' riboseq <- read_bwg(
#'   system.file(
#'     "extdata",
#'     c("gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"),
#'     package = "GeneTrackR"
#'   ),
#'   format = "bedgraph",
#'   sample_names = c("Ribo_seq_plus", "Ribo_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#' plot_signal_transcript(
#'   signal = riboseq,
#'   annotation = gp,
#'   transcript_id = "TxA1",
#'   coordinate = "transcript",
#'   plot_type = "frame"
#' )
#' }
#' @export
plot_signal_transcript <- function(
  signal,
  annotation,
  transcript_id,
  samples = NULL,
  sample_groups = NULL,
  signal_color_by = c("sample", "group"),
  signal_summary = c("none", "mean", "median", "sum"),
  coordinate = c("transcript", "genomic"),
  plot_type = c("bar", "line", "heatmap", "frame"),
  strand = c("auto", "+", "-", "both", "ignore"),
  bin_size = NULL,
  highlight = NULL,
  show_gene_model = TRUE,
  signal_track_height = 3,
  gene_track_height = 1,
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL,
  frame_palette = "Paired",
  frame_colors = NULL,
  signal_transform = c("none", "log2", "log10", "sqrt"),
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  heatmap_bin_size = NULL,
  heatmap_max_bins = 800L,
  heatmap_summary = c("mean", "max", "sum", "median"),
  cds_height = 0.50,
  utr_height = 0.25,
  direction_mode = c("transcript", "gene", "end", "none"),
  label_position = c("axis", "feature"),
  label_by = c("gene", "transcript"),
  text_size = 14,
  signal_y_limits = NULL,
  signal_alpha = 0.85,
  signal_bar_width = 1,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
) {
  stop_if_not(
    inherits(signal, "BwgTrack"),
    "`signal` must be a BwgTrack object."
  )
  stop_if_not(
    inherits(annotation, "GenePred"),
    "`annotation` must be a GenePred object."
  )
  coordinate <- match.arg(coordinate)
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  signal_track_height <- normalize_track_height(
    signal_track_height,
    "signal_track_height",
    3
  )
  gene_track_height <- normalize_track_height(
    gene_track_height,
    "gene_track_height",
    1
  )
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  text_color <- "black"
  grid_linewidth <- NULL

  transcript_id_value <- as.character(transcript_id)
  tx_all <- annotation$transcripts
  ex_all <- annotation$exons
  tx <- tx_all[tx_all[["transcript_id"]] == transcript_id_value]
  stop_if_not(nrow(tx) > 0L, "Transcript ID was not found.")
  ex <- ex_all[ex_all[["transcript_id"]] == transcript_id_value]
  selected_strand <- if (strand == "auto") tx$strand[1] else strand
  expected_samples <- get_expected_signal_samples(
    signal,
    samples = samples,
    strand = selected_strand
  )

  dt <- retrieve_bwg(
    signal,
    tx$chrom[1],
    tx$tx_start[1],
    tx$tx_end[1],
    samples = samples,
    strand = selected_strand
  )
  if (coordinate == "transcript") {
    query_start <- 1L
    query_end <- as.integer(sum(
      as.integer(ex[["exon_end"]]) - as.integer(ex[["exon_start"]]) + 1L,
      na.rm = TRUE
    ))
    if (plot_type != "frame") {
      dt <- map_signal_to_exons(dt, ex)
    }
  } else {
    query_start <- as.integer(tx$tx_start[1])
    query_end <- as.integer(tx$tx_end[1])
  }
  if (plot_type == "frame") {
    if (!is.null(bin_size) && !identical(as.integer(bin_size)[1L], 1L)) {
      message(
        "[GeneTrackR] `plot_type = 'frame'` requires base-resolution signal; ",
        "`bin_size` was ignored."
      )
    }
    frame_anno <- build_transcript_frame_annotation(
      annotation = annotation,
      transcript_id = transcript_id_value,
      coordinate = coordinate
    )
    p_signal <- plot_signal_frame_core(
      dt = dt,
      frame_annotation = frame_anno,
      expected_samples = expected_samples,
      sample_groups = sample_groups,
      signal_summary = signal_summary,
      signal_transform = signal_transform,
      signal_y_scale = signal_y_scale,
      signal_y_ticks = signal_y_ticks,
      signal_y_limits = signal_y_limits,
      signal_alpha = signal_alpha,
      signal_bar_width = signal_bar_width,
      frame_palette = frame_palette,
      frame_colors = frame_colors,
      x_label = ifelse(
        coordinate == "transcript",
        "Transcript coordinate",
        paste0("Chromosome ", as.character(tx$chrom[1]), " position (bp)")
      ),
      text_size = text_size,
      grid_linewidth = grid_linewidth,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      highlight = highlight
    )
  } else {
    if (!is.null(bin_size)) {
      dt <- bin_bwg(dt, bin_size = bin_size)
    }
    dt <- complete_empty_signal_tracks(
      dt,
      sample_ids = expected_samples,
      chrom = tx$chrom[1],
      start = query_start,
      end = query_end,
      strand = selected_strand
    )

    p_signal <- plot_signal_core(
      dt,
      plot_type = plot_type,
      highlight = highlight,
      x_label = ifelse(
        coordinate == "transcript",
        "Transcript coordinate",
        paste0("Chromosome ", as.character(tx$chrom[1]), " position (bp)")
      ),
      signal_palette = signal_palette,
      signal_palette_direction = signal_palette_direction,
      signal_colors = signal_colors,
      sample_groups = sample_groups,
      signal_color_by = signal_color_by,
      signal_summary = signal_summary,
      signal_transform = signal_transform,
      signal_y_scale = signal_y_scale,
      signal_y_ticks = signal_y_ticks,
      signal_y_limits = signal_y_limits,
      signal_alpha = signal_alpha,
      signal_bar_width = signal_bar_width,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      text_size = text_size,
      heatmap_bin_size = heatmap_bin_size,
      heatmap_max_bins = heatmap_max_bins,
      heatmap_summary = heatmap_summary
    )
  }

  # The gene model track already uses the correct transcript-level coordinate
  # range and its own ggplot expansion. Do not modify it here. Instead, only
  # force the signal panel to use the full transcript range, because signal data
  # often cover only a subset of the transcript and ggplot would otherwise shrink
  # the x-axis to the covered intervals. Keeping the default x-scale expansion
  # makes the signal panel align with the unmodified gene model panel.
  shared_x_limits <- c(as.numeric(query_start), as.numeric(query_end))
  p_signal <- p_signal + ggplot2::scale_x_continuous(limits = shared_x_limits)

  if (!show_gene_model) {
    return(p_signal)
  }
  p_model <- plot_transcript(
    annotation,
    transcript_id = transcript_id,
    coordinate = coordinate,
    show_cds = TRUE,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    highlight = highlight,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  )
  patchwork::wrap_plots(
    p_signal,
    p_model,
    ncol = 1,
    heights = c(signal_track_height, gene_track_height)
  )
}

#' Plot signal over a gene
#'
#' @param signal A BwgTrack object.
#' @param annotation A GenePred object.
#' @param gene_id Gene ID.
#' @param samples Optional sample IDs to plot. If NULL, all samples are used.
#' @param sample_groups Optional sample group mapping for group-level coloring or replicate summaries. Use a named character vector, a data frame with `sample_id` and `group`, or an unnamed vector with one group per selected sample.
#' @param signal_color_by Color signal tracks by `sample` or `group`.
#' @param signal_summary Replicate summary mode. Use `none` to plot individual samples, or `mean`, `median`, or `sum` to summarize samples within each group. Summary is performed on the current signal intervals, so using `bin_size` is recommended when raw interval boundaries differ among samples.
#' @param plot_type Signal plot type: `bar`, `line`, or `heatmap`.
#' @param strand Strand selector. Use `auto` to use the gene strand.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates.
#' @param show_gene_model Whether to append the gene model track.
#' @param signal_track_height Relative height of the signal panel when the gene model is shown. Default 3.
#' @param gene_track_height Relative height of the gene model panel when it is shown. Default 1.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. Discrete sample/group colors are assigned in the standard RColorBrewer class order; heatmaps use the corresponding continuous gradient.
#' @param signal_palette_direction Direction for generated signal colors. Use `1` for the standard palette order and `-1` for the reversed palette order. Discrete sample/group colors preserve level-to-color order; heatmap gradients reverse continuously.
#' @param signal_colors Optional named or unnamed vector of explicit colors for samples. If supplied, it overrides `signal_palette`; explicit colors are not modified by `signal_palette_direction`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param signal_y_limits Optional two-element numeric vector giving the plotted y-axis limits after `signal_transform`. Supplying limits changes `signal_y_scale` to `fixed`.
#' @param signal_alpha Signal geometry transparency from 0 to 1.
#' @param signal_bar_width Relative width of bar intervals from greater than 0 to 1. A value below 1 creates proportional gaps without changing interval centers.
#' @param cds_height Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis and `feature` draws labels at the center of each gene/transcript structure.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw panel borders. `NULL` preserves the selected theme default.
#' @details
#' `samples` selects the samples to draw. `sample_groups` can be used for
#' group-level coloring and replicate summaries. When raw intervals differ among
#' samples, set `bin_size` before using `signal_summary` so the summary is made
#' on comparable bins. `signal_transform` changes the plotted y value only; the
#' queried signal table is not modified. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot or patchwork object.
#' @examples
#' \dontrun{
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt",
#'   verbose = FALSE
#' )
#' rnaseq <- read_bwg(
#'   system.file(
#'     "extdata",
#'     c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
#'     package = "GeneTrackR"
#'   ),
#'   format = "bedgraph",
#'   sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
#'   strand = c("+", "-"),
#'   mode = "memory"
#' )
#' plot_signal_gene(
#'   signal = rnaseq,
#'   annotation = gp,
#'   gene_id = "GeneA",
#'   plot_type = "bar",
#'   signal_palette = "Paired",
#'   signal_palette_direction = -1,
#'   signal_y_scale = "fixed",
#'   signal_y_ticks = "pretty",
#'   signal_track_height = 4,
#'   gene_track_height = 1
#' )
#' }
#' @export
plot_signal_gene <- function(
  signal,
  annotation,
  gene_id,
  samples = NULL,
  sample_groups = NULL,
  signal_color_by = c("sample", "group"),
  signal_summary = c("none", "mean", "median", "sum"),
  plot_type = c("bar", "line", "heatmap"),
  strand = c("auto", "+", "-", "both", "ignore"),
  bin_size = NULL,
  highlight = NULL,
  show_gene_model = TRUE,
  signal_track_height = 3,
  gene_track_height = 1,
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL,
  signal_transform = c("none", "log2", "log10", "sqrt"),
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  heatmap_bin_size = NULL,
  heatmap_max_bins = 800L,
  heatmap_summary = c("mean", "max", "sum", "median"),
  cds_height = 0.50,
  utr_height = 0.25,
  direction_mode = c("transcript", "gene", "end", "none"),
  label_position = c("axis", "feature"),
  label_by = c("gene", "transcript"),
  text_size = 14,
  signal_y_limits = NULL,
  signal_alpha = 0.85,
  signal_bar_width = 1,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
) {
  stop_if_not(
    inherits(signal, "BwgTrack"),
    "`signal` must be a BwgTrack object."
  )
  stop_if_not(
    inherits(annotation, "GenePred"),
    "`annotation` must be a GenePred object."
  )
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  signal_track_height <- normalize_track_height(
    signal_track_height,
    "signal_track_height",
    3
  )
  gene_track_height <- normalize_track_height(
    gene_track_height,
    "gene_track_height",
    1
  )
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  text_color <- "black"
  grid_linewidth <- NULL

  gene_id_value <- as.character(gene_id)
  tx_all <- annotation$transcripts
  tx <- tx_all[tx_all[["gene_id"]] == gene_id_value]
  stop_if_not(nrow(tx) > 0L, "Gene ID was not found.")
  gene <- build_gene_table(tx)
  selected_strand <- if (strand == "auto") gene$strand[1] else strand
  expected_samples <- get_expected_signal_samples(
    signal,
    samples = samples,
    strand = selected_strand
  )

  dt <- retrieve_bwg(
    signal,
    gene$chrom[1],
    gene$gene_start[1],
    gene$gene_end[1],
    samples = samples,
    strand = selected_strand
  )
  if (!is.null(bin_size)) {
    dt <- bin_bwg(dt, bin_size = bin_size)
  }
  dt <- complete_empty_signal_tracks(
    dt,
    sample_ids = expected_samples,
    chrom = gene$chrom[1],
    start = gene$gene_start[1],
    end = gene$gene_end[1],
    strand = selected_strand
  )
  p_signal <- plot_signal_core(
    dt,
    plot_type = plot_type,
    highlight = highlight,
    x_label = paste0(
      "Chromosome ",
      as.character(gene$chrom[1]),
      " position (bp)"
    ),
    signal_palette = signal_palette,
    signal_palette_direction = signal_palette_direction,
    signal_colors = signal_colors,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    signal_y_limits = signal_y_limits,
    signal_alpha = signal_alpha,
    signal_bar_width = signal_bar_width,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    text_size = text_size,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary
  )
  # Keep the gene model track unchanged and use it as the coordinate reference.
  # Only expand/constrain the signal panel to the complete gene range; otherwise
  # ggplot may shrink the x-axis to the covered signal intervals.
  p_signal <- p_signal +
    ggplot2::scale_x_continuous(
      limits = c(
        as.numeric(gene[["gene_start"]][1L]),
        as.numeric(gene[["gene_end"]][1L])
      )
    )

  if (!show_gene_model) {
    return(p_signal)
  }
  p_model <- plot_gene(
    annotation,
    gene_id = gene_id,
    collapse = "none",
    coordinate = "genomic",
    highlight = highlight,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  )
  patchwork::wrap_plots(
    p_signal,
    p_model,
    ncol = 1,
    heights = c(signal_track_height, gene_track_height)
  )
}

#' Plot signal over a genomic region
#'
#' @param signal A BwgTrack object.
#' @param chrom Chromosome name.
#' @param start Region start.
#' @param end Region end.
#' @param samples Optional sample IDs to plot. If NULL, all samples are used.
#' @param sample_groups Optional sample group mapping for group-level coloring or replicate summaries. Use a named character vector, a data frame with `sample_id` and `group`, or an unnamed vector with one group per selected sample.
#' @param signal_color_by Color signal tracks by `sample` or `group`.
#' @param signal_summary Replicate summary mode. Use `none` to plot individual samples, or `mean`, `median`, or `sum` to summarize samples within each group. Summary is performed on the current signal intervals, so using `bin_size` is recommended when raw interval boundaries differ among samples.
#' @param plot_type Signal plot type: `bar`, `line`, or `heatmap`.
#' @param strand Strand selector.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates.
#' @param annotation Optional GenePred object.
#' @param show_gene_model Whether to append a gene model track. Default TRUE.
#' @param signal_track_height Relative height of the signal panel when the gene model is shown. Default 3.
#' @param gene_track_height Relative height of the gene model panel when it is shown. Default 1.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. Discrete sample/group colors are assigned in the standard RColorBrewer class order; heatmaps use the corresponding continuous gradient.
#' @param signal_palette_direction Direction for generated signal colors. Use `1` for the standard palette order and `-1` for the reversed palette order. Discrete sample/group colors preserve level-to-color order; heatmap gradients reverse continuously.
#' @param signal_colors Optional named or unnamed vector of explicit colors for samples. If supplied, it overrides `signal_palette`; explicit colors are not modified by `signal_palette_direction`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param signal_y_limits Optional two-element numeric vector giving the plotted y-axis limits after `signal_transform`. Supplying limits changes `signal_y_scale` to `fixed`.
#' @param signal_alpha Signal geometry transparency from 0 to 1.
#' @param signal_bar_width Relative width of bar intervals from greater than 0 to 1. A value below 1 creates proportional gaps without changing interval centers.
#' @param cds_height Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis and `feature` draws labels at the center of each gene/transcript structure.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @param plot_theme Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw panel borders. `NULL` preserves the selected theme default.
#' @details
#' If `annotation` is supplied and `show_gene_model = TRUE`, a gene model track
#' is appended below the signal panel. `signal_colors` can be an unnamed vector
#' or a named vector. Named values are matched to sample IDs when
#' `signal_color_by = "sample"`, and to group names when
#' `signal_color_by = "group"`. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot or patchwork object.
#' @examples
#' \dontrun{
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt",
#'   verbose = FALSE
#' )
#' signal_all <- read_bwg(
#'   system.file(
#'     "extdata",
#'     c(
#'       "gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph",
#'       "gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"
#'     ),
#'     package = "GeneTrackR"
#'   ),
#'   format = "bedgraph",
#'   sample_names = c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus"),
#'   strand = c("+", "-", "+", "-"),
#'   mode = "memory"
#' )
#' plot_signal_region(
#'   signal = signal_all,
#'   annotation = gp,
#'   chrom = "chr1",
#'   start = 12339001,
#'   end = 12374500,
#'   strand = "both",
#'   plot_type = "bar"
#' )
#' }
#' @export
plot_signal_region <- function(
  signal,
  chrom,
  start,
  end,
  samples = NULL,
  sample_groups = NULL,
  signal_color_by = c("sample", "group"),
  signal_summary = c("none", "mean", "median", "sum"),
  plot_type = c("bar", "line", "heatmap"),
  strand = c("ignore", "+", "-", "both"),
  bin_size = NULL,
  highlight = NULL,
  annotation = NULL,
  show_gene_model = TRUE,
  signal_track_height = 3,
  gene_track_height = 1,
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL,
  signal_transform = c("none", "log2", "log10", "sqrt"),
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  heatmap_bin_size = NULL,
  heatmap_max_bins = 800L,
  heatmap_summary = c("mean", "max", "sum", "median"),
  cds_height = 0.50,
  utr_height = 0.25,
  direction_mode = c("transcript", "gene", "end", "none"),
  label_position = c("axis", "feature"),
  label_by = c("gene", "transcript"),
  text_size = 14,
  signal_y_limits = NULL,
  signal_alpha = 0.85,
  signal_bar_width = 1,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
) {
  stop_if_not(
    inherits(signal, "BwgTrack"),
    "`signal` must be a BwgTrack object."
  )
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  signal_track_height <- normalize_track_height(
    signal_track_height,
    "signal_track_height",
    3
  )
  gene_track_height <- normalize_track_height(
    gene_track_height,
    "gene_track_height",
    1
  )
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  text_color <- "black"
  grid_linewidth <- NULL

  expected_samples <- get_expected_signal_samples(
    signal,
    samples = samples,
    strand = strand
  )
  dt <- retrieve_bwg(
    signal,
    chrom,
    start,
    end,
    samples = samples,
    strand = strand
  )
  if (!is.null(bin_size)) {
    dt <- bin_bwg(dt, bin_size = bin_size)
  }
  dt <- complete_empty_signal_tracks(
    dt,
    sample_ids = expected_samples,
    chrom = chrom,
    start = start,
    end = end,
    strand = strand
  )
  p_signal <- plot_signal_core(
    dt,
    plot_type = plot_type,
    highlight = highlight,
    x_label = paste0("Chromosome ", as.character(chrom), " position (bp)"),
    signal_palette = signal_palette,
    signal_palette_direction = signal_palette_direction,
    signal_colors = signal_colors,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    signal_y_limits = signal_y_limits,
    signal_alpha = signal_alpha,
    signal_bar_width = signal_bar_width,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    text_size = text_size,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary
  ) +
    ggplot2::coord_cartesian(xlim = c(start, end))

  if (is.null(annotation) || !show_gene_model) {
    return(p_signal)
  }
  p_model <- plot_region(
    annotation,
    chrom,
    start,
    end,
    mode = "overlap",
    collapse = "none",
    highlight = highlight,
    cds_height = cds_height,
    utr_height = utr_height,
    direction_mode = direction_mode,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    label_position = label_position,
    label_by = label_by,
    text_size = text_size
  )
  patchwork::wrap_plots(
    p_signal,
    p_model,
    ncol = 1,
    heights = c(signal_track_height, gene_track_height)
  )
}


get_expected_signal_samples <- function(
  signal,
  samples = NULL,
  strand = "ignore",
  strand_policy = "ignore_unstranded"
) {
  sample_tbl <- data.table::copy(data.table::as.data.table(signal$samples))
  if (!"has_strand" %in% names(sample_tbl)) {
    sample_tbl[, "has_strand" := FALSE]
  }
  if (!"strand" %in% names(sample_tbl)) {
    sample_tbl[, "strand" := "*"]
  }
  if (!is.null(samples)) {
    samples_value <- as.character(samples)
    sample_tbl <- sample_tbl[sample_tbl[["sample_id"]] %in% samples_value]
  }
  sample_tbl <- filter_sample_table_by_strand(
    sample_tbl = sample_tbl,
    strand = strand,
    strand_policy = strand_policy
  )
  unique(as.character(sample_tbl[["sample_id"]]))
}

complete_empty_signal_tracks <- function(
  dt,
  sample_ids,
  chrom,
  start,
  end,
  strand = "*"
) {
  sample_ids <- unique(as.character(sample_ids))
  if (length(sample_ids) == 0L) {
    return(data.table::as.data.table(dt))
  }

  chrom_value <- as.character(chrom)[1L]
  start_value <- as.integer(start)[1L]
  end_value <- as.integer(end)[1L]
  strand_value <- as.character(strand)[1L]
  if (is.na(strand_value) || strand_value %in% c("ignore", "both", "auto")) {
    strand_value <- "*"
  }

  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L) {
    present <- character()
    base_cols <- c("sample_id", "chrom", "start", "end", "value", "strand")
    dt <- data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    )
  } else {
    dt <- ensure_signal_strand_column(dt)
    present <- unique(as.character(dt[["sample_id"]]))
  }

  missing_samples <- setdiff(sample_ids, present)
  if (length(missing_samples) == 0L) {
    dt[,
      "sample_id" := factor(
        as.character(dt[["sample_id"]]),
        levels = sample_ids
      )
    ]
    return(dt[])
  }

  empty_dt <- data.table::data.table(
    sample_id = missing_samples,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
    value = 0,
    strand = strand_value
  )

  out <- data.table::rbindlist(list(dt, empty_dt), fill = TRUE)
  out[,
    "sample_id" := factor(as.character(out[["sample_id"]]), levels = sample_ids)
  ]
  out[]
}


aggregate_signal_for_heatmap <- function(
  dt,
  heatmap_bin_size = NULL,
  heatmap_max_bins = 800L,
  heatmap_summary = c("mean", "max", "sum", "median")
) {
  heatmap_summary <- match.arg(heatmap_summary)
  text_color <- "black"
  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L) {
    return(dt[])
  }

  x_min <- suppressWarnings(min(as.integer(dt[["start"]]), na.rm = TRUE))
  x_max <- suppressWarnings(max(as.integer(dt[["end"]]), na.rm = TRUE))
  if (!is.finite(x_min) || !is.finite(x_max) || x_max < x_min) {
    return(dt[])
  }

  region_width <- as.integer(x_max - x_min + 1L)
  if (is.null(heatmap_bin_size)) {
    heatmap_max_bins <- suppressWarnings(as.integer(heatmap_max_bins)[1L])
    if (!is.finite(heatmap_max_bins) || heatmap_max_bins <= 0L) {
      heatmap_max_bins <- 800L
    }
    interval_count <- length(unique(paste(
      dt[["sample_id"]],
      dt[["start"]],
      dt[["end"]],
      sep = "\r"
    )))
    if (
      region_width <= heatmap_max_bins &&
        interval_count <= heatmap_max_bins * length(unique(dt[["sample_id"]]))
    ) {
      return(dt[])
    }
    heatmap_bin_size <- max(
      1L,
      as.integer(ceiling(region_width / heatmap_max_bins))
    )
  } else {
    heatmap_bin_size <- suppressWarnings(as.integer(heatmap_bin_size)[1L])
    if (!is.finite(heatmap_bin_size) || heatmap_bin_size <= 0L) {
      return(dt[])
    }
  }

  dt[, heatmap_bin := as.integer(floor((mid - x_min) / heatmap_bin_size))]
  dt[heatmap_bin < 0L, heatmap_bin := 0L]

  group_cols <- c("sample_id", "heatmap_bin")
  if ("sample_group" %in% names(dt)) {
    group_cols <- c("sample_group", group_cols)
  }
  if ("chrom" %in% names(dt)) {
    group_cols <- c("chrom", group_cols)
  }
  if ("strand" %in% names(dt)) {
    group_cols <- c("strand", group_cols)
  }

  value_fun <- switch(
    heatmap_summary,
    mean = function(x) mean(x, na.rm = TRUE),
    max = function(x) max(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE)
  )

  out <- dt[,
    .(
      start = x_min + min(heatmap_bin) * heatmap_bin_size,
      end = pmin(
        x_max,
        x_min + (max(heatmap_bin) + 1L) * heatmap_bin_size - 1L
      ),
      value = value_fun(value),
      plot_value = value_fun(plot_value)
    ),
    by = group_cols
  ]
  out[, mid := (start + end) / 2]
  out[, heatmap_bin := NULL]
  out[]
}


expand_bwg_to_positions <- function(dt) {
  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L) {
    out <- data.table::copy(dt)
    out[, "genomic_pos" := integer()]
    return(out)
  }
  dt <- dt[!is.na(start) & !is.na(end) & as.integer(start) <= as.integer(end)]
  if (nrow(dt) == 0L) {
    out <- data.table::copy(dt)
    out[, "genomic_pos" := integer()]
    return(out)
  }
  widths <- as.integer(dt[["end"]]) - as.integer(dt[["start"]]) + 1L
  idx <- rep(seq_len(nrow(dt)), widths)
  pos <- unlist(
    Map(seq.int, as.integer(dt[["start"]]), as.integer(dt[["end"]])),
    use.names = FALSE
  )
  out <- dt[idx]
  out[, "genomic_pos" := as.integer(pos)]
  out[, `:=`(start = as.integer(genomic_pos), end = as.integer(genomic_pos))]
  out[]
}

build_transcript_frame_annotation <- function(
  annotation,
  transcript_id,
  coordinate = c("transcript", "genomic")
) {
  coordinate <- match.arg(coordinate)
  stop_if_not(
    inherits(annotation, "GenePred"),
    "`annotation` must be a GenePred object."
  )
  query_transcript_id <- as.character(transcript_id)[1L]
  stop_if_not(
    !is.na(query_transcript_id) && nzchar(query_transcript_id),
    "`transcript_id` must be a non-empty transcript ID."
  )

  tx <- data.table::as.data.table(annotation$transcripts)
  ex <- data.table::as.data.table(annotation$exons)

  # Avoid data.table non-standard evaluation collisions between the function
  # argument `transcript_id` and the column named `transcript_id`.
  tx <- tx[as.character(tx[["transcript_id"]]) == query_transcript_id]
  ex <- ex[as.character(ex[["transcript_id"]]) == query_transcript_id]

  if (nrow(tx) != 1L) {
    matched_ids <- unique(as.character(tx[["transcript_id"]]))
    msg <- if (nrow(tx) == 0L) {
      paste0("Transcript ID was not found: ", query_transcript_id)
    } else {
      paste0(
        "`transcript_id` must match exactly one transcript. Matched ",
        nrow(tx),
        " records for: ",
        query_transcript_id,
        ". Please check whether the annotation contains duplicated transcript IDs."
      )
    }
    stop(msg, call. = FALSE)
  }
  stop_if_not(
    nrow(ex) > 0L,
    "No exon records were found for the selected transcript."
  )

  strand_value <- as.character(tx[["strand"]][1L])
  cds_start <- suppressWarnings(as.integer(tx[["cds_start"]][1L]))
  cds_end <- suppressWarnings(as.integer(tx[["cds_end"]][1L]))
  has_cds <- !is.na(cds_start) && !is.na(cds_end) && cds_start <= cds_end

  data.table::setorder(ex, exon_start, exon_end)
  if (identical(strand_value, "-")) {
    ex_order <- ex[order(-as.integer(exon_end), -as.integer(exon_start))]
  } else {
    ex_order <- ex[order(as.integer(exon_start), as.integer(exon_end))]
  }

  # Build biological transcript order for CDS frame assignment. For negative
  # strand transcripts, the CDS frame must be counted from high genomic
  # coordinates to low genomic coordinates.
  pieces <- vector("list", nrow(ex_order))
  tx_offset <- 0L
  for (i in seq_len(nrow(ex_order))) {
    exon_start <- as.integer(ex_order[["exon_start"]][i])
    exon_end <- as.integer(ex_order[["exon_end"]][i])
    if (identical(strand_value, "-")) {
      gpos <- seq.int(exon_end, exon_start)
    } else {
      gpos <- seq.int(exon_start, exon_end)
    }
    n_pos <- length(gpos)
    pieces[[i]] <- data.table::data.table(
      chrom = as.character(ex_order[["chrom"]][i]),
      genomic_pos = as.integer(gpos),
      biological_transcript_pos = as.integer(tx_offset + seq_len(n_pos))
    )
    tx_offset <- tx_offset + n_pos
  }
  anno <- data.table::rbindlist(pieces, use.names = TRUE)
  anno[, "region" := "UTR"]
  anno[, "frame" := NA_character_]
  anno[, "cds_pos" := NA_integer_]

  if (isTRUE(has_cds)) {
    cds_idx <- which(
      anno[["genomic_pos"]] >= cds_start & anno[["genomic_pos"]] <= cds_end
    )
    if (length(cds_idx) > 0L) {
      anno[cds_idx, "region" := "CDS"]
      anno[cds_idx, "cds_pos" := seq_along(cds_idx)]
      anno[cds_idx, "frame" := paste0("frame", (seq_along(cds_idx) - 1L) %% 3L)]
    }
  }

  # Build plotting transcript coordinates with the same convention as
  # plot_transcript(coordinate = "transcript"). The existing gene model track
  # maps exons from low genomic coordinates to high genomic coordinates for the
  # x-axis, regardless of strand. Keeping this separate from the biological CDS
  # frame direction prevents negative-strand signal tracks from being reversed
  # relative to the gene model track.
  ex_plot <- data.table::copy(ex)
  data.table::setorderv(ex_plot, c("exon_start", "exon_end"))
  ex_plot[, "exon_width" := as.integer(exon_end - exon_start + 1L)]
  ex_plot[,
    "plot_offset" := as.integer(cumsum(data.table::shift(
      exon_width,
      fill = 0L
    )))
  ]

  plot_pieces <- vector("list", nrow(ex_plot))
  for (i in seq_len(nrow(ex_plot))) {
    exon_start <- as.integer(ex_plot[["exon_start"]][i])
    exon_end <- as.integer(ex_plot[["exon_end"]][i])
    gpos <- seq.int(exon_start, exon_end)
    plot_pieces[[i]] <- data.table::data.table(
      genomic_pos = as.integer(gpos),
      transcript_pos = as.integer(ex_plot[["plot_offset"]][i] + seq_along(gpos))
    )
  }
  plot_map <- data.table::rbindlist(plot_pieces, use.names = TRUE)
  anno <- merge(anno, plot_map, by = "genomic_pos", all.x = TRUE, sort = FALSE)

  if (coordinate == "transcript") {
    anno[, `:=`(
      plot_pos = as.integer(transcript_pos),
      start = as.integer(transcript_pos),
      end = as.integer(transcript_pos)
    )]
  } else {
    anno[, `:=`(
      plot_pos = as.integer(genomic_pos),
      start = as.integer(genomic_pos),
      end = as.integer(genomic_pos)
    )]
  }
  anno[]
}

make_frame_colors <- function(frame_palette = "Paired", frame_colors = NULL) {
  frame_levels <- c("frame0", "frame1", "frame2")
  if (!is.null(frame_colors)) {
    frame_colors <- as.character(frame_colors)
    if (
      is.null(names(frame_colors)) ||
        !all(frame_levels %in% names(frame_colors))
    ) {
      if (length(frame_colors) < 3L) {
        frame_colors <- grDevices::colorRampPalette(frame_colors)(3L)
      }
      frame_colors <- frame_colors[seq_len(3L)]
      names(frame_colors) <- frame_levels
    }
    return(frame_colors[frame_levels])
  }
  cols <- normalize_signal_colors(
    sample_ids = frame_levels,
    signal_palette = frame_palette,
    signal_colors = NULL
  )
  cols[frame_levels]
}

plot_signal_frame_core <- function(
  dt,
  frame_annotation,
  expected_samples,
  sample_groups = NULL,
  signal_summary = c("none", "mean", "median", "sum"),
  signal_transform = c("none", "log2", "log10", "sqrt"),
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  signal_y_limits = NULL,
  signal_alpha = 0.85,
  signal_bar_width = 1,
  frame_palette = "Paired",
  frame_colors = NULL,
  x_label = "Transcript coordinate",
  text_color = "black",
  text_size = 14,
  grid_linewidth = NULL,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL,
  highlight = NULL
) {
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)

  dt <- expand_bwg_to_positions(dt)
  frame_annotation <- data.table::as.data.table(frame_annotation)
  if (nrow(dt) == 0L) {
    dt <- data.table::data.table(
      sample_id = expected_samples,
      chrom = as.character(frame_annotation[["chrom"]][1L]),
      genomic_pos = as.integer(frame_annotation[["genomic_pos"]][1L]),
      start = as.integer(frame_annotation[["genomic_pos"]][1L]),
      end = as.integer(frame_annotation[["genomic_pos"]][1L]),
      value = 0,
      strand = "*"
    )
  }

  dt <- merge(
    dt,
    frame_annotation[, .(chrom, genomic_pos, plot_pos, region, frame)],
    by = c("chrom", "genomic_pos"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  if (nrow(dt) == 0L) {
    stop(
      "No signal records overlapped the exon positions of the selected transcript.",
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    as.character(expected_samples),
    unique(as.character(dt[["sample_id"]]))
  )
  if (length(missing_samples) > 0L) {
    filler <- data.table::CJ(
      sample_id = missing_samples,
      row_id = seq_len(nrow(frame_annotation)),
      unique = TRUE
    )
    filler <- merge(
      filler,
      frame_annotation[, .(
        row_id = seq_len(.N),
        chrom,
        genomic_pos,
        plot_pos,
        region,
        frame
      )],
      by = "row_id",
      all.x = TRUE,
      sort = FALSE
    )
    filler[, `:=`(
      start = as.integer(genomic_pos),
      end = as.integer(genomic_pos),
      value = 0,
      strand = "*"
    )]
    filler[, "row_id" := NULL]
    dt <- data.table::rbindlist(list(dt, filler), fill = TRUE)
  }

  dt[, `:=`(
    start = as.integer(plot_pos),
    end = as.integer(plot_pos),
    mid = as.numeric(plot_pos)
  )]

  dt <- apply_signal_grouping(
    dt,
    sample_groups = sample_groups,
    signal_summary = signal_summary
  )
  if (!"region" %in% names(dt) || !"frame" %in% names(dt)) {
    dt <- merge(
      dt,
      unique(frame_annotation[, .(
        chrom,
        start = plot_pos,
        end = plot_pos,
        region,
        frame
      )]),
      by = c("chrom", "start", "end"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  dt[, "mid" := as.numeric((start + end) / 2)]
  dt[, "plot_value" := transform_signal_value(value, signal_transform)]
  dt[,
    "frame" := factor(
      as.character(frame),
      levels = c("frame0", "frame1", "frame2")
    )
  ]

  sample_ids <- get_ordered_signal_ids(dt, "sample_id")
  dt[, "sample_id" := factor(as.character(sample_id), levels = sample_ids)]
  facet_scales <- if (signal_y_scale == "free") "free_y" else "fixed"
  y_scale <- make_signal_y_scale(
    values = dt[["plot_value"]],
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    signal_y_limits = signal_y_limits
  )
  y_anchor_dt <- make_signal_y_anchor(
    dt,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks
  )
  frame_cols <- make_frame_colors(
    frame_palette = frame_palette,
    frame_colors = frame_colors
  )

  y_label <- switch(
    signal_transform,
    none = "Signal",
    log2 = "Signal (signed log2(x + 1))",
    log10 = "Signal (signed log10(x + 1))",
    sqrt = "Signal (signed sqrt)"
  )

  base_theme <- make_track_theme(
    plot_theme = plot_theme,
    show_panel_border = show_panel_border
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.text.y = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.x = ggplot2::element_text(
        color = text_color,
        size = text_size
      ),
      axis.title.y = ggplot2::element_text(
        color = text_color,
        size = text_size
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(color = text_color, size = text_size),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = text_color, size = text_size),
      strip.text.y = ggplot2::element_text(
        color = text_color,
        size = text_size,
        angle = 0
      ),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  utr_dt <- dt[is.na(frame) | as.character(region) != "CDS"]
  cds_dt <- dt[!is.na(frame) & as.character(region) == "CDS"]
  p <- ggplot2::ggplot() +
    ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
    ggplot2::labs(x = x_label, y = y_label, fill = NULL) +
    y_scale +
    base_theme
  if (!is.null(y_anchor_dt)) {
    p <- p +
      ggplot2::geom_blank(
        data = y_anchor_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        inherit.aes = FALSE
      )
  }
  if (nrow(utr_dt) > 0L) {
    p <- p +
      ggplot2::geom_col(
        data = utr_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        fill = "grey80",
        color = NA,
        width = signal_bar_width,
        alpha = signal_alpha
      )
  }
  if (nrow(cds_dt) > 0L) {
    p <- p +
      ggplot2::geom_col(
        data = cds_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value, fill = .data$frame),
        width = signal_bar_width,
        alpha = signal_alpha
      ) +
      ggplot2::scale_fill_manual(
        values = frame_cols,
        breaks = c("frame0", "frame1", "frame2"),
        labels = c("frame0", "frame1", "frame2"),
        name = NULL,
        drop = FALSE
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))
  }
  add_highlight_layer(p, highlight)
}

plot_signal_core <- function(
  dt,
  plot_type = c("bar", "line", "heatmap"),
  highlight = NULL,
  x_label = "Genomic coordinate",
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL,
  sample_groups = NULL,
  signal_color_by = c("sample", "group"),
  signal_summary = c("none", "mean", "median", "sum"),
  signal_transform = c("none", "log2", "log10", "sqrt"),
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  signal_y_limits = NULL,
  signal_alpha = 0.85,
  signal_bar_width = 1,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL,
  text_size = 14,
  heatmap_bin_size = NULL,
  heatmap_max_bins = 800L,
  heatmap_summary = c("mean", "max", "sum", "median")
) {
  plot_type <- match.arg(plot_type)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  heatmap_summary <- match.arg(heatmap_summary)
  text_color <- "black"
  dt <- data.table::as.data.table(dt)
  stop_if_not(
    nrow(dt) > 0L,
    "No signal records were found in the specified region."
  )
  dt[, mid := (start + end) / 2]
  dt <- apply_signal_grouping(
    dt,
    sample_groups = sample_groups,
    signal_summary = signal_summary
  )
  dt[, mid := (start + end) / 2]
  dt[, plot_value := transform_signal_value(value, signal_transform)]
  sample_ids <- get_ordered_signal_ids(dt, "sample_id")
  color_ids <- if (
    signal_color_by == "group" && "sample_group" %in% names(dt)
  ) {
    get_ordered_signal_ids(dt, "sample_group")
  } else {
    sample_ids
  }
  dt[,
    "sample_id" := factor(as.character(dt[["sample_id"]]), levels = sample_ids)
  ]
  if ("sample_group" %in% names(dt)) {
    group_levels <- get_ordered_signal_ids(dt, "sample_group")
    dt[,
      "sample_group" := factor(
        as.character(dt[["sample_group"]]),
        levels = group_levels
      )
    ]
  }
  discrete_signal_colors <- normalize_signal_colors(
    sample_ids = color_ids,
    signal_palette = signal_palette,
    signal_palette_direction = signal_palette_direction,
    signal_colors = signal_colors
  )
  color_aes <- if (
    signal_color_by == "group" && "sample_group" %in% names(dt)
  ) {
    "sample_group"
  } else {
    "sample_id"
  }

  y_label <- switch(
    signal_transform,
    none = "Signal",
    log2 = "Signal (signed log2(x + 1))",
    log10 = "Signal (signed log10(x + 1))",
    sqrt = "Signal (signed sqrt)"
  )
  facet_scales <- if (signal_y_scale == "free") "free_y" else "fixed"
  y_scale <- make_signal_y_scale(
    values = dt[["plot_value"]],
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    signal_y_limits = signal_y_limits
  )
  y_anchor_dt <- make_signal_y_anchor(
    dt,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks
  )

  base_theme <- make_track_theme(
    plot_theme = plot_theme,
    show_panel_border = show_panel_border
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.text.y = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.x = ggplot2::element_text(
        color = text_color,
        size = text_size
      ),
      axis.title.y = ggplot2::element_text(
        color = text_color,
        size = text_size
      ),
      legend.position = "none",
      legend.text = ggplot2::element_text(color = text_color, size = text_size),
      legend.title = ggplot2::element_text(
        color = text_color,
        size = text_size
      ),
      strip.text = ggplot2::element_text(color = text_color, size = text_size),
      strip.text.y = ggplot2::element_text(
        color = text_color,
        size = text_size,
        angle = 0
      ),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (plot_type == "heatmap") {
    dt <- aggregate_signal_for_heatmap(
      dt,
      heatmap_bin_size = heatmap_bin_size,
      heatmap_max_bins = heatmap_max_bins,
      heatmap_summary = heatmap_summary
    )
    p <- ggplot2::ggplot(
      dt,
      ggplot2::aes(
        x = mid,
        y = sample_id,
        fill = plot_value,
        width = pmax(end - start + 1L, 1L)
      )
    ) +
      ggplot2::geom_tile(height = 0.90, alpha = signal_alpha) +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::labs(x = x_label, y = NULL, fill = y_label) +
      apply_signal_continuous_fill_scale(
        signal_palette = signal_palette,
        signal_palette_direction = signal_palette_direction,
        signal_colors = signal_colors
      ) +
      base_theme +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          title.position = "top",
          barwidth = grid::unit(4.5, "cm"),
          barheight = grid::unit(0.35, "cm")
        )
      ) +
      ggplot2::theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(
          color = text_color,
          size = text_size
        ),
        legend.text = ggplot2::element_text(
          color = text_color,
          size = text_size
        ),
        axis.text.y.right = ggplot2::element_text(
          color = text_color,
          size = text_size
        ),
        axis.text.y.left = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
    return(add_highlight_layer(p, highlight))
  }

  if (plot_type == "line") {
    dt <- fill_signal_line_gaps_with_zero(dt, color_aes = color_aes)
    p <- ggplot2::ggplot(
      dt,
      ggplot2::aes(
        x = mid,
        y = plot_value,
        color = .data[[color_aes]],
        group = .data$line_group
      )
    ) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, color = "Sample") +
      ggplot2::scale_color_manual(values = discrete_signal_colors) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p +
        ggplot2::geom_blank(
          data = y_anchor_dt,
          ggplot2::aes(x = .data$mid, y = .data$plot_value),
          inherit.aes = FALSE
        )
    }
    p <- p + ggplot2::geom_line(
      linewidth = 0.4,
      alpha = signal_alpha,
      na.rm = TRUE
    )
  } else if (plot_type == "bar") {
    dt[, "bar_ymin" := data.table::fifelse(plot_value >= 0, 0, plot_value)]
    dt[, "bar_ymax" := data.table::fifelse(plot_value >= 0, plot_value, 0)]
    dt[, "bar_width" := pmax(as.numeric(end) - as.numeric(start) + 1, 1)]
    dt[, "bar_xmid" := (as.numeric(start) + as.numeric(end)) / 2]
    dt[, "bar_xmin" := bar_xmid - bar_width * signal_bar_width / 2]
    dt[, "bar_xmax" := bar_xmid + bar_width * signal_bar_width / 2]
    p <- ggplot2::ggplot(dt, ggplot2::aes(fill = .data[[color_aes]])) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, fill = "Sample") +
      ggplot2::scale_fill_manual(values = discrete_signal_colors) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p +
        ggplot2::geom_blank(
          data = y_anchor_dt,
          ggplot2::aes(x = .data$mid, y = .data$plot_value),
          inherit.aes = FALSE
        )
    }
    p <- p +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$bar_xmin,
          xmax = .data$bar_xmax,
          ymin = .data$bar_ymin,
          ymax = .data$bar_ymax
        ),
        color = NA,
        alpha = signal_alpha
      )
  }

  add_highlight_layer(p, highlight)
}


fill_signal_line_gaps_with_zero <- function(dt, color_aes = "sample_id") {
  dt <- data.table::copy(data.table::as.data.table(dt))
  if (nrow(dt) == 0L) {
    dt[, "line_group" := character()]
    return(dt)
  }

  if (!"sample_id" %in% names(dt)) {
    dt[, "sample_id" := "sample"]
  }
  if (!color_aes %in% names(dt)) {
    color_aes <- "sample_id"
  }

  # Line tracks should not connect two non-adjacent covered intervals directly.
  # Instead of splitting the line into disconnected segments, add zero-valued
  # anchors at both sides of each uncovered gap. This makes the curve return to
  # the baseline across no-coverage regions, which is closer to genome-browser
  # coverage behavior and visually less abrupt than broken line segments.
  split_cols <- unique(c("sample_id", color_aes))
  data.table::setorderv(dt, c(split_cols, "start", "end"))

  expanded <- dt[,
    {
      d <- data.table::copy(.SD)
      s <- suppressWarnings(as.integer(d[["start"]]))
      e <- suppressWarnings(as.integer(d[["end"]]))
      prev_end <- data.table::shift(e, type = "lag")
      has_gap <- !is.na(prev_end) & !is.na(s) & s > (prev_end + 1L)
      gap_idx <- which(has_gap)

      if (length(gap_idx) == 0L) {
        d[,
          "line_group" := paste0(
            as.character(sample_id[1L]),
            "_",
            as.character(.GRP)
          )
        ]
        d
      } else {
        zero_rows <- vector("list", length(gap_idx) * 2L)
        k <- 1L
        for (idx in gap_idx) {
          prev_row <- data.table::copy(d[idx - 1L])
          next_row <- data.table::copy(d[idx])

          # The two zero anchors are placed at the uncovered interval boundaries.
          # They use half-base positions so the line drops after the previous
          # covered interval and rises immediately before the next covered interval.
          prev_row[, "start" := as.integer(prev_end[idx] + 1L)]
          prev_row[, "end" := as.integer(prev_end[idx] + 1L)]
          prev_row[, "mid" := as.numeric(prev_end[idx]) + 0.5]
          prev_row[, "value" := 0]
          prev_row[, "plot_value" := 0]

          next_row[, "start" := as.integer(s[idx] - 1L)]
          next_row[, "end" := as.integer(s[idx] - 1L)]
          next_row[, "mid" := as.numeric(s[idx]) - 0.5]
          next_row[, "value" := 0]
          next_row[, "plot_value" := 0]

          zero_rows[[k]] <- prev_row
          zero_rows[[k + 1L]] <- next_row
          k <- k + 2L
        }

        out <- data.table::rbindlist(
          c(list(d), zero_rows),
          fill = TRUE,
          use.names = TRUE
        )
        data.table::setorderv(out, c("mid", "start", "end"))
        out[,
          "line_group" := paste0(
            as.character(sample_id[1L]),
            "_",
            as.character(.GRP)
          )
        ]
        out
      }
    },
    by = split_cols
  ]

  expanded[]
}

make_signal_y_scale <- function(
  values = NULL,
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty"),
  signal_y_limits = NULL
) {
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)

  make_numeric_limits <- function(
    x,
    force_zero_baseline = TRUE,
    add_padding = TRUE
  ) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NULL)
    }

    lower <- min(x, na.rm = TRUE)
    upper <- max(x, na.rm = TRUE)

    if (isTRUE(force_zero_baseline)) {
      if (lower >= 0) {
        lower <- 0
      }
      if (upper <= 0) upper <- 0
    }

    if (!is.finite(lower) || !is.finite(upper)) {
      return(NULL)
    }

    if (identical(lower, upper)) {
      if (lower == 0) {
        upper <- 1
      } else if (lower > 0) {
        lower <- 0
      } else {
        upper <- 0
      }
    }

    if (isTRUE(add_padding)) {
      pad <- (upper - lower) * 0.03
      if (!is.finite(pad) || pad <= 0) {
        pad <- max(abs(upper), abs(lower), 1) * 0.03
      }
      if (upper >= 0) {
        upper <- upper + pad
      } else {
        lower <- lower - pad
      }
    }

    c(lower, upper)
  }

  make_range_breaks <- function(lim) {
    lim <- as.numeric(lim)
    lim <- lim[is.finite(lim)]
    if (length(lim) == 0L) {
      return(NULL)
    }
    lower <- min(lim, na.rm = TRUE)
    upper <- max(lim, na.rm = TRUE)
    if (!is.finite(lower) || !is.finite(upper)) {
      return(NULL)
    }
    if (identical(lower, upper)) {
      if (lower == 0) {
        upper <- 1
      } else if (lower > 0) {
        lower <- 0
      } else {
        upper <- 0
      }
    }
    unique(c(lower, upper))
  }

  signal_label <- function(x) {
    x <- as.numeric(x)
    vapply(
      x,
      function(v) {
        if (!is.finite(v)) {
          return(NA_character_)
        }
        if (abs(v) < .Machine$double.eps^0.5) {
          return("0")
        }
        if (abs(v) >= 1000) {
          return(formatC(v, format = "fg", digits = 4, big.mark = ","))
        }
        if (abs(v) < 0.001) {
          return(formatC(v, format = "e", digits = 2))
        }
        formatC(v, format = "fg", digits = 4, big.mark = ",")
      },
      character(1L)
    )
  }

  squish_to_limits <- function(x, range) {
    x <- as.numeric(x)
    x[x < range[1L]] <- range[1L]
    x[x > range[2L]] <- range[2L]
    x
  }

  if (signal_y_scale == "fixed") {
    lim <- if (is.null(signal_y_limits)) {
      make_numeric_limits(
        values,
        force_zero_baseline = TRUE,
        add_padding = TRUE
      )
    } else {
      signal_y_limits
    }
    if (is.null(lim)) {
      return(ggplot2::scale_y_continuous())
    }
    if (!is.null(signal_y_limits)) {
      if (signal_y_ticks == "range") {
        brks <- make_range_breaks(lim)
        return(ggplot2::scale_y_continuous(
          limits = lim,
          breaks = brks,
          labels = signal_label,
          expand = c(0, 0),
          oob = squish_to_limits
        ))
      }
      return(ggplot2::scale_y_continuous(
        limits = lim,
        labels = signal_label,
        expand = c(0, 0),
        oob = squish_to_limits
      ))
    }
    if (signal_y_ticks == "range") {
      brks <- make_range_breaks(lim)
      return(ggplot2::scale_y_continuous(
        limits = lim,
        breaks = brks,
        labels = signal_label,
        expand = c(0, 0)
      ))
    }
    return(ggplot2::scale_y_continuous(
      limits = lim,
      labels = signal_label,
      expand = c(0, 0)
    ))
  }

  if (signal_y_ticks == "pretty") {
    return(ggplot2::scale_y_continuous(
      labels = signal_label,
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ))
  }

  ggplot2::scale_y_continuous(
    breaks = function(x) make_range_breaks(x),
    labels = signal_label,
    expand = c(0, 0)
  )
}

make_signal_y_anchor <- function(
  dt,
  signal_y_scale = c("free", "fixed"),
  signal_y_ticks = c("range", "pretty")
) {
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  if (signal_y_scale != "free" || signal_y_ticks != "range") {
    return(NULL)
  }

  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L || !"plot_value" %in% names(dt)) {
    return(NULL)
  }

  anchor <- dt[,
    {
      values <- as.numeric(plot_value)
      values <- values[is.finite(values)]
      if (length(values) == 0L) {
        .(plot_value = numeric())
      } else {
        lower <- min(values, na.rm = TRUE)
        upper <- max(values, na.rm = TRUE)
        if (lower >= 0) {
          lower <- 0
        }
        if (upper <= 0) {
          upper <- 0
        }
        if (identical(lower, upper)) {
          if (lower == 0) {
            upper <- 1
          } else if (lower > 0) {
            lower <- 0
          } else {
            upper <- 0
          }
        }
        pad <- (upper - lower) * 0.03
        if (!is.finite(pad) || pad <= 0) {
          pad <- max(abs(upper), abs(lower), 1) * 0.03
        }
        if (upper >= 0) {
          upper <- upper + pad
        } else {
          lower <- lower - pad
        }
        .(plot_value = c(lower, upper))
      }
    },
    by = "sample_id"
  ]

  if (nrow(anchor) == 0L) {
    return(NULL)
  }

  mid_value <- stats::median(as.numeric(dt[["mid"]]), na.rm = TRUE)
  if (!is.finite(mid_value)) {
    mid_value <- 0
  }
  anchor[, "mid" := mid_value]
  anchor
}

transform_signal_value <- function(
  value,
  signal_transform = c("none", "log2", "log10", "sqrt")
) {
  signal_transform <- match.arg(signal_transform)
  value <- as.numeric(value)
  if (signal_transform == "none") {
    return(value)
  }
  if (signal_transform == "log2") {
    return(sign(value) * log2(abs(value) + 1))
  }
  if (signal_transform == "log10") {
    return(sign(value) * log10(abs(value) + 1))
  }
  sign(value) * sqrt(abs(value))
}

make_signal_palette <- function(
  n,
  signal_palette = "Paired",
  signal_palette_direction = 1
) {
  n <- max(1L, as.integer(n))
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  signal_palette <- as.character(signal_palette)[1L]
  if (is.na(signal_palette) || !nzchar(signal_palette)) {
    signal_palette <- "Paired"
  }

  if (
    requireNamespace("RColorBrewer", quietly = TRUE) &&
      signal_palette %in% rownames(RColorBrewer::brewer.pal.info)
  ) {
    pal_info <- RColorBrewer::brewer.pal.info[signal_palette, , drop = FALSE]
    max_colors <- as.integer(pal_info[["maxcolors"]][1L])

    # Discrete signal mappings must follow the published palette order. Ask
    # RColorBrewer for the required class count instead of interpolating from
    # the full palette, because interpolation makes two or three tracks jump
    # toward the palette endpoints rather than use colors in sequence.
    if (n <= max_colors) {
      request_n <- max(3L, n)
      cols <- RColorBrewer::brewer.pal(request_n, signal_palette)
      if (signal_palette_direction == -1L) {
        cols <- rev(cols)
      }
      return(cols[seq_len(n)])
    }

    base <- RColorBrewer::brewer.pal(max_colors, signal_palette)
    if (signal_palette_direction == -1L) {
      base <- rev(base)
    }
    return(grDevices::colorRampPalette(base)(n))
  }

  predefined <- list(
    Paired = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
    Blues = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
    Reds = c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
    RdBu = c("#B2182B", "#EF8A62", "#FDDDBC", "#D1E5F0", "#67A9CF", "#2166AC")
  )
  base <- predefined[[signal_palette]]
  if (is.null(base)) {
    warning(
      sprintf(
        "Unknown `signal_palette`: %s. Falling back to 'Paired'.",
        signal_palette
      ),
      call. = FALSE
    )
    base <- predefined[["Paired"]]
  }
  if (signal_palette_direction == -1L) {
    base <- rev(base)
  }
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::colorRampPalette(base)(n)
}

make_signal_continuous_palette <- function(
  n = 256L,
  signal_palette = "Paired",
  signal_palette_direction = 1
) {
  n <- max(2L, as.integer(n))
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  signal_palette <- as.character(signal_palette)[1L]
  if (is.na(signal_palette) || !nzchar(signal_palette)) {
    signal_palette <- "Paired"
  }

  if (
    requireNamespace("RColorBrewer", quietly = TRUE) &&
      signal_palette %in% rownames(RColorBrewer::brewer.pal.info)
  ) {
    max_colors <- as.integer(
      RColorBrewer::brewer.pal.info[signal_palette, "maxcolors"]
    )
    base <- RColorBrewer::brewer.pal(max_colors, signal_palette)
  } else {
    predefined <- list(
      Paired = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
      Blues = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
      Reds = c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
      RdBu = c("#B2182B", "#EF8A62", "#FDDDBC", "#D1E5F0", "#67A9CF", "#2166AC")
    )
    base <- predefined[[signal_palette]]
    if (is.null(base)) {
      warning(
        sprintf(
          "Unknown `signal_palette`: %s. Falling back to 'Paired'.",
          signal_palette
        ),
        call. = FALSE
      )
      base <- predefined[["Paired"]]
    }
  }

  if (signal_palette_direction == -1L) {
    base <- rev(base)
  }
  grDevices::colorRampPalette(base)(n)
}

normalize_palette_direction <- function(direction = 1) {
  direction <- suppressWarnings(as.integer(direction[1L]))
  if (is.na(direction) || !direction %in% c(1L, -1L)) {
    warning(
      "`signal_palette_direction` must be 1 or -1. Falling back to 1.",
      call. = FALSE
    )
    direction <- 1L
  }
  direction
}

normalize_signal_colors <- function(
  sample_ids,
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL
) {
  sample_ids <- unique(as.character(sample_ids))
  n <- length(sample_ids)
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  if (!is.null(signal_colors)) {
    cols <- as.character(signal_colors)
    if (!is.null(names(cols)) && all(sample_ids %in% names(cols))) {
      return(cols[sample_ids])
    }
    if (length(cols) < n) {
      cols <- rep(cols, length.out = n)
    }
    cols <- cols[seq_len(n)]
    names(cols) <- sample_ids
    return(cols)
  }
  cols <- make_signal_palette(
    n,
    signal_palette,
    signal_palette_direction = signal_palette_direction
  )
  names(cols) <- sample_ids
  cols
}

apply_signal_continuous_fill_scale <- function(
  signal_palette = "Paired",
  signal_palette_direction = 1,
  signal_colors = NULL
) {
  signal_palette_direction <- normalize_palette_direction(
    signal_palette_direction
  )
  if (!is.null(signal_colors)) {
    cols <- as.character(signal_colors)
  } else {
    cols <- make_signal_continuous_palette(
      256L,
      signal_palette = signal_palette,
      signal_palette_direction = signal_palette_direction
    )
  }
  ggplot2::scale_fill_gradientn(colors = cols)
}

get_ordered_signal_ids <- function(dt, column = "sample_id") {
  x <- dt[[column]]
  x_chr <- as.character(x)
  if (is.factor(x)) {
    lev <- levels(x)
    lev <- lev[lev %in% x_chr]
    return(lev)
  }
  unique(x_chr)
}

normalize_sample_groups <- function(sample_ids, sample_groups = NULL) {
  sample_ids <- unique(as.character(sample_ids))
  if (is.null(sample_groups)) {
    return(stats::setNames(sample_ids, sample_ids))
  }
  if (is.data.frame(sample_groups)) {
    stop_if_not(
      all(c("sample_id", "group") %in% names(sample_groups)),
      "`sample_groups` data frame must contain `sample_id` and `group` columns."
    )
    map <- stats::setNames(
      as.character(sample_groups[["group"]]),
      as.character(sample_groups[["sample_id"]])
    )
    missing <- setdiff(sample_ids, names(map))
    if (length(missing) > 0L) {
      map[missing] <- missing
    }
    return(map[sample_ids])
  }
  groups <- as.character(sample_groups)
  if (!is.null(names(groups)) && all(sample_ids %in% names(groups))) {
    return(groups[sample_ids])
  }
  stop_if_not(
    length(groups) == length(sample_ids),
    "Unnamed `sample_groups` must have the same length as selected samples."
  )
  stats::setNames(groups, sample_ids)
}

apply_signal_grouping <- function(
  dt,
  sample_groups = NULL,
  signal_summary = c("none", "mean", "median", "sum")
) {
  signal_summary <- match.arg(signal_summary)
  dt <- data.table::copy(data.table::as.data.table(dt))
  sample_ids <- get_ordered_signal_ids(dt, "sample_id")
  group_map <- normalize_sample_groups(sample_ids, sample_groups)
  group_levels <- unique(as.character(group_map[sample_ids]))
  dt[,
    "sample_group" := factor(
      as.character(group_map[as.character(dt[["sample_id"]])]),
      levels = group_levels
    )
  ]
  if (signal_summary == "none") {
    return(dt[])
  }
  summary_fun <- switch(
    signal_summary,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
  )
  out <- dt[,
    .(
      value = summary_fun(value),
      strand = as.character(strand[1L])
    ),
    by = .(sample_group, chrom, start, end)
  ]
  out[, "sample_id" := as.character(out[["sample_group"]])]
  out[,
    "sample_id" := factor(
      as.character(out[["sample_id"]]),
      levels = group_levels
    )
  ]
  out[,
    "sample_group" := factor(
      as.character(out[["sample_group"]]),
      levels = group_levels
    )
  ]
  data.table::setcolorder(
    out,
    c("sample_id", "sample_group", "chrom", "start", "end", "value", "strand")
  )
  out[]
}

map_signal_to_exons <- function(signal_dt, exons) {
  dt <- data.table::as.data.table(signal_dt)
  ex <- data.table::copy(exons)
  if (nrow(dt) == 0L || nrow(ex) == 0L) {
    return(dt[0])
  }

  data.table::setorder(ex, transcript_id, exon_start, exon_end)
  ex[, exon_width := exon_end - exon_start + 1L]
  ex[,
    exon_offset := cumsum(data.table::shift(exon_width, fill = 0L)),
    by = transcript_id
  ]
  ex <- ex[, .(chrom, exon_start, exon_end, exon_offset)]

  mapped <- list()
  for (i in seq_len(nrow(ex))) {
    piece <- dt[
      chrom == ex$chrom[i] & start <= ex$exon_end[i] & end >= ex$exon_start[i]
    ]
    if (nrow(piece) == 0L) {
      next
    }
    piece[, `:=`(
      start = pmax(start, ex$exon_start[i]),
      end = pmin(end, ex$exon_end[i])
    )]
    piece[, `:=`(
      start = ex$exon_offset[i] + (start - ex$exon_start[i]) + 1L,
      end = ex$exon_offset[i] + (end - ex$exon_start[i]) + 1L
    )]
    mapped[[length(mapped) + 1L]] <- piece
  }
  data.table::rbindlist(mapped, fill = TRUE)
}
