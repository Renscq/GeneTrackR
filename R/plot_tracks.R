# Author: Rensc
# Date: 2026-08-31
# Version: dev009
# Function: Combined genome-browser-like track plotting
# Input: GenePred and BwgTrack objects
# Output: Combined patchwork track figure

#' Plot combined signal and gene model tracks
#'
#' @description
#' Draw a genome-browser-like figure containing a gene model track alone or a
#' signal track combined with gene, feature, and variant tracks. The genomic
#' interval can be specified in three ways: by `gene_id`, by `transcript_id`,
#' or by explicit `chrom`, `start`, and `end` coordinates.
#'
#' @param annotation A GenePred object or a standardized Feature object with transcript/exon records.
#' @param signal Optional BwgTrack object. If NULL, only annotation/feature/variant tracks are drawn.
#' @param features Optional FeatureTrack object or named list of FeatureTrack objects from `read_bed()`, `read_gff()`, or `read_gtf()`.
#' @param variants Optional VariantTrack object or named list of VariantTrack objects from `read_vcf()`.
#' @param chrom Chromosome name. Required when `gene_id` and `transcript_id` are not supplied.
#' @param start Region start in 1-based closed coordinates. Required with `chrom`/`end`.
#' @param end Region end in 1-based closed coordinates. Required with `chrom`/`start`.
#' @param gene_id Optional gene ID. If supplied, the plotting interval is inferred from the gene locus.
#' @param transcript_id Optional transcript ID. If supplied, the plotting interval is inferred from the transcript locus.
#' @param samples Optional signal sample IDs.
#' @param sample_groups Optional sample group mapping for group-level coloring or replicate summaries. Use a named character vector, a data frame with `sample_id` and `group`, or an unnamed vector with one group per selected sample.
#' @param signal_color_by Color signal tracks by `sample` or `group`.
#' @param signal_summary Replicate summary mode. Use `none` to plot individual samples, or `mean`, `median`, or `sum` to summarize samples within each group.
#' @param signal_type Signal plot type for non-Ribo-seq browser signal tracks.
#' @param signal_palette Signal color palette for signal tracks. Any palette name from `RColorBrewer::brewer.pal.info` can be used. Discrete signal sample/group colors follow the standard RColorBrewer class order; heatmaps use the corresponding continuous gradient.
#' @param signal_palette_direction Direction for generated signal colors. Use `1` for the standard palette order and `-1` for the reversed palette order. Discrete sample/group colors preserve level-to-color order; heatmap gradients reverse continuously.
#' @param signal_colors Optional named or unnamed vector of explicit colors for signal samples. Explicit colors override `signal_palette` and are not modified by `signal_palette_direction`.
#' @param gene_palette RColorBrewer palette name used for gene model feature fills.
#' @param gene_colors Optional custom fill colors for gene model features. Use a named vector such as `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")`.
#' @param gene_border_color Optional rectangle border color for gene model features. Use `NA` to hide borders.
#' @param feature_color_by Feature-track color grouping. `auto` chooses a compact grouping automatically; other common choices are `feature_group`, `type`, `source`, `name`, and `strand`.
#' @param feature_palette RColorBrewer palette name used for feature-track fills.
#' @param feature_colors Optional explicit feature-track fill colors.
#' @param feature_border_color Optional rectangle border color for feature tracks. Use `NA` to hide borders.
#' @param feature_max_legend_levels Maximum number of legend groups shown for each feature track.
#' @param variant_palette RColorBrewer palette name used for variant-track colors.
#' @param variant_colors Optional named or unnamed explicit colors passed to `plot_variant()`. Named vectors are recommended when stable SNP/INS/DEL colors are required.
#' @param ribo_signal_type How Ribo-seq samples should be displayed in transcript-centered browser plots. `auto` uses `frame` for samples whose IDs look like Ribo-seq/RPF tracks when `transcript_id` is supplied, otherwise it falls back to standard genomic signal tracks. `bar` always uses the standard signal geometry, and `frame` forces frame rendering when transcript-centered plotting is possible.
#' @param frame_palette RColorBrewer palette name used for Ribo-seq frame plots.
#' @param frame_colors Optional explicit colors for `frame0`, `frame1`, and `frame2`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for independent sample-specific y-axis ranges or `fixed` for a shared y-axis range across samples.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer minimum and maximum limits or `pretty` for default-style breaks.
#' @param signal_y_limits Optional two-element numeric vector giving the plotted y-axis limits after `signal_transform`. Supplying limits changes `signal_y_scale` to `fixed`.
#' @param signal_alpha Signal geometry transparency from 0 to 1.
#' @param signal_bar_width Relative width of bar intervals from greater than 0 to 1. A value below 1 creates proportional gaps without changing interval centers.
#' @param collapse Gene model collapse mode for region-level plotting.
#' @param strand Signal strand selector.
#' @param bin_size Optional signal bin size.
#' @param heatmap_bin_size Optional bin size used only for heatmap rendering.
#' `NULL` uses adaptive binning.
#' @param heatmap_max_bins Maximum number of heatmap bins when
#' `heatmap_bin_size = NULL`.
#' @param heatmap_summary Summary statistic used to aggregate signal into
#' heatmap bins: `mean`, `max`, `sum`, or `median`.
#' @param highlight Optional data frame used to shade intervals on signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates. Optional columns are allowed but ignored by the default renderer.
#' @param layout Track layout. Use `signal_top` to place signal above gene model, or `gene_top` to place gene model above signal.
#' @param heights Relative panel heights. Must contain at least `signal`, `gene`, `feature`, and `variant` names when those tracks are used.
#' @param cds_height Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_height Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws one short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis and `feature` draws labels near the feature.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param plot_theme Base ggplot2 theme used by all standard track panels. Use `bw`, `classic`, `light`, or `minimal`.
#' @param show_panel_border Whether to draw panel borders. `NULL` preserves the selected theme default.
#' @param text_size Text size in points for axis text, axis titles, legends, and facet labels.
#'
#' @return A patchwork object or ggplot object.
#' @examples
#' gp_file <- system.file(
#'   "extdata", "gtr_demo.genePredExt", package = "GeneTrackR"
#' )
#' gp <- read_genepred(
#'   gp_file, format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' plot_tracks(annotation = gp, gene_id = "GeneA")
#' @export
plot_tracks <- function(annotation,
                        signal = NULL,
                        features = NULL,
                        variants = NULL,
                        chrom = NULL,
                        start = NULL,
                        end = NULL,
                        gene_id = NULL,
                        transcript_id = NULL,
                        samples = NULL,
                        sample_groups = NULL,
                        signal_color_by = c("sample", "group"),
                        signal_summary = c("none", "mean", "median", "sum"),
                        signal_type = c("bar", "line", "heatmap"),
                        signal_palette = "Paired",
                        signal_palette_direction = 1,
                        signal_colors = NULL,
                        gene_palette = "Paired",
                        gene_colors = NULL,
                        gene_border_color = NA,
                        feature_color_by = "auto",
                        feature_palette = "Paired",
                        feature_colors = NULL,
                        feature_border_color = NA,
                        feature_max_legend_levels = 5,
                        variant_palette = "Paired",
                        variant_colors = NULL,
                        ribo_signal_type = c("auto", "bar", "frame"),
                        frame_palette = "Paired",
                        frame_colors = NULL,
                        signal_transform = c("none", "log2", "log10", "sqrt"),
                        signal_y_scale = c("free", "fixed"),
                        signal_y_ticks = c("range", "pretty"),
                        heatmap_bin_size = NULL,
                        heatmap_max_bins = 800L,
                        heatmap_summary = c("mean", "max", "sum", "median"),
                        collapse = c("none", "union_exon", "longest"),
                        strand = c("ignore", "+", "-", "both"),
                        bin_size = NULL,
                        highlight = NULL,
                        layout = c("signal_top", "gene_top"),
                        heights = c(signal = 3, gene = 1, feature = 0.8, variant = 0.7),
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
                        show_panel_border = NULL) {
  stop_if_not(is_gene_model_feature(annotation), "`annotation` must be a GenePred object or a Feature object with transcript/exon records.")
  annotation <- as_genepred(annotation)
  signal_type <- match.arg(signal_type)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  collapse <- match.arg(collapse)
  strand <- match.arg(strand)
  layout <- match.arg(layout)
  ribo_signal_type <- match.arg(ribo_signal_type)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  signal_y_limits <- normalize_signal_y_limits(signal_y_limits)
  signal_y_scale <- resolve_signal_y_scale(signal_y_scale, signal_y_limits)
  signal_alpha <- normalize_signal_alpha(signal_alpha)
  signal_bar_width <- normalize_signal_bar_width(signal_bar_width)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  heatmap_summary <- match.arg(heatmap_summary)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  signal_palette_direction <- normalize_palette_direction(signal_palette_direction)
  gene_border_color <- normalize_border_color(gene_border_color)
  feature_border_color <- normalize_border_color(feature_border_color)

  has_gene_id <- !is.null(gene_id)
  has_transcript_id <- !is.null(transcript_id)
  has_region <- !is.null(chrom) || !is.null(start) || !is.null(end)
  n_locator <- sum(c(has_gene_id, has_transcript_id, has_region))
  stop_if_not(n_locator == 1L, "Specify exactly one locator: `gene_id`, `transcript_id`, or `chrom` + `start` + `end`.")

  tx_all <- data.table::as.data.table(annotation$transcripts)
  transcript_id_value <- NULL

  if (has_gene_id) {
    gene_id_value <- as.character(gene_id)[1L]
    tx <- tx_all[tx_all[["gene_id"]] == gene_id_value]
    stop_if_not(nrow(tx) > 0L, "Gene ID was not found.")
    gene <- build_gene_table(tx)
    chrom_value <- as.character(gene[["chrom"]][1L])
    start_value <- as.integer(gene[["gene_start"]][1L])
    end_value <- as.integer(gene[["gene_end"]][1L])
    p_gene <- plot_gene(
      annotation,
      gene_id = gene_id_value,
      collapse = collapse,
      coordinate = "genomic",
      highlight = highlight,
      cds_height = cds_height,
      utr_height = utr_height,
      direction_mode = direction_mode,
      gene_palette = gene_palette,
      gene_colors = gene_colors,
      gene_border_color = gene_border_color,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      label_position = label_position,
      label_by = label_by,
      text_size = text_size
    )
  } else if (has_transcript_id) {
    transcript_id_value <- as.character(transcript_id)[1L]
    tx <- tx_all[tx_all[["transcript_id"]] == transcript_id_value]
    stop_if_not(nrow(tx) > 0L, "Transcript ID was not found.")
    chrom_value <- as.character(tx[["chrom"]][1L])
    start_value <- as.integer(tx[["tx_start"]][1L])
    end_value <- as.integer(tx[["tx_end"]][1L])
    p_gene <- plot_transcript(
      annotation,
      transcript_id = transcript_id_value,
      coordinate = "genomic",
      highlight = highlight,
      cds_height = cds_height,
      utr_height = utr_height,
      direction_mode = direction_mode,
      gene_palette = gene_palette,
      gene_colors = gene_colors,
      gene_border_color = gene_border_color,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      label_position = label_position,
      label_by = label_by,
      text_size = text_size
    )
  } else {
    stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end), "`chrom`, `start`, and `end` are required for region-level plotting.")
    check_region(chrom, start, end)
    chrom_value <- as.character(chrom)[1L]
    start_value <- as.integer(start)[1L]
    end_value <- as.integer(end)[1L]
    p_gene <- plot_region(
      annotation,
      chrom = chrom_value,
      start = start_value,
      end = end_value,
      mode = "overlap",
      collapse = collapse,
      highlight = highlight,
      cds_height = cds_height,
      utr_height = utr_height,
      direction_mode = direction_mode,
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

  plot_list <- list()
  height_list <- numeric()

  if (!is.null(signal)) {
    stop_if_not(inherits(signal, "BwgTrack"), "`signal` must be a BwgTrack object.")
    p_signal <- make_browser_signal_plot(
      signal = signal,
      annotation = annotation,
      chrom = chrom_value,
      start = start_value,
      end = end_value,
      transcript_id = transcript_id_value,
      samples = samples,
      sample_groups = sample_groups,
      signal_color_by = signal_color_by,
      signal_summary = signal_summary,
      signal_type = signal_type,
      ribo_signal_type = ribo_signal_type,
      strand = strand,
      bin_size = bin_size,
      highlight = highlight,
      signal_palette = signal_palette,
      signal_palette_direction = signal_palette_direction,
      signal_colors = signal_colors,
      frame_palette = frame_palette,
      frame_colors = frame_colors,
      signal_transform = signal_transform,
      signal_y_scale = signal_y_scale,
      signal_y_ticks = signal_y_ticks,
      signal_y_limits = signal_y_limits,
      signal_alpha = signal_alpha,
      signal_bar_width = signal_bar_width,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      heatmap_bin_size = heatmap_bin_size,
      heatmap_max_bins = heatmap_max_bins,
      heatmap_summary = heatmap_summary,
      text_size = text_size
    )
    plot_list$signal <- p_signal
    height_list <- c(height_list, signal = get_track_height(heights, "signal", 3))
  }

  feature_plots <- make_feature_track_plots(
    features = features,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
    color_by = feature_color_by,
    feature_palette = feature_palette,
    feature_colors = feature_colors,
    feature_border_color = feature_border_color,
    feature_max_legend_levels = feature_max_legend_levels,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    text_size = text_size
  )
  if (length(feature_plots) > 0L) {
    for (nm in names(feature_plots)) {
      plot_list[[paste0("feature_", nm)]] <- feature_plots[[nm]]
      height_list <- c(height_list, feature = get_track_height(heights, "feature", 0.8))
    }
  }

  variant_plots <- make_variant_track_plots(
    variants = variants,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    text_size = text_size
  )
  if (length(variant_plots) > 0L) {
    for (nm in names(variant_plots)) {
      plot_list[[paste0("variant_", nm)]] <- variant_plots[[nm]]
      height_list <- c(height_list, variant = get_track_height(heights, "variant", 0.7))
    }
  }

  if (layout == "signal_top") {
    plot_list$gene <- p_gene
    height_list <- c(height_list, gene = get_track_height(heights, "gene", 1))
  } else {
    plot_list <- c(list(gene = p_gene), plot_list)
    height_list <- c(gene = get_track_height(heights, "gene", 1), height_list)
  }

  if (length(plot_list) == 1L) {
    return(plot_list[[1L]])
  }
  patchwork::wrap_plots(plot_list, ncol = 1) + patchwork::plot_layout(heights = unname(height_list))
}

make_browser_signal_plot <- function(signal,
                                     annotation,
                                     chrom,
                                     start,
                                     end,
                                     transcript_id = NULL,
                                     samples = NULL,
                                     sample_groups = NULL,
                                     signal_color_by = c("sample", "group"),
                                     signal_summary = c("none", "mean", "median", "sum"),
                                     signal_type = c("bar", "line", "heatmap"),
                                     ribo_signal_type = c("auto", "bar", "frame"),
                                     strand = c("ignore", "+", "-", "both"),
                                     bin_size = NULL,
                                     highlight = NULL,
                                     signal_palette = "Paired",
                                     signal_palette_direction = 1,
                                     signal_colors = NULL,
                                     frame_palette = "Paired",
                                     frame_colors = NULL,
                                     signal_transform = c("none", "log2", "log10", "sqrt"),
                                     signal_y_scale = c("free", "fixed"),
                                     signal_y_ticks = c("range", "pretty"),
                                     signal_y_limits = NULL,
                                     signal_alpha = 0.85,
                                     signal_bar_width = 1,
                                     plot_theme = c("bw", "classic", "light", "minimal"),
                                     show_panel_border = NULL,
                                     heatmap_bin_size = NULL,
                                     heatmap_max_bins = 800L,
                                     heatmap_summary = c("mean", "max", "sum", "median"),
                                     text_size = 14) {
  signal_type <- match.arg(signal_type)
  ribo_signal_type <- match.arg(ribo_signal_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_ticks <- match.arg(signal_y_ticks)
  plot_theme <- normalize_plot_theme(plot_theme)
  show_panel_border <- normalize_show_panel_border(show_panel_border)
  heatmap_summary <- match.arg(heatmap_summary)

  region_plot <- plot_signal_region(
    signal = signal,
    chrom = chrom,
    start = start,
    end = end,
    samples = samples,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    plot_type = signal_type,
    strand = strand,
    bin_size = bin_size,
    highlight = highlight,
    annotation = NULL,
    show_gene_model = FALSE,
    signal_palette = signal_palette,
    signal_palette_direction = signal_palette_direction,
    signal_colors = signal_colors,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    signal_y_limits = signal_y_limits,
    signal_alpha = signal_alpha,
    signal_bar_width = signal_bar_width,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary,
    text_size = text_size
  )

  if (is.null(transcript_id)) {
    if (identical(ribo_signal_type, "frame")) {
      message("[GeneTrackR] `ribo_signal_type = 'frame'` requires `transcript_id`; falling back to standard genomic browser signal tracks.")
    }
    return(region_plot)
  }

  selected_samples <- get_expected_signal_samples(
    signal,
    samples = samples,
    strand = strand
  )
  ribo_samples <- detect_ribo_signal_samples(selected_samples)
  if (length(ribo_samples) == 0L || identical(ribo_signal_type, "bar")) {
    return(region_plot)
  }

  non_ribo_samples <- setdiff(selected_samples, ribo_samples)
  plot_list <- list()
  height_list <- numeric()

  if (length(non_ribo_samples) > 0L) {
    plot_list$signal_non_ribo <- plot_signal_transcript(
      signal = signal,
      annotation = annotation,
      transcript_id = transcript_id,
      samples = non_ribo_samples,
      sample_groups = sample_groups,
      signal_color_by = signal_color_by,
      signal_summary = signal_summary,
      coordinate = "genomic",
      plot_type = signal_type,
      strand = strand,
      bin_size = bin_size,
      highlight = highlight,
      show_gene_model = FALSE,
      signal_palette = signal_palette,
      signal_palette_direction = signal_palette_direction,
      signal_colors = signal_colors,
      signal_transform = signal_transform,
      signal_y_scale = signal_y_scale,
      signal_y_ticks = signal_y_ticks,
      heatmap_bin_size = heatmap_bin_size,
      heatmap_max_bins = heatmap_max_bins,
      heatmap_summary = heatmap_summary,
      text_size = text_size,
      signal_y_limits = signal_y_limits,
      signal_alpha = signal_alpha,
      signal_bar_width = signal_bar_width,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border
    )
    height_list <- c(height_list, max(1, length(non_ribo_samples)))
  }

  plot_list$signal_ribo <- plot_signal_transcript(
    signal = signal,
    annotation = annotation,
    transcript_id = transcript_id,
    samples = ribo_samples,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    coordinate = "genomic",
    plot_type = "frame",
    strand = strand,
    bin_size = NULL,
    highlight = highlight,
    show_gene_model = FALSE,
    signal_palette = signal_palette,
    signal_palette_direction = signal_palette_direction,
    signal_colors = signal_colors,
    frame_palette = frame_palette,
    frame_colors = frame_colors,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary,
    text_size = text_size,
    signal_y_limits = signal_y_limits,
    signal_alpha = signal_alpha,
    signal_bar_width = signal_bar_width,
    plot_theme = plot_theme,
    show_panel_border = show_panel_border
  )
  height_list <- c(height_list, max(1, length(ribo_samples)))

  if (length(plot_list) == 1L) {
    return(plot_list[[1L]])
  }
  patchwork::wrap_plots(plot_list, ncol = 1) + patchwork::plot_layout(heights = unname(height_list))
}

detect_ribo_signal_samples <- function(sample_ids) {
  sample_ids <- unique(as.character(sample_ids))
  sample_ids[grepl("ribo|rpf|footprint", sample_ids, ignore.case = TRUE)]
}

make_feature_track_plots <- function(features,
                                     chrom,
                                     start,
                                     end,
                                     color_by = "auto",
                                     feature_palette = "Paired",
                                     feature_colors = NULL,
                                     feature_border_color = NA,
                                     feature_max_legend_levels = 5,
                                     plot_theme = "bw",
                                     show_panel_border = NULL,
                                     text_size = 14) {
  if (is.null(features)) return(list())
  if (inherits(features, "FeatureTrack")) features <- list(Feature = features)
  stop_if_not(is.list(features), "`features` must be a FeatureTrack object or a list of FeatureTrack objects.")
  if (is.null(names(features)) || any(!nzchar(names(features)))) names(features) <- paste0("Feature", seq_along(features))
  out <- list()
  for (nm in names(features)) {
    stop_if_not(inherits(features[[nm]], "FeatureTrack"), "All `features` entries must be FeatureTrack objects.")
    out[[nm]] <- plot_feature_track(
      features[[nm]],
      chrom = chrom,
      start = start,
      end = end,
      mode = "trim",
      color_by = color_by,
      label_by = "none",
      feature_palette = feature_palette,
      feature_colors = feature_colors,
      feature_border_color = feature_border_color,
      max_legend_levels = feature_max_legend_levels,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      text_size = text_size
    )
  }
  out
}

make_variant_track_plots <- function(variants, chrom, start, end, variant_palette = "Paired", variant_colors = NULL, plot_theme = "bw", show_panel_border = NULL, text_size = 14) {
  if (is.null(variants)) return(list())
  if (inherits(variants, "VariantTrack")) variants <- list(Variant = variants)
  stop_if_not(is.list(variants), "`variants` must be a VariantTrack object or a list of VariantTrack objects.")
  if (is.null(names(variants)) || any(!nzchar(names(variants)))) names(variants) <- paste0("Variant", seq_along(variants))
  out <- list()
  for (nm in names(variants)) {
    stop_if_not(inherits(variants[[nm]], "VariantTrack"), "All `variants` entries must be VariantTrack objects.")
    out[[nm]] <- plot_variant(
      variants[[nm]],
      chrom = chrom,
      start = start,
      end = end,
      label_by = "none",
      variant_palette = variant_palette,
      variant_colors = variant_colors,
      plot_theme = plot_theme,
      show_panel_border = show_panel_border,
      text_size = text_size
    )
  }
  out
}

get_track_height <- function(heights, name, default) {
  if (is.null(heights) || is.null(names(heights)) || !name %in% names(heights)) {
    return(as.numeric(default))
  }
  val <- suppressWarnings(as.numeric(heights[[name]]))
  if (!is.finite(val) || val <= 0) {
    return(as.numeric(default))
  }
  val
}
