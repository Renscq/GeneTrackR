# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.32
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
#' @param plot_type Signal plot type: `bar`, `line`, `area`, `heatmap`, or `frame_bar`. `frame_bar` colors per-base CDS signal by CDS frame and is useful for Ribo-seq tracks.
#' @param strand Strand selector. Use `auto` to use the transcript strand.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns. For `coordinate = "genomic"`, these are genomic coordinates; for `coordinate = "transcript"`, these are spliced transcript coordinates.
#' @param show_gene_model Whether to append the transcript gene model track. Default TRUE.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. If the number of samples or groups exceeds the palette maximum, colors are automatically interpolated.
#' @param signal_colors Optional named or unnamed vector of colors for samples. If supplied, it overrides `signal_palette`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param grid_linewidth Grid line width in the signal panel.
#' @param cds_width Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_width Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param show_direction Whether to draw direction arrows in the gene model track.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, and `end` draws one short arrow at the directional end of each gene.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis, `feature` draws labels on the model, and `none` hides labels.
#' @param label_by Which identifier to use for gene model labels.
#' @param label_side Label placement when `label_position = "feature"`. Use `above`, `below`, or `center`.
#' @param label_offset Vertical offset used for feature labels when `label_side` is `above` or `below`.
#' @param text_color Text color for signal and gene model text.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @param label_size Text size in points for labels drawn in the gene model track.
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
#' groups <- c(sampleA = "WT", sampleB = "KO")
#' plot_signal_transcript(
#'   signal = bg,
#'   annotation = gp,
#'   transcript_id = "TxA1",
#'   samples = c("sampleA", "sampleB"),
#'   sample_groups = groups,
#'   signal_color_by = "group",
#'   signal_summary = "mean",
#'   signal_palette = "Set1",
#'   signal_transform = "sqrt",
#'   signal_y_scale = "fixed",
#'   signal_y_ticks = "range",
#'   show_gene_model = TRUE,
#'   highlight = data.frame(start = 100, end = 250)
#' )
#' }
#' @export
plot_signal_transcript <- function(signal, annotation, transcript_id, samples = NULL, sample_groups = NULL, signal_color_by = c("sample", "group"), signal_summary = c("none", "mean", "median", "sum"), coordinate = c("transcript", "genomic"), plot_type = c("bar", "line", "area", "heatmap", "frame_bar"), strand = c("auto", "+", "-", "both", "ignore"), bin_size = NULL, highlight = NULL, show_gene_model = TRUE, signal_palette = "Blues", signal_colors = NULL, signal_transform = c("none", "log2", "log10", "sqrt"), signal_y_scale = c("free", "fixed"), signal_y_ticks = c("range", "pretty"), grid_linewidth = 0.25, heatmap_bin_size = NULL, heatmap_max_bins = 800L, heatmap_summary = c("mean", "max", "sum", "median"), frame_palette = "Set1", frame_colors = NULL, show_frame_legend = TRUE, cds_width = 0.50, utr_width = 0.25, show_direction = TRUE, direction_mode = c("transcript", "gene", "end"), label_position = c("axis", "feature", "none"), label_by = c("gene", "transcript"), label_side = c("above", "below", "center"), label_offset = 0.45, text_color = "black", text_size = 14, label_size = 12) {
  stop_if_not(inherits(signal, "BwgTrack"), "`signal` must be a BwgTrack object.")
  stop_if_not(inherits(annotation, "GenePred"), "`annotation` must be a GenePred object.")
  coordinate <- match.arg(coordinate)
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  label_side <- match.arg(label_side)

  transcript_id_value <- as.character(transcript_id)
  tx_all <- annotation$transcripts
  ex_all <- annotation$exons
  tx <- tx_all[tx_all[["transcript_id"]] == transcript_id_value]
  stop_if_not(nrow(tx) > 0L, "Transcript ID was not found.")
  ex <- ex_all[ex_all[["transcript_id"]] == transcript_id_value]
  selected_strand <- if (strand == "auto") tx$strand[1] else strand
  expected_samples <- get_expected_signal_samples(signal, samples = samples, strand = selected_strand)

  dt <- retrieve_bwg(signal, tx$chrom[1], tx$tx_start[1], tx$tx_end[1], samples = samples, strand = selected_strand)
  if (coordinate == "transcript") {
    dt <- map_signal_to_exons(dt, ex)
    query_start <- 1L
    query_end <- as.integer(sum(as.integer(ex[["exon_end"]]) - as.integer(ex[["exon_start"]]) + 1L, na.rm = TRUE))
  } else {
    query_start <- as.integer(tx$tx_start[1])
    query_end <- as.integer(tx$tx_end[1])
  }
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

  if (plot_type == "frame_bar") {
    if (coordinate != "genomic") {
      warning("`plot_type = 'frame_bar'` uses genomic CDS coordinates. `coordinate` was reset to 'genomic'.", call. = FALSE)
      dt <- retrieve_bwg(signal, tx$chrom[1], tx$tx_start[1], tx$tx_end[1], samples = samples, strand = selected_strand)
      if (!is.null(bin_size)) dt <- bin_bwg(dt, bin_size = bin_size)
      dt <- complete_empty_signal_tracks(
        dt,
        sample_ids = expected_samples,
        chrom = tx$chrom[1],
        start = as.integer(tx$tx_start[1]),
        end = as.integer(tx$tx_end[1]),
        strand = selected_strand
      )
      query_start <- as.integer(tx$tx_start[1])
      query_end <- as.integer(tx$tx_end[1])
    }
    frame_dt <- build_cds_frame_table(annotation, transcript_id = transcript_id_value)
    dt <- map_signal_to_cds_frame(dt, frame_dt)
  }

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
    signal_colors = signal_colors,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    text_color = text_color,
    text_size = text_size,
    grid_linewidth = grid_linewidth,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary,
    frame_palette = frame_palette,
    frame_colors = frame_colors,
    show_frame_legend = show_frame_legend
  )

  if (!show_gene_model) return(p_signal)
  p_model <- plot_transcript(
    annotation,
    transcript_id = transcript_id,
    coordinate = coordinate,
    show_cds = TRUE,
    cds_width = cds_width,
    utr_width = utr_width,
    show_direction = show_direction,
    direction_mode = direction_mode,
    highlight = highlight,
    label_position = label_position,
    label_by = label_by,
    label_side = label_side,
    label_offset = label_offset,
    text_color = text_color,
    text_size = text_size,
    label_size = label_size
  )
  p_signal / p_model + patchwork::plot_layout(heights = c(3, 1))
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
#' @param plot_type Signal plot type: `bar`, `line`, `area`, `heatmap`, or `frame_bar`. `frame_bar` colors per-base CDS signal by CDS frame and is useful for Ribo-seq tracks.
#' @param strand Strand selector. Use `auto` to use the gene strand.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates.
#' @param show_gene_model Whether to append the gene model track.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. If the number of samples or groups exceeds the palette maximum, colors are automatically interpolated.
#' @param signal_colors Optional named or unnamed vector of colors for samples. If supplied, it overrides `signal_palette`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param grid_linewidth Grid line width in the signal panel.
#' @param cds_width Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_width Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param show_direction Whether to draw direction arrows in the gene model track.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis, `feature` draws labels at the center of each gene/transcript structure, and `none` hides labels.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_color Text color for signal and gene model text.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @param label_size Text size in points for gene model labels drawn when `label_position = "feature"`.
#' @param label_side Label placement when `label_position = "feature"`. Use `above`, `below`, or `center`.
#' @param label_offset Vertical offset used for feature labels when `label_side` is `above` or `below`.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, and `end` draws one short arrow at the directional end of each gene.
#' @details
#' `samples` selects the samples to draw. `sample_groups` can be used for
#' group-level coloring and replicate summaries. When raw intervals differ among
#' samples, set `bin_size` before using `signal_summary` so the summary is made
#' on comparable bins. `signal_transform` changes the plotted y value only; the
#' queried signal table is not modified. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot or patchwork object.
#' @examples
#' \dontrun{
#' groups <- c(KO1 = "KO", KO2 = "KO", WT1 = "WT", WT2 = "WT")
#' plot_signal_gene(
#'   signal = bg,
#'   annotation = gp,
#'   gene_id = "GeneA",
#'   samples = c("KO1", "KO2", "WT1", "WT2"),
#'   sample_groups = groups,
#'   signal_color_by = "group",
#'   signal_summary = "mean",
#'   bin_size = 50,
#'   signal_palette = "RdBu",
#'   signal_transform = "log2",
#'   signal_y_scale = "fixed",
#'   signal_y_ticks = "range",
#'   label_position = "feature"
#' )
#' }
#' @export
plot_signal_gene <- function(signal, annotation, gene_id, samples = NULL, sample_groups = NULL, signal_color_by = c("sample", "group"), signal_summary = c("none", "mean", "median", "sum"), plot_type = c("bar", "line", "area", "heatmap", "frame_bar"), strand = c("auto", "+", "-", "both", "ignore"), bin_size = NULL, highlight = NULL, show_gene_model = TRUE, signal_palette = "Blues", signal_colors = NULL, signal_transform = c("none", "log2", "log10", "sqrt"), signal_y_scale = c("free", "fixed"), signal_y_ticks = c("range", "pretty"), grid_linewidth = 0.25, heatmap_bin_size = NULL, heatmap_max_bins = 800L, heatmap_summary = c("mean", "max", "sum", "median"), frame_palette = "Set1", frame_colors = NULL, show_frame_legend = TRUE, cds_width = 0.50, utr_width = 0.25, show_direction = TRUE, direction_mode = c("transcript", "gene", "end"), label_position = c("axis", "feature", "none"), label_by = c("gene", "transcript"), label_side = c("above", "below", "center"), label_offset = 0.45, text_color = "black", text_size = 14, label_size = 12) {
  stop_if_not(inherits(signal, "BwgTrack"), "`signal` must be a BwgTrack object.")
  stop_if_not(inherits(annotation, "GenePred"), "`annotation` must be a GenePred object.")
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  label_side <- match.arg(label_side)

  gene_id_value <- as.character(gene_id)
  tx_all <- annotation$transcripts
  tx <- tx_all[tx_all[["gene_id"]] == gene_id_value]
  stop_if_not(nrow(tx) > 0L, "Gene ID was not found.")
  gene <- build_gene_table(tx)
  selected_strand <- if (strand == "auto") gene$strand[1] else strand
  expected_samples <- get_expected_signal_samples(signal, samples = samples, strand = selected_strand)

  dt <- retrieve_bwg(signal, gene$chrom[1], gene$gene_start[1], gene$gene_end[1], samples = samples, strand = selected_strand)
  if (!is.null(bin_size)) dt <- bin_bwg(dt, bin_size = bin_size)
  dt <- complete_empty_signal_tracks(
    dt,
    sample_ids = expected_samples,
    chrom = gene$chrom[1],
    start = gene$gene_start[1],
    end = gene$gene_end[1],
    strand = selected_strand
  )
  if (plot_type == "frame_bar") {
    frame_dt <- build_cds_frame_table(annotation, gene_id = gene_id_value)
    dt <- map_signal_to_cds_frame(dt, frame_dt)
  }
  p_signal <- plot_signal_core(
    dt,
    plot_type = plot_type,
    highlight = highlight,
    x_label = paste0("Chromosome ", as.character(gene$chrom[1]), " position (bp)"),
    signal_palette = signal_palette,
    signal_colors = signal_colors,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    text_color = text_color,
    text_size = text_size,
    grid_linewidth = grid_linewidth,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary,
    frame_palette = frame_palette,
    frame_colors = frame_colors,
    show_frame_legend = show_frame_legend
  )

  if (!show_gene_model) return(p_signal)
  p_model <- plot_gene(
    annotation,
    gene_id = gene_id,
    collapse = "none",
    coordinate = "genomic",
    highlight = highlight,
    cds_width = cds_width,
    utr_width = utr_width,
    show_direction = show_direction,
    direction_mode = direction_mode,
    label_position = label_position,
    label_by = label_by,
    label_side = label_side,
    label_offset = label_offset,
    text_color = text_color,
    text_size = text_size,
    label_size = label_size
  )
  p_signal / p_model + patchwork::plot_layout(heights = c(3, 1))
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
#' @param plot_type Signal plot type: `bar`, `line`, `area`, `heatmap`, or `frame_bar`. `frame_bar` colors per-base CDS signal by CDS frame and is useful for Ribo-seq tracks.
#' @param strand Strand selector.
#' @param bin_size Optional bin size for signal aggregation.
#' @param highlight Optional data frame used to shade intervals on the signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates.
#' @param annotation Optional GenePred object.
#' @param show_gene_model Whether to append a gene model track. Default TRUE.
#' @param signal_palette Signal color palette. Any palette name from `RColorBrewer::brewer.pal.info` can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`, `Dark2`, `YlGnBu`, or `Spectral`. If the number of samples or groups exceeds the palette maximum, colors are automatically interpolated.
#' @param signal_colors Optional named or unnamed vector of colors for samples. If supplied, it overrides `signal_palette`.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`. Log transforms use signed log1p-style transformation to tolerate zero values.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for each sample to have its own y-axis range, or `fixed` to force all samples to share the same y-axis range.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer axis limits as the minimum and maximum ticks, or `pretty` to use ggplot2 default-style breaks.
#' @param grid_linewidth Grid line width in the signal panel.
#' @param cds_width Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_width Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param show_direction Whether to draw direction arrows in the gene model track.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis, `feature` draws labels at the center of each gene/transcript structure, and `none` hides labels.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param text_color Text color for signal and gene model text.
#' @param text_size Text size in points for signal and gene model axis text, axis titles, facet strips, and legends.
#' @param label_size Text size in points for gene model labels drawn when `label_position = "feature"`.
#' @param label_side Label placement when `label_position = "feature"`. Use `above`, `below`, or `center`.
#' @param label_offset Vertical offset used for feature labels when `label_side` is `above` or `below`.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, and `end` draws one short arrow at the directional end of each gene.
#' @details
#' If `annotation` is supplied and `show_gene_model = TRUE`, a gene model track
#' is appended below the signal panel. `signal_colors` can be an unnamed vector
#' or a named vector. Named values are matched to sample IDs when
#' `signal_color_by = "sample"`, and to group names when
#' `signal_color_by = "group"`. See also [GeneTrackR-advanced-parameters].
#' @return A ggplot or patchwork object.
#' @examples
#' \dontrun{
#' plot_signal_region(
#'   signal = bg,
#'   annotation = gp,
#'   chrom = "chr1",
#'   start = 1,
#'   end = 1000,
#'   samples = c("sampleA", "sampleB"),
#'   signal_colors = c(sampleA = "#2166AC", sampleB = "#B2182B"),
#'   plot_type = "bar",
#'   signal_y_scale = "free",
#'   signal_y_ticks = "range",
#'   show_gene_model = TRUE,
#'   label_by = "gene"
#' )
#' }
#' @export
plot_signal_region <- function(signal, chrom, start, end, samples = NULL, sample_groups = NULL, signal_color_by = c("sample", "group"), signal_summary = c("none", "mean", "median", "sum"), plot_type = c("bar", "line", "area", "heatmap", "frame_bar"), strand = c("ignore", "+", "-", "both"), bin_size = NULL, highlight = NULL, annotation = NULL, show_gene_model = TRUE, signal_palette = "Blues", signal_colors = NULL, signal_transform = c("none", "log2", "log10", "sqrt"), signal_y_scale = c("free", "fixed"), signal_y_ticks = c("range", "pretty"), grid_linewidth = 0.25, heatmap_bin_size = NULL, heatmap_max_bins = 800L, heatmap_summary = c("mean", "max", "sum", "median"), frame_palette = "Set1", frame_colors = NULL, show_frame_legend = TRUE, cds_width = 0.50, utr_width = 0.25, show_direction = TRUE, direction_mode = c("transcript", "gene", "end"), label_position = c("axis", "feature", "none"), label_by = c("gene", "transcript"), label_side = c("above", "below", "center"), label_offset = 0.45, text_color = "black", text_size = 14, label_size = 12) {
  stop_if_not(inherits(signal, "BwgTrack"), "`signal` must be a BwgTrack object.")
  plot_type <- match.arg(plot_type)
  strand <- match.arg(strand)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  label_side <- match.arg(label_side)

  expected_samples <- get_expected_signal_samples(signal, samples = samples, strand = strand)
  dt <- retrieve_bwg(signal, chrom, start, end, samples = samples, strand = strand)
  if (!is.null(bin_size)) dt <- bin_bwg(dt, bin_size = bin_size)
  dt <- complete_empty_signal_tracks(
    dt,
    sample_ids = expected_samples,
    chrom = chrom,
    start = start,
    end = end,
    strand = strand
  )
  if (plot_type == "frame_bar") {
    stop_if_not(!is.null(annotation), "`annotation` is required for `plot_type = 'frame_bar'` in region-level plotting.")
    frame_dt <- build_cds_frame_table(annotation, chrom = chrom, start = start, end = end)
    dt <- map_signal_to_cds_frame(dt, frame_dt)
  }
  p_signal <- plot_signal_core(
    dt,
    plot_type = plot_type,
    highlight = highlight,
    x_label = paste0("Chromosome ", as.character(chrom), " position (bp)"),
    signal_palette = signal_palette,
    signal_colors = signal_colors,
    sample_groups = sample_groups,
    signal_color_by = signal_color_by,
    signal_summary = signal_summary,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    text_color = text_color,
    text_size = text_size,
    grid_linewidth = grid_linewidth,
    heatmap_bin_size = heatmap_bin_size,
    heatmap_max_bins = heatmap_max_bins,
    heatmap_summary = heatmap_summary,
    frame_palette = frame_palette,
    frame_colors = frame_colors,
    show_frame_legend = show_frame_legend
  ) +
    ggplot2::coord_cartesian(xlim = c(start, end))

  if (is.null(annotation) || !show_gene_model) return(p_signal)
  p_model <- plot_region(
    annotation,
    chrom,
    start,
    end,
    mode = "overlap",
    collapse = "none",
    highlight = highlight,
    cds_width = cds_width,
    utr_width = utr_width,
    show_direction = show_direction,
    direction_mode = direction_mode,
    label_position = label_position,
    label_by = label_by,
    label_side = label_side,
    label_offset = label_offset,
    text_color = text_color,
    text_size = text_size,
    label_size = label_size
  )
  p_signal / p_model + patchwork::plot_layout(heights = c(3, 1))
}


get_expected_signal_samples <- function(signal, samples = NULL, strand = "ignore", strand_policy = "ignore_unstranded") {
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

complete_empty_signal_tracks <- function(dt, sample_ids, chrom, start, end, strand = "*") {
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
    dt[, "sample_id" := factor(as.character(dt[["sample_id"]]), levels = sample_ids)]
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
  out[, "sample_id" := factor(as.character(out[["sample_id"]]), levels = sample_ids)]
  out[]
}


aggregate_signal_for_heatmap <- function(dt, heatmap_bin_size = NULL, heatmap_max_bins = 800L, heatmap_summary = c("mean", "max", "sum", "median")) {
  heatmap_summary <- match.arg(heatmap_summary)
  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L) return(dt[])

  x_min <- suppressWarnings(min(as.integer(dt[["start"]]), na.rm = TRUE))
  x_max <- suppressWarnings(max(as.integer(dt[["end"]]), na.rm = TRUE))
  if (!is.finite(x_min) || !is.finite(x_max) || x_max < x_min) return(dt[])

  region_width <- as.integer(x_max - x_min + 1L)
  if (is.null(heatmap_bin_size)) {
    heatmap_max_bins <- suppressWarnings(as.integer(heatmap_max_bins)[1L])
    if (!is.finite(heatmap_max_bins) || heatmap_max_bins <= 0L) heatmap_max_bins <- 800L
    interval_count <- length(unique(paste(dt[["sample_id"]], dt[["start"]], dt[["end"]], sep = "\r")))
    if (region_width <= heatmap_max_bins && interval_count <= heatmap_max_bins * length(unique(dt[["sample_id"]]))) {
      return(dt[])
    }
    heatmap_bin_size <- max(1L, as.integer(ceiling(region_width / heatmap_max_bins)))
  } else {
    heatmap_bin_size <- suppressWarnings(as.integer(heatmap_bin_size)[1L])
    if (!is.finite(heatmap_bin_size) || heatmap_bin_size <= 0L) {
      return(dt[])
    }
  }

  dt[, heatmap_bin := as.integer(floor((mid - x_min) / heatmap_bin_size))]
  dt[heatmap_bin < 0L, heatmap_bin := 0L]

  group_cols <- c("sample_id", "heatmap_bin")
  if ("sample_group" %in% names(dt)) group_cols <- c("sample_group", group_cols)
  if ("chrom" %in% names(dt)) group_cols <- c("chrom", group_cols)
  if ("strand" %in% names(dt)) group_cols <- c("strand", group_cols)

  value_fun <- switch(
    heatmap_summary,
    mean = function(x) mean(x, na.rm = TRUE),
    max = function(x) max(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE)
  )

  out <- dt[, .(
    start = x_min + min(heatmap_bin) * heatmap_bin_size,
    end = pmin(x_max, x_min + (max(heatmap_bin) + 1L) * heatmap_bin_size - 1L),
    value = value_fun(value),
    plot_value = value_fun(plot_value)
  ), by = group_cols]
  out[, mid := (start + end) / 2]
  out[, heatmap_bin := NULL]
  out[]
}

plot_signal_core <- function(dt, plot_type = c("bar", "line", "area", "heatmap", "frame_bar"), highlight = NULL, x_label = "Genomic coordinate", signal_palette = "Blues", signal_colors = NULL, sample_groups = NULL, signal_color_by = c("sample", "group"), signal_summary = c("none", "mean", "median", "sum"), signal_transform = c("none", "log2", "log10", "sqrt"), signal_y_scale = c("free", "fixed"), signal_y_ticks = c("range", "pretty"), text_color = "black", text_size = 14, grid_linewidth = 0.25, heatmap_bin_size = NULL, heatmap_max_bins = 800L, heatmap_summary = c("mean", "max", "sum", "median"), frame_palette = "Set1", frame_colors = NULL, show_frame_legend = TRUE) {
  plot_type <- match.arg(plot_type)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  signal_transform <- match.arg(signal_transform)
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  heatmap_summary <- match.arg(heatmap_summary)
  dt <- data.table::as.data.table(dt)
  stop_if_not(nrow(dt) > 0L, "No signal records were found in the specified region.")
  dt[, mid := (start + end) / 2]
  dt <- apply_signal_grouping(
    dt,
    sample_groups = sample_groups,
    signal_summary = signal_summary
  )
  dt[, mid := (start + end) / 2]
  dt[, plot_value := transform_signal_value(value, signal_transform)]
  sample_ids <- get_ordered_signal_ids(dt, "sample_id")
  color_ids <- if (signal_color_by == "group" && "sample_group" %in% names(dt)) {
    get_ordered_signal_ids(dt, "sample_group")
  } else {
    sample_ids
  }
  dt[, "sample_id" := factor(as.character(dt[["sample_id"]]), levels = sample_ids)]
  if ("sample_group" %in% names(dt)) {
    group_levels <- get_ordered_signal_ids(dt, "sample_group")
    dt[, "sample_group" := factor(as.character(dt[["sample_group"]]), levels = group_levels)]
  }
  discrete_signal_colors <- normalize_signal_colors(
    sample_ids = color_ids,
    signal_palette = signal_palette,
    signal_colors = signal_colors
  )
  color_aes <- if (signal_color_by == "group" && "sample_group" %in% names(dt)) "sample_group" else "sample_id"

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
    signal_y_ticks = signal_y_ticks
  )
  y_anchor_dt <- make_signal_y_anchor(dt, signal_y_scale = signal_y_scale, signal_y_ticks = signal_y_ticks)

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.text.y = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.x = ggplot2::element_text(color = text_color, size = text_size),
      axis.title.y = ggplot2::element_text(color = text_color, size = text_size),
      legend.text = ggplot2::element_text(color = text_color, size = text_size),
      legend.title = ggplot2::element_text(color = text_color, size = text_size),
      strip.text = ggplot2::element_text(color = text_color, size = text_size),
      strip.text.y = ggplot2::element_text(color = text_color, size = text_size, angle = 0),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = grid_linewidth),
      panel.grid.minor = ggplot2::element_line(linewidth = grid_linewidth)
    )

  if (plot_type == "heatmap") {
    dt <- aggregate_signal_for_heatmap(
      dt,
      heatmap_bin_size = heatmap_bin_size,
      heatmap_max_bins = heatmap_max_bins,
      heatmap_summary = heatmap_summary
    )
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = mid, y = sample_id, fill = plot_value, width = pmax(end - start + 1L, 1L))) +
      ggplot2::geom_tile(height = 0.90) +
      ggplot2::labs(x = x_label, y = NULL, fill = y_label) +
      apply_signal_continuous_fill_scale(signal_palette = signal_palette, signal_colors = signal_colors) +
      base_theme +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    return(add_highlight_layer(p, highlight))
  }

  if (plot_type == "line") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = mid, y = plot_value, color = .data[[color_aes]], group = sample_id)) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, color = "Sample") +
      ggplot2::scale_color_manual(values = discrete_signal_colors) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p + ggplot2::geom_blank(
        data = y_anchor_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::geom_line(linewidth = 0.4)
  } else if (plot_type == "bar") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = mid, y = plot_value, fill = .data[[color_aes]])) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, fill = "Sample") +
      ggplot2::scale_fill_manual(values = discrete_signal_colors) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p + ggplot2::geom_blank(
        data = y_anchor_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::geom_col(width = pmax(dt$end - dt$start + 1L, 1), alpha = 0.85)
  } else if (plot_type == "frame_bar") {
    stop_if_not("cds_frame" %in% names(dt), "`cds_frame` was not found. Use `plot_type = 'frame_bar'` only with CDS frame-mapped signal.")
    frame_levels <- c("frame0", "frame1", "frame2")
    dt[, "cds_frame" := factor(as.character(cds_frame), levels = frame_levels)]
    frame_values <- normalize_signal_colors(
      sample_ids = frame_levels,
      signal_palette = frame_palette,
      signal_colors = frame_colors
    )
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = mid, y = plot_value, fill = .data$cds_frame)) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, fill = "CDS frame") +
      ggplot2::scale_fill_manual(values = frame_values, drop = FALSE) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p + ggplot2::geom_blank(
        data = y_anchor_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::geom_col(width = pmax(dt$end - dt$start + 1L, 1), alpha = 0.88)
    if (!isTRUE(show_frame_legend)) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
  } else if (plot_type == "area") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = mid, y = plot_value, fill = .data[[color_aes]], group = sample_id)) +
      ggplot2::facet_grid(sample_id ~ ., scales = facet_scales) +
      ggplot2::labs(x = x_label, y = y_label, fill = "Sample") +
      ggplot2::scale_fill_manual(values = discrete_signal_colors) +
      y_scale +
      base_theme
    if (!is.null(y_anchor_dt)) {
      p <- p + ggplot2::geom_blank(
        data = y_anchor_dt,
        ggplot2::aes(x = .data$mid, y = .data$plot_value),
        inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::geom_area(alpha = 0.75)
  }

  add_highlight_layer(p, highlight)
}


build_cds_frame_table <- function(annotation,
                                  gene_id = NULL,
                                  transcript_id = NULL,
                                  chrom = NULL,
                                  start = NULL,
                                  end = NULL) {
  stop_if_not(is_gene_model_feature(annotation), "`annotation` must contain transcript and exon tables for frame-aware plotting.")
  annotation <- as_genepred(annotation)
  tx <- data.table::as.data.table(annotation$transcripts)
  ex <- data.table::as.data.table(annotation$exons)
  stop_if_not(nrow(tx) > 0L && nrow(ex) > 0L, "No transcript/exon records are available for frame-aware plotting.")

  if (!is.null(transcript_id)) {
    query_transcript_id <- as.character(transcript_id)[1L]
    tx <- tx[as.character(tx[["transcript_id"]]) == query_transcript_id]
  } else if (!is.null(gene_id)) {
    query_gene_id <- as.character(gene_id)[1L]
    tx <- tx[as.character(tx[["gene_id"]]) == query_gene_id]
  } else {
    stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end), "Specify `gene_id`, `transcript_id`, or `chrom/start/end` for frame-aware plotting.")
    query_chrom <- as.character(chrom)[1L]
    query_start <- as.integer(start)[1L]
    query_end <- as.integer(end)[1L]
    tx <- tx[
      as.character(tx[["chrom"]]) == query_chrom &
        as.integer(tx[["tx_start"]]) <= query_end &
        as.integer(tx[["tx_end"]]) >= query_start
    ]
  }

  stop_if_not(nrow(tx) > 0L, "No transcript was found for frame-aware plotting.")
  tx <- tx[!is.na(cds_start) & !is.na(cds_end) & as.integer(cds_start) <= as.integer(cds_end)]
  stop_if_not(nrow(tx) > 0L, "No coding transcript was found for frame-aware plotting.")

  # Gene-level queries can contain multiple isoforms. Use the longest CDS
  # transcript as the frame reference to avoid conflicting frame assignments.
  tx[, "cds_len_ref" := as.integer(cds_end) - as.integer(cds_start) + 1L]
  data.table::setorderv(tx, c("cds_len_ref", "transcript_id"), order = c(-1L, 1L))
  tx <- tx[1L]

  ex <- ex[as.character(ex[["transcript_id"]]) == as.character(tx$transcript_id[1L])]
  stop_if_not(nrow(ex) > 0L, "No exons were found for the selected coding transcript.")

  cds_start <- as.integer(tx$cds_start[1L])
  cds_end <- as.integer(tx$cds_end[1L])
  strand <- as.character(tx$strand[1L])
  tx_gene_id <- as.character(tx$gene_id[1L])
  tx_transcript_id <- as.character(tx$transcript_id[1L])
  cds <- ex[, .(
    chrom = as.character(chrom),
    start = pmax(as.integer(exon_start), cds_start),
    end = pmin(as.integer(exon_end), cds_end),
    gene_id = tx_gene_id,
    transcript_id = tx_transcript_id,
    strand = strand
  )]
  cds <- cds[start <= end]
  stop_if_not(nrow(cds) > 0L, "No CDS-overlapping exon segment was found for frame-aware plotting.")

  if (identical(strand, "-")) {
    data.table::setorderv(cds, c("end", "start"), order = c(-1L, -1L))
  } else {
    data.table::setorderv(cds, c("start", "end"), order = c(1L, 1L))
  }

  pieces <- vector("list", nrow(cds))
  offset <- 0L
  for (i in seq_len(nrow(cds))) {
    if (identical(strand, "-")) {
      pos <- seq.int(cds$end[i], cds$start[i], by = -1L)
    } else {
      pos <- seq.int(cds$start[i], cds$end[i], by = 1L)
    }
    n <- length(pos)
    cds_offset <- offset + seq.int(0L, n - 1L)
    pieces[[i]] <- data.table::data.table(
      chrom = cds$chrom[i],
      pos = as.integer(pos),
      gene_id = cds$gene_id[i],
      transcript_id = cds$transcript_id[i],
      strand = strand,
      cds_pos = as.integer(cds_offset + 1L),
      cds_frame = paste0("frame", cds_offset %% 3L)
    )
    offset <- offset + n
  }

  out <- data.table::rbindlist(pieces)
  data.table::setorderv(out, c("chrom", "pos"))
  out[]
}

map_signal_to_cds_frame <- function(signal_dt, frame_dt, max_expand = 2000000L) {
  dt <- data.table::as.data.table(signal_dt)
  frame_dt <- data.table::as.data.table(frame_dt)
  stop_if_not(nrow(frame_dt) > 0L, "No CDS frame records were available for frame-aware plotting.")
  if (nrow(dt) == 0L) return(dt[0])

  dt[, "start" := as.integer(start)]
  dt[, "end" := as.integer(end)]
  dt[is.na(end) | end < start, "end" := start]
  widths <- as.integer(dt$end - dt$start + 1L)
  total_width <- sum(widths, na.rm = TRUE)
  stop_if_not(total_width <= as.integer(max_expand), paste0(
    "The signal region is too large to expand for `frame_bar` (", total_width,
    " positions). Use a smaller gene/region or pre-filter the signal."
  ))

  idx <- rep(seq_len(nrow(dt)), widths)
  pos <- unlist(mapply(seq.int, dt$start, dt$end, SIMPLIFY = FALSE), use.names = FALSE)
  expanded <- dt[idx]
  expanded[, `:=`(start = as.integer(pos), end = as.integer(pos), mid = as.numeric(pos))]

  frame_key <- frame_dt[, .(chrom = as.character(chrom), start = as.integer(pos), cds_pos, cds_frame)]
  expanded[, "chrom" := as.character(chrom)]
  mapped <- merge(expanded, frame_key, by = c("chrom", "start"), all.x = FALSE, sort = FALSE)
  stop_if_not(nrow(mapped) > 0L, "No signal records overlapped CDS frame positions.")
  mapped[, "end" := as.integer(start)]
  mapped[, "mid" := as.numeric(start)]
  mapped[]
}

make_signal_y_scale <- function(values = NULL,
                                signal_y_scale = c("free", "fixed"),
                                signal_y_ticks = c("range", "pretty")) {
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  if (signal_y_ticks == "pretty") {
    return(ggplot2::scale_y_continuous())
  }

  make_integer_limits <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NULL)
    }
    lim <- c(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)))
    if (identical(lim[1L], lim[2L])) {
      if (lim[1L] == 0) {
        lim <- c(0, 1)
      } else if (lim[1L] > 0) {
        lim <- c(0, lim[2L])
      } else {
        lim <- c(lim[1L], 0)
      }
    }
    lim
  }

  integer_label <- function(x) {
    formatC(as.integer(round(x)), format = "d", big.mark = ",")
  }

  if (signal_y_scale == "fixed") {
    lim <- make_integer_limits(values)
    if (is.null(lim)) {
      return(ggplot2::scale_y_continuous())
    }
    return(
      ggplot2::scale_y_continuous(
        limits = lim,
        breaks = lim,
        labels = integer_label,
        expand = ggplot2::expansion(mult = c(0, 0.03))
      )
    )
  }

  ggplot2::scale_y_continuous(
    breaks = function(x) {
      lim <- make_integer_breaks_from_panel_limits(x)
      if (is.null(lim)) {
        return(NULL)
      }
      lim
    },
    labels = integer_label,
    expand = ggplot2::expansion(mult = c(0, 0.03))
  )
}

make_integer_breaks_from_panel_limits <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NULL)
  }
  lower_raw <- min(x, na.rm = TRUE)
  upper_raw <- max(x, na.rm = TRUE)
  lower <- if (lower_raw < 0) ceiling(lower_raw) else floor(lower_raw)
  upper <- if (upper_raw < 0) ceiling(upper_raw) else floor(upper_raw)
  if (!is.finite(lower) || !is.finite(upper)) {
    return(NULL)
  }
  if (lower == upper) {
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

make_signal_y_anchor <- function(dt, signal_y_scale = c("free", "fixed"), signal_y_ticks = c("range", "pretty")) {
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  if (signal_y_scale != "free" || signal_y_ticks != "range") {
    return(NULL)
  }
  dt <- data.table::as.data.table(dt)
  if (nrow(dt) == 0L) {
    return(NULL)
  }
  anchor <- dt[
    , {
      lim <- c(floor(min(as.numeric(plot_value), na.rm = TRUE)), ceiling(max(as.numeric(plot_value), na.rm = TRUE)))
      if (identical(lim[1L], lim[2L])) {
        if (lim[1L] == 0) {
          lim <- c(0, 1)
        } else if (lim[1L] > 0) {
          lim <- c(0, lim[2L])
        } else {
          lim <- c(lim[1L], 0)
        }
      }
      .(plot_value = as.numeric(lim))
    },
    by = "sample_id"
  ]
  mid_value <- stats::median(as.numeric(dt[["mid"]]), na.rm = TRUE)
  if (!is.finite(mid_value)) {
    mid_value <- 0
  }
  anchor[, "mid" := mid_value]
  anchor
}

transform_signal_value <- function(value, signal_transform = c("none", "log2", "log10", "sqrt")) {
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

make_signal_palette <- function(n, signal_palette = "Blues") {
  n <- max(1L, as.integer(n))
  signal_palette <- as.character(signal_palette)[1L]
  if (is.na(signal_palette) || !nzchar(signal_palette)) {
    signal_palette <- "Blues"
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      signal_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal_info <- RColorBrewer::brewer.pal.info[signal_palette, , drop = FALSE]
    max_colors <- as.integer(pal_info[["maxcolors"]][1L])
    base_n <- max(3L, min(max_colors, max(3L, n)))
    base <- RColorBrewer::brewer.pal(base_n, signal_palette)
    if (n > max_colors) {
      base <- grDevices::colorRampPalette(base)(n)
      return(base)
    }
  } else {
    predefined <- list(
      Blues = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
      Reds = c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
      RdBu = c("#B2182B", "#EF8A62", "#FDDDBC", "#D1E5F0", "#67A9CF", "#2166AC")
    )
    base <- predefined[[signal_palette]]
    if (is.null(base)) {
      warning(
        sprintf("Unknown `signal_palette`: %s. Falling back to 'Blues'.", signal_palette),
        call. = FALSE
      )
      base <- predefined[["Blues"]]
    }
  }
  grDevices::colorRampPalette(base)(n)
}

normalize_signal_colors <- function(sample_ids, signal_palette = "Blues", signal_colors = NULL) {
  sample_ids <- unique(as.character(sample_ids))
  n <- length(sample_ids)
  if (!is.null(signal_colors)) {
    cols <- as.character(signal_colors)
    if (!is.null(names(cols)) && all(sample_ids %in% names(cols))) {
      return(cols[sample_ids])
    }
    if (length(cols) < n) {
      cols <- rep(cols, length.out = n)
    }
    names(cols) <- sample_ids
    return(cols)
  }
  cols <- make_signal_palette(n, signal_palette)
  names(cols) <- sample_ids
  cols
}

apply_signal_discrete_fill_scale <- function(sample_ids, signal_palette = "Blues", signal_colors = NULL) {
  cols <- normalize_signal_colors(sample_ids, signal_palette = signal_palette, signal_colors = signal_colors)
  ggplot2::scale_fill_manual(values = cols)
}

apply_signal_discrete_color_scale <- function(sample_ids, signal_palette = "Blues", signal_colors = NULL) {
  cols <- normalize_signal_colors(sample_ids, signal_palette = signal_palette, signal_colors = signal_colors)
  ggplot2::scale_color_manual(values = cols)
}

apply_signal_continuous_fill_scale <- function(signal_palette = "Blues", signal_colors = NULL) {
  if (!is.null(signal_colors)) {
    cols <- as.character(signal_colors)
  } else {
    cols <- make_signal_palette(256L, signal_palette = signal_palette)
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
    stop_if_not(all(c("sample_id", "group") %in% names(sample_groups)), "`sample_groups` data frame must contain `sample_id` and `group` columns.")
    map <- stats::setNames(as.character(sample_groups[["group"]]), as.character(sample_groups[["sample_id"]]))
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
  stop_if_not(length(groups) == length(sample_ids), "Unnamed `sample_groups` must have the same length as selected samples.")
  stats::setNames(groups, sample_ids)
}

apply_signal_grouping <- function(dt, sample_groups = NULL, signal_summary = c("none", "mean", "median", "sum")) {
  signal_summary <- match.arg(signal_summary)
  dt <- data.table::copy(data.table::as.data.table(dt))
  sample_ids <- get_ordered_signal_ids(dt, "sample_id")
  group_map <- normalize_sample_groups(sample_ids, sample_groups)
  dt[, "sample_group" := as.character(group_map[as.character(dt[["sample_id"]])])]
  if (signal_summary == "none") {
    return(dt)
  }
  summary_fun <- switch(
    signal_summary,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
  )
  group_cols <- c("sample_group", "chrom", "start", "end")
  if ("cds_frame" %in% names(dt)) {
    group_cols <- c(group_cols, "cds_frame")
  }
  out <- dt[, .(
    value = summary_fun(value),
    strand = as.character(strand[1L])
  ), by = group_cols]
  group_levels <- unique(as.character(group_map[sample_ids]))
  out[, "sample_id" := as.character(out[["sample_group"]])]
  out[, "sample_id" := factor(as.character(out[["sample_id"]]), levels = group_levels)]
  out[, "sample_group" := factor(as.character(out[["sample_group"]]), levels = group_levels)]
  data.table::setcolorder(out, c("sample_id", "sample_group", "chrom", "start", "end", "value", "strand"))
  out[]
}

map_signal_to_exons <- function(signal_dt, exons) {
  dt <- data.table::as.data.table(signal_dt)
  ex <- data.table::copy(exons)
  if (nrow(dt) == 0L || nrow(ex) == 0L) return(dt[0])

  data.table::setorder(ex, transcript_id, exon_start, exon_end)
  ex[, exon_width := exon_end - exon_start + 1L]
  ex[, exon_offset := cumsum(data.table::shift(exon_width, fill = 0L)), by = transcript_id]
  ex <- ex[, .(chrom, exon_start, exon_end, exon_offset)]

  mapped <- list()
  for (i in seq_len(nrow(ex))) {
    piece <- dt[chrom == ex$chrom[i] & start <= ex$exon_end[i] & end >= ex$exon_start[i]]
    if (nrow(piece) == 0L) next
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
