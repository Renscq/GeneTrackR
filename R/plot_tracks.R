# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.26
# Function: Combined genome-browser-like track plotting
# Input: GenePred and BwgTrack objects
# Output: Combined patchwork track figure

#' Plot combined signal and gene model tracks
#'
#' @description
#' Draw a genome-browser-like figure containing a gene model track alone or a
#' signal track combined with a gene model track. The genomic interval can be
#' specified in three ways: by `gene_id`, by `transcript_id`, or by explicit
#' `chrom`, `start`, and `end` coordinates.
#'
#' @param annotation A GenePred object.
#' @param signal Optional BwgTrack object. If NULL, only the gene model track is drawn.
#' @param chrom Chromosome name. Required when `gene_id` and `transcript_id` are not supplied.
#' @param start Region start in 1-based closed coordinates. Required with `chrom`/`end`.
#' @param end Region end in 1-based closed coordinates. Required with `chrom`/`start`.
#' @param gene_id Optional gene ID. If supplied, the plotting interval is inferred from the gene locus.
#' @param transcript_id Optional transcript ID. If supplied, the plotting interval is inferred from the transcript locus.
#' @param samples Optional signal sample IDs.
#' @param sample_groups Optional sample group mapping for group-level coloring or replicate summaries. Use a named character vector, a data frame with `sample_id` and `group`, or an unnamed vector with one group per selected sample.
#' @param signal_color_by Color signal tracks by `sample` or `group`.
#' @param signal_summary Replicate summary mode. Use `none` to plot individual samples, or `mean`, `median`, or `sum` to summarize samples within each group.
#' @param signal_type Signal plot type.
#' @param signal_palette Signal color palette for signal tracks. Any palette name from `RColorBrewer::brewer.pal.info` can be used.
#' @param signal_colors Optional named or unnamed vector of colors for signal samples.
#' @param signal_transform Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`.
#' @param signal_y_scale Signal y-axis scale mode. Use `free` for independent sample-specific y-axis ranges or `fixed` for a shared y-axis range across samples.
#' @param signal_y_ticks Signal y-axis tick mode. Use `range` to show only integer minimum and maximum limits or `pretty` for default-style breaks.
#' @param grid_linewidth Grid line width in the signal panel.
#' @param collapse Gene model collapse mode for region-level plotting.
#' @param strand Signal strand selector.
#' @param bin_size Optional signal bin size.
#' @param highlight Optional data frame used to shade intervals on signal and gene model tracks. It must contain `start` and `end` columns in genomic coordinates. Optional columns are allowed but ignored by the default renderer.
#' @param layout Track layout. Use `signal_top` to place signal above gene model, or `gene_top` to place gene model above signal.
#' @param heights Relative panel heights. Must contain `signal` and `gene` names.
#' @param cds_width Vertical thickness of CDS rectangles in the gene model track.
#' @param utr_width Vertical thickness of UTR/non-coding exon rectangles in the gene model track.
#' @param show_direction Whether to draw direction arrows in the gene model track.
#' @param direction_mode Direction-arrow style for the gene model track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, and `end` draws one short arrow at the directional end of each gene.
#' @param label_position Where to draw gene model labels. `axis` draws labels on the y axis, `feature` draws labels near the feature, and `none` hides labels.
#' @param label_by Which identifier to use for gene model labels. Use `gene` for gene IDs or `transcript` for transcript IDs.
#' @param label_side Label placement when `label_position = "feature"`. Use `above`, `below`, or `center`.
#' @param label_offset Vertical offset used for feature labels when `label_side` is `above` or `below`.
#' @param text_color Text color for axis text, axis titles, legends, and in-plot labels.
#' @param text_size Text size in points for axis text, axis titles, legends, and facet labels.
#' @param label_size Text size in points for gene model labels drawn when `label_position = "feature"`.
#'
#' @return A patchwork object or ggplot object.
#'
#' @examples
#' \dontrun{
#' gp <- read_genepred(
#'   system.file("extdata", "example.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt"
#' )
#' bg <- read_bwg(
#'   system.file("extdata", c("example_signal_A.bedgraph", "example_signal_B.bedgraph"), package = "GeneTrackR"),
#'   format = "bedgraph"
#' )
#'
#' plot_tracks(annotation = gp, gene_id = "GeneA")
#' plot_tracks(annotation = gp, signal = bg, gene_id = "GeneA")
#' plot_tracks(annotation = gp, signal = bg, transcript_id = "TxA1")
#' plot_tracks(annotation = gp, signal = bg, chrom = "chr1", start = 1, end = 1200)
#' }
#' @export
plot_tracks <- function(annotation,
                        signal = NULL,
                        chrom = NULL,
                        start = NULL,
                        end = NULL,
                        gene_id = NULL,
                        transcript_id = NULL,
                        samples = NULL,
                        sample_groups = NULL,
                        signal_color_by = c("sample", "group"),
                        signal_summary = c("none", "mean", "median", "sum"),
                        signal_type = c("bar", "line", "area", "heatmap"),
                        signal_palette = "Blues",
                        signal_colors = NULL,
                        signal_transform = c("none", "log2", "log10", "sqrt"),
                        signal_y_scale = c("free", "fixed"),
                        signal_y_ticks = c("range", "pretty"),
                        grid_linewidth = 0.25,
                        collapse = c("none", "union_exon", "longest"),
                        strand = c("ignore", "+", "-", "both"),
                        bin_size = NULL,
                        highlight = NULL,
                        layout = c("signal_top", "gene_top"),
                        heights = c(signal = 3, gene = 1),
                        cds_width = 0.50,
                        utr_width = 0.25,
                        show_direction = TRUE,
                        direction_mode = c("transcript", "gene", "end"),
                        label_position = c("axis", "feature", "none"),
                        label_by = c("gene", "transcript"),
                        label_side = c("above", "below", "center"),
                        label_offset = 0.45,
                        text_color = "black",
                        text_size = 14,
                        label_size = 12) {
  stop_if_not(inherits(annotation, "GenePred"), "`annotation` must be a GenePred object.")
  signal_type <- match.arg(signal_type)
  signal_color_by <- match.arg(signal_color_by)
  signal_summary <- match.arg(signal_summary)
  collapse <- match.arg(collapse)
  strand <- match.arg(strand)
  layout <- match.arg(layout)
  signal_transform <- match.arg(signal_transform)
  signal_y_scale <- match.arg(signal_y_scale)
  signal_y_ticks <- match.arg(signal_y_ticks)
  direction_mode <- match.arg(direction_mode)
  label_position <- match.arg(label_position)
  label_by <- match.arg(label_by)
  label_side <- match.arg(label_side)

  has_gene_id <- !is.null(gene_id)
  has_transcript_id <- !is.null(transcript_id)
  has_region <- !is.null(chrom) || !is.null(start) || !is.null(end)
  n_locator <- sum(c(has_gene_id, has_transcript_id, has_region))
  stop_if_not(n_locator == 1L, "Specify exactly one locator: `gene_id`, `transcript_id`, or `chrom` + `start` + `end`.")

  tx_all <- data.table::as.data.table(annotation$transcripts)

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
  }

  if (is.null(signal)) {
    return(p_gene)
  }
  stop_if_not(inherits(signal, "BwgTrack"), "`signal` must be a BwgTrack object.")
  p_signal <- plot_signal_region(
    signal = signal,
    chrom = chrom_value,
    start = start_value,
    end = end_value,
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
    signal_colors = signal_colors,
    signal_transform = signal_transform,
    signal_y_scale = signal_y_scale,
    signal_y_ticks = signal_y_ticks,
    grid_linewidth = grid_linewidth,
    text_color = text_color,
    text_size = text_size
  )

  if (layout == "signal_top") {
    p_signal / p_gene + patchwork::plot_layout(heights = c(heights["signal"], heights["gene"]))
  } else {
    p_gene / p_signal + patchwork::plot_layout(heights = c(heights["gene"], heights["signal"]))
  }
}
