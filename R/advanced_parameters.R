# Author: Rensc
# Date: 2026-05-27
# Version: 0.1.29
# Function: Shared documentation for advanced plotting and query parameters
# Input: None
# Output: R documentation topic

#' Advanced parameter guide for GeneTrackR
#'
#' @description
#' This help topic summarizes frequently used advanced parameters that are shared
#' by several GeneTrackR functions. It is intended as a quick reference for
#' custom colors, highlighted intervals, signal transformations, sample grouping,
#' and strand handling.
#'
#' @section Gene model colors:
#' Gene model functions such as [plot_gene()], [plot_transcript()], and
#' [plot_region()] use `color_palette` and `fill_colors` for exon/CDS/UTR fills.
#' `color_palette` accepts any palette from `RColorBrewer::brewer.pal.info`. If
#' the number of discrete groups exceeds the palette limit, colors are
#' interpolated automatically.
#'
#' `fill_colors` can be either unnamed or named. For gene models, the most stable
#' named form is:
#'
#' ```r
#' fill_colors = c(
#'   CDS = "#33a02c",
#'   UTR = "#b2df8a",
#'   exon = "#fb9a99"
#' )
#' ```
#'
#' Unnamed colors are matched to the observed discrete groups in plotting order.
#'
#' @section Highlight intervals:
#' `highlight` must be a data frame with at least `start` and `end` columns:
#'
#' ```r
#' highlight = data.frame(start = 105000, end = 106000)
#' ```
#'
#' In genomic plots, `start` and `end` are genomic coordinates. In transcript
#' coordinate plots, they are spliced transcript coordinates.
#'
#' @section Signal colors:
#' Signal plotting functions use `signal_palette`, `signal_colors`,
#' `sample_groups`, and `signal_color_by`.
#'
#' ```r
#' signal_palette = "YlGnBu"
#' signal_colors = c(sampleA = "#2166AC", sampleB = "#B2182B")
#' ```
#'
#' Group-level coloring is controlled by a named vector or a two-column data
#' frame:
#'
#' ```r
#' sample_groups = c(KO1 = "KO", KO2 = "KO", WT1 = "WT", WT2 = "WT")
#' signal_color_by = "group"
#' ```
#'
#' @section Signal y-axis:
#' Signal plots support transformed values and per-sample or shared y-axis
#' ranges:
#'
#' ```r
#' signal_transform = "sqrt"  # none, log2, log10, or sqrt
#' signal_y_scale = "free"   # free or fixed
#' signal_y_ticks = "range"  # range or pretty
#' ```
#'
#' `signal_y_ticks = "range"` displays integer y-axis limits only.
#'
#' @section Strand handling:
#' Standard bigWig and wig tracks are unstranded. For such tracks,
#' `strand_policy = "ignore_unstranded"` lets `strand = "+"` or `"-"` return the
#' same unstranded signal. Use `strand_policy = "strict"` when you only want
#' records with explicit strand information, such as plus/minus bedGraph tracks.
#'
#' @seealso [plot_gene()], [plot_transcript()], [plot_region()],
#' [plot_signal_gene()], [plot_signal_transcript()], [plot_signal_region()],
#' [plot_tracks()], [query_bwg()], [read_bwg()]
#' @name GeneTrackR-advanced-parameters
NULL
