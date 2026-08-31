# Author: Rensc
# Date: 2026-09-01
# Version: dev005
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
#' [plot_region()] use `gene_palette` and `gene_colors` for exon/CDS/UTR fills.
#' `gene_palette` accepts any palette from `RColorBrewer::brewer.pal.info`. If
#' the number of discrete groups exceeds the palette limit, colors are
#' interpolated automatically.
#'
#' `gene_colors` can be either unnamed or named. For gene models, the most stable
#' named form is:
#'
#' ```r
#' gene_colors = c(
#'   CDS = "#33a02c",
#'   UTR = "#b2df8a",
#'   exon = "#fb9a99"
#' )
#' ```
#'
#' Unnamed colors are matched to the fixed gene-model levels UTR, CDS, and exon, so colors remain stable even when a plot contains only one feature type.
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
#' Signal plotting functions use `signal_palette`, `signal_palette_direction`, `signal_colors`,
#' `sample_groups`, and `signal_color_by`.
#'
#' ```r
#' signal_palette = "YlGnBu"
#' signal_palette_direction = -1
#' signal_colors = c(sampleA = "#2166AC", sampleB = "#B2182B")
#' ```
#'
#' `signal_palette_direction` reverses the standard discrete palette order. Discrete
#' sample/group/frame colors otherwise follow the RColorBrewer class order exactly,
#' while heatmaps use a continuous gradient. Explicit `signal_colors` are treated
#' as exact user mappings and are not reversed.
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
#' signal_y_scale = "fixed"   # free or fixed
#' signal_y_ticks = "range"   # range or pretty
#' signal_y_limits = c(0, 20) # limits after signal_transform
#' ```
#'
#' `signal_y_ticks = "range"` displays integer y-axis limits only.
#' Supplying `signal_y_limits` changes `signal_y_scale` to `"fixed"`.
#'
#' @section Plot appearance:
#' Standard track panels use `plot_theme`, `show_panel_border`, and black text.
#' Signal geometry can be adjusted with `signal_alpha` and
#' `signal_bar_width`:
#'
#' ```r
#' plot_theme = "classic"      # bw, classic, light, or minimal
#' show_panel_border = FALSE
#' signal_alpha = 0.80
#' signal_bar_width = 0.85
#' signal_track_height = 4
#' gene_track_height = 1
#' ```
#'
#' `signal_track_height` and `gene_track_height` control the relative panel
#' heights in `plot_signal_transcript()`, `plot_signal_gene()`, and
#' `plot_signal_region()`. `plot_tracks()` uses its named `heights` vector for
#' the same purpose.
#'
#' `signal_bar_width` is relative to each BedGraph interval. It changes the
#' visual bar width around the interval center without changing genomic
#' positions.
#'
#' @section Strand handling:
#' Standard bigWig and wig tracks are unstranded. For such tracks,
#' `strand_policy = "ignore_unstranded"` lets `strand = "+"` or `"-"` return the
#' same unstranded signal. Use `strand_policy = "strict"` when you only want
#' records with explicit strand information, such as plus/minus bedGraph tracks.
#'
#' @return No value is returned; this is a documentation-only help topic.
#'
#' @seealso [plot_gene()], [plot_transcript()], [plot_region()],
#' [plot_signal_gene()], [plot_signal_transcript()], [plot_signal_region()],
#' [plot_tracks()], [retrieve_bwg()], [read_bwg()]
#' @name GeneTrackR-advanced-parameters
NULL
