# Plot signal over a genomic region

Plot signal over a genomic region

## Usage

``` r
plot_signal_region(
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
  cds_height = 0.5,
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
)
```

## Arguments

- signal:

  A BwgTrack object.

- chrom:

  Chromosome name.

- start:

  Region start.

- end:

  Region end.

- samples:

  Optional sample IDs to plot. If NULL, all samples are used.

- sample_groups:

  Optional sample group mapping for group-level coloring or replicate
  summaries. Use a named character vector, a data frame with `sample_id`
  and `group`, or an unnamed vector with one group per selected sample.

- signal_color_by:

  Color signal tracks by `sample` or `group`.

- signal_summary:

  Replicate summary mode. Use `none` to plot individual samples, or
  `mean`, `median`, or `sum` to summarize samples within each group.
  Summary is performed on the current signal intervals, so using
  `bin_size` is recommended when raw interval boundaries differ among
  samples.

- plot_type:

  Signal plot type: `bar`, `line`, or `heatmap`.

- strand:

  Strand selector.

- bin_size:

  Optional bin size for signal aggregation.

- highlight:

  Optional data frame used to shade intervals on the signal and gene
  model tracks. It must contain `start` and `end` columns in genomic
  coordinates.

- annotation:

  Optional GenePred object.

- show_gene_model:

  Whether to append a gene model track. Default TRUE.

- signal_track_height:

  Relative height of the signal panel when the gene model is shown.
  Default 3.

- gene_track_height:

  Relative height of the gene model panel when it is shown. Default 1.

- signal_palette:

  Signal color palette. Any palette name from
  [`RColorBrewer::brewer.pal.info`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  can be used, such as `Blues`, `Reds`, `RdBu`, `Paired`, `Set1`,
  `Dark2`, `YlGnBu`, or `Spectral`. Discrete sample/group colors are
  assigned in the standard RColorBrewer class order; heatmaps use the
  corresponding continuous gradient.

- signal_palette_direction:

  Direction for generated signal colors. Use `1` for the standard
  palette order and `-1` for the reversed palette order. Discrete
  sample/group colors preserve level-to-color order; heatmap gradients
  reverse continuously.

- signal_colors:

  Optional named or unnamed vector of explicit colors for samples. If
  supplied, it overrides `signal_palette`; explicit colors are not
  modified by `signal_palette_direction`.

- signal_transform:

  Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`.
  Log transforms use signed log1p-style transformation to tolerate zero
  values.

- signal_y_scale:

  Signal y-axis scale mode. Use `free` for each sample to have its own
  y-axis range, or `fixed` to force all samples to share the same y-axis
  range.

- signal_y_ticks:

  Signal y-axis tick mode. Use `range` to show only integer axis limits
  as the minimum and maximum ticks, or `pretty` to use ggplot2
  default-style breaks.

- cds_height:

  Vertical thickness of CDS rectangles in the gene model track.

- utr_height:

  Vertical thickness of UTR/non-coding exon rectangles in the gene model
  track.

- direction_mode:

  Direction-arrow style for the gene model track. `transcript` draws one
  arrow per transcript, `gene` draws one arrow per gene, `end` draws one
  short arrow at the directional end of each gene, and `none` hides
  direction arrows.

- label_position:

  Where to draw gene model labels. `axis` draws labels on the y axis and
  `feature` draws labels at the center of each gene/transcript
  structure.

- label_by:

  Which identifier to use for gene model labels. Use `gene` for gene IDs
  or `transcript` for transcript IDs.

- text_size:

  Text size in points for signal and gene model axis text, axis titles,
  facet strips, and legends.

- signal_y_limits:

  Optional two-element numeric vector giving the plotted y-axis limits
  after `signal_transform`. Supplying limits changes `signal_y_scale` to
  `fixed`.

- signal_alpha:

  Signal geometry transparency from 0 to 1.

- signal_bar_width:

  Relative width of bar intervals from greater than 0 to 1. A value
  below 1 creates proportional gaps without changing interval centers.

- plot_theme:

  Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw panel borders. `NULL` preserves the selected theme
  default.

## Value

A ggplot or patchwork object.

## Details

If `annotation` is supplied and `show_gene_model = TRUE`, a gene model
track is appended below the signal panel. `signal_colors` can be an
unnamed vector or a named vector. Named values are matched to sample IDs
when `signal_color_by = "sample"`, and to group names when
`signal_color_by = "group"`. See also
[GeneTrackR-advanced-parameters](https://renscq.github.io/GeneTrackR/reference/GeneTrackR-advanced-parameters.md).

## Examples

``` r
if (FALSE) { # \dontrun{
gp <- read_genepred(
  system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
  format = "genePredExt",
  verbose = FALSE
)
signal_all <- read_bwg(
  system.file(
    "extdata",
    c(
      "gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph",
      "gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"
    ),
    package = "GeneTrackR"
  ),
  format = "bedgraph",
  sample_names = c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus"),
  strand = c("+", "-", "+", "-"),
  mode = "memory"
)
plot_signal_region(
  signal = signal_all,
  annotation = gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  strand = "both",
  plot_type = "bar"
)
} # }
```
