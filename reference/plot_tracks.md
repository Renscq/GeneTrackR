# Plot combined signal and gene model tracks

Draw a genome-browser-like figure containing a gene model track alone or
a signal track combined with gene, feature, and variant tracks. The
genomic interval can be specified in three ways: by `gene_id`, by
`transcript_id`, or by explicit `chrom`, `start`, and `end` coordinates.

## Usage

``` r
plot_tracks(
  annotation,
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

- annotation:

  A GenePred object or a standardized Feature object with
  transcript/exon records.

- signal:

  Optional BwgTrack object. If NULL, only annotation/feature/variant
  tracks are drawn.

- features:

  Optional FeatureTrack object or named list of FeatureTrack objects
  from
  [`read_bed()`](https://renscq.github.io/GeneTrackR/reference/read_bed.md),
  [`read_gff()`](https://renscq.github.io/GeneTrackR/reference/read_gff.md),
  or
  [`read_gtf()`](https://renscq.github.io/GeneTrackR/reference/read_gtf.md).

- variants:

  Optional VariantTrack object or named list of VariantTrack objects
  from
  [`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md).

- chrom:

  Chromosome name. Required when `gene_id` and `transcript_id` are not
  supplied.

- start:

  Region start in 1-based closed coordinates. Required with
  `chrom`/`end`.

- end:

  Region end in 1-based closed coordinates. Required with
  `chrom`/`start`.

- gene_id:

  Optional gene ID. If supplied, the plotting interval is inferred from
  the gene locus.

- transcript_id:

  Optional transcript ID. If supplied, the plotting interval is inferred
  from the transcript locus.

- samples:

  Optional signal sample IDs.

- sample_groups:

  Optional sample group mapping for group-level coloring or replicate
  summaries. Use a named character vector, a data frame with `sample_id`
  and `group`, or an unnamed vector with one group per selected sample.

- signal_color_by:

  Color signal tracks by `sample` or `group`.

- signal_summary:

  Replicate summary mode. Use `none` to plot individual samples, or
  `mean`, `median`, or `sum` to summarize samples within each group.

- signal_type:

  Signal plot type for non-Ribo-seq browser signal tracks.

- signal_palette:

  Signal color palette for signal tracks. Any palette name from
  [`RColorBrewer::brewer.pal.info`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  can be used. Discrete signal sample/group colors follow the standard
  RColorBrewer class order; heatmaps use the corresponding continuous
  gradient.

- signal_palette_direction:

  Direction for generated signal colors. Use `1` for the standard
  palette order and `-1` for the reversed palette order. Discrete
  sample/group colors preserve level-to-color order; heatmap gradients
  reverse continuously.

- signal_colors:

  Optional named or unnamed vector of explicit colors for signal
  samples. Explicit colors override `signal_palette` and are not
  modified by `signal_palette_direction`.

- gene_palette:

  RColorBrewer palette name used for gene model feature fills.

- gene_colors:

  Optional custom fill colors for gene model features. Use a named
  vector such as
  `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")`.

- gene_border_color:

  Optional rectangle border color for gene model features. Use `NA` to
  hide borders.

- feature_color_by:

  Feature-track color grouping. `auto` chooses a compact grouping
  automatically; other common choices are `feature_group`, `type`,
  `source`, `name`, and `strand`.

- feature_palette:

  RColorBrewer palette name used for feature-track fills.

- feature_colors:

  Optional explicit feature-track fill colors.

- feature_border_color:

  Optional rectangle border color for feature tracks. Use `NA` to hide
  borders.

- feature_max_legend_levels:

  Maximum number of legend groups shown for each feature track.

- variant_palette:

  RColorBrewer palette name used for variant-track colors.

- variant_colors:

  Optional named or unnamed explicit colors passed to
  [`plot_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_variant.md).
  Named vectors are recommended when stable SNP/INS/DEL colors are
  required.

- ribo_signal_type:

  How Ribo-seq samples should be displayed in transcript-centered
  browser plots. `auto` uses `frame` for samples whose IDs look like
  Ribo-seq/RPF tracks when `transcript_id` is supplied, otherwise it
  falls back to standard genomic signal tracks. `bar` always uses the
  standard signal geometry, and `frame` forces frame rendering when
  transcript-centered plotting is possible.

- frame_palette:

  RColorBrewer palette name used for Ribo-seq frame plots.

- frame_colors:

  Optional explicit colors for `frame0`, `frame1`, and `frame2`.

- signal_transform:

  Signal-axis transformation. Use `none`, `log2`, `log10`, or `sqrt`.

- signal_y_scale:

  Signal y-axis scale mode. Use `free` for independent sample-specific
  y-axis ranges or `fixed` for a shared y-axis range across samples.

- signal_y_ticks:

  Signal y-axis tick mode. Use `range` to show only integer minimum and
  maximum limits or `pretty` for default-style breaks.

- collapse:

  Gene model collapse mode for region-level plotting.

- strand:

  Signal strand selector.

- bin_size:

  Optional signal bin size.

- highlight:

  Optional data frame used to shade intervals on signal and gene model
  tracks. It must contain `start` and `end` columns in genomic
  coordinates. Optional columns are allowed but ignored by the default
  renderer.

- layout:

  Track layout. Use `signal_top` to place signal above gene model, or
  `gene_top` to place gene model above signal.

- heights:

  Relative panel heights. Must contain at least `signal`, `gene`,
  `feature`, and `variant` names when those tracks are used.

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
  `feature` draws labels near the feature.

- label_by:

  Which identifier to use for gene model labels. Use `gene` for gene IDs
  or `transcript` for transcript IDs.

- text_size:

  Text size in points for axis text, axis titles, legends, and facet
  labels.

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

  Base ggplot2 theme used by all standard track panels. Use `bw`,
  `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw panel borders. `NULL` preserves the selected theme
  default.

## Value

A patchwork object or ggplot object.
