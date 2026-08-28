# Plot variants in a genomic region

Draws VCF-derived variant records in a genomic region. This function is
dedicated to variant visualization and can be used directly or through
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md).

## Usage

``` r
plot_variant(
  variant,
  chrom = NULL,
  start = NULL,
  end = NULL,
  variant_type = NULL,
  color_by = c("variant_type", "filter", "none"),
  label_by = c("none", "variant_id", "variant_type", "ref", "alt"),
  variant_shape = c("lollipop", "point", "rug"),
  variant_palette = "Paired",
  variant_colors = NULL,
  point_size = 2,
  line_width = 0.35,
  text_size = 14,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
)
```

## Arguments

- variant:

  A VariantTrack object.

- chrom:

  Chromosome name. If NULL, all chromosomes in `variant` are used.

- start:

  Region start in 1-based coordinates. If NULL, the minimum variant
  position is used.

- end:

  Region end in 1-based coordinates. If NULL, the maximum variant
  position is used.

- variant_type:

  Optional variant type filter, such as `SNP`, `INS`, `DEL`, or `MNV`.

- color_by:

  Column used for variant colors. Use `none` to disable grouping.

- label_by:

  Column used for labels. Use `none` to hide labels.

- variant_shape:

  Plot geometry. `lollipop` draws vertical stems and points, `point`
  draws points only, and `rug` draws bottom ticks.

- variant_palette:

  RColorBrewer palette name for variant colors.

- variant_colors:

  Optional named or unnamed variant color vector.

- point_size:

  Point size.

- line_width:

  Line width for lollipop/rug stems.

- text_size:

  Text size.

- plot_theme:

  Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw the panel border. `NULL` preserves the selected theme
  default.

## Value

A ggplot object.
