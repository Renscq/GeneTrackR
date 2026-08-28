# Plot a VariantTrack object

Backward-compatible alias of
[`plot_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_variant.md).

## Usage

``` r
plot_variant_track(
  track,
  chrom,
  start,
  end,
  color_by = c("variant_type", "filter", "none"),
  label_by = c("none", "variant_id", "variant_type", "ref", "alt"),
  variant_palette = "Paired",
  variant_colors = NULL,
  text_size = 14,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
)
```

## Arguments

- chrom:

  Chromosome name. If NULL, all chromosomes in `variant` are used.

- start:

  Region start in 1-based coordinates. If NULL, the minimum variant
  position is used.

- end:

  Region end in 1-based coordinates. If NULL, the maximum variant
  position is used.

- color_by:

  Column used for variant colors. Use `none` to disable grouping.

- label_by:

  Column used for labels. Use `none` to hide labels.

- variant_palette:

  RColorBrewer palette name for variant colors.

- variant_colors:

  Optional named or unnamed variant color vector.

- text_size:

  Text size.

- plot_theme:

  Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw the panel border. `NULL` preserves the selected theme
  default.

## Value

A ggplot object.
