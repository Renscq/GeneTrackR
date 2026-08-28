# Plot a FeatureTrack object

Draws BED/GFF/GTF-derived interval features in a genomic region. This
function is designed to be compatible with
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md).

## Usage

``` r
plot_feature_track(
  track,
  chrom,
  start,
  end,
  mode = c("overlap", "within", "trim"),
  color_by = "auto",
  label_by = c("none", "name", "feature_id", "type"),
  feature_palette = "Paired",
  feature_colors = NULL,
  feature_border_color = NA,
  max_legend_levels = 5,
  text_size = 14,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
)
```

## Arguments

- track:

  A FeatureTrack object.

- chrom:

  Chromosome name.

- start:

  Region start.

- end:

  Region end.

- mode:

  `overlap`, `within`, or `trim`.

- color_by:

  Column used for fill colors. Use `auto` to choose a compact,
  informative grouping automatically. Common manual choices are
  `feature_group`, `type`, `source`, `name`, and `strand`.

- label_by:

  Column used for labels. Use `none` to hide labels.

- feature_palette:

  RColorBrewer palette name for feature fills.

- feature_colors:

  Optional named or unnamed feature fill color vector.

- feature_border_color:

  Rectangle border color. Use NA to hide borders.

- max_legend_levels:

  Maximum number of displayed legend groups when the selected color
  attribute contains many categories.

- text_size:

  Text size.

- plot_theme:

  Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw the panel border. `NULL` preserves the selected theme
  default.

## Value

A ggplot object.
