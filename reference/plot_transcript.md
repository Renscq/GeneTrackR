# Plot a transcript structure

Plot a transcript structure

## Usage

``` r
plot_transcript(
  object,
  transcript_id,
  coordinate = c("transcript", "genomic"),
  show_cds = TRUE,
  cds_height = 0.5,
  utr_height = 0.25,
  direction_mode = c("transcript", "gene", "end", "none"),
  highlight = NULL,
  gene_palette = "Paired",
  gene_colors = NULL,
  gene_border_color = NA,
  label_position = c("axis", "feature"),
  label_by = c("gene", "transcript"),
  text_size = 14,
  plot_theme = c("bw", "classic", "light", "minimal"),
  show_panel_border = NULL
)
```

## Arguments

- object:

  A GenePred object.

- transcript_id:

  Transcript ID.

- coordinate:

  transcript removes introns, genomic keeps genomic scale.

- show_cds:

  Whether to distinguish CDS and UTR segments.

- cds_height:

  Vertical thickness of CDS rectangles.

- utr_height:

  Vertical thickness of UTR/non-coding exon rectangles.

- direction_mode:

  Direction-arrow style. `transcript` draws one arrow per transcript,
  `gene` draws one arrow per gene, `end` draws one short arrow at the
  directional end of each gene, and `none` hides direction arrows.

- highlight:

  Optional data frame used to shade genomic or transcript intervals. It
  must contain `start` and `end` columns. Optional columns such as
  `label` or `group` are allowed but are not required. In
  `coordinate = "genomic"`, `start` and `end` are genomic positions; in
  `coordinate = "transcript"`, they are spliced transcript coordinates.

- gene_palette:

  RColorBrewer palette name used for discrete fills. If the number of
  discrete groups exceeds the palette maximum, colors are automatically
  interpolated.

- gene_colors:

  Optional custom fill colors for gene model features. Use a named
  vector such as `c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99")`
  to map colors explicitly. Unnamed colors are matched to the fixed
  gene-model levels `UTR`, `CDS`, and `exon`, so colors remain stable
  even when only one feature type is present.

- gene_border_color:

  Optional rectangle border color. Use `NA` to hide borders, or a color
  such as `"black"` or `"grey30"` to draw feature outlines.

- label_position:

  Where to draw feature labels. `axis` draws labels on the y axis and
  `feature` draws labels at the center of each gene/transcript
  structure.

- label_by:

  Which identifier to use for feature labels. Use `gene` for gene IDs or
  `transcript` for transcript IDs.

- text_size:

  Text size in points for axis text, axis titles, and legends.

- plot_theme:

  Base ggplot2 theme. Use `bw`, `classic`, `light`, or `minimal`.

- show_panel_border:

  Whether to draw the panel border. `NULL` preserves the selected theme
  default.

## Value

A ggplot object.

## Details

Common feature names used by `gene_colors` are `UTR`, `CDS`, and `exon`.
`highlight` must contain `start` and `end` columns. In genomic
coordinate plots these are genomic positions; in transcript coordinate
plots they are spliced transcript positions. See also
[GeneTrackR-advanced-parameters](https://renscq.github.io/GeneTrackR/reference/GeneTrackR-advanced-parameters.md).

## Examples

``` r
if (FALSE) { # \dontrun{
gp <- read_genepred("annotation.genePredExt", format = "genePredExt")
plot_transcript(gp, transcript_id = "EAA28783")
plot_transcript(
  gp,
  transcript_id = "EAA28783",
  coordinate = "transcript",
  gene_colors = c(UTR = "#b2df8a", CDS = "#33a02c", exon = "#fb9a99"),
  gene_border_color = "black",
  highlight = data.frame(start = 100, end = 250),
  label_position = "feature",
)
} # }
```
