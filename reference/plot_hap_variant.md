# Plot gene variants and haplotype table

Draws a haplotype-variant figure with a compact gene model track,
variant markers, connector lines, and a haplotype genotype table.
Genomic gaps are compressed by mapping variant positions to evenly
spaced table columns while retaining the relative order of gene model
features.

## Usage

``` r
plot_hap_variant(
  hap,
  annotation = NULL,
  show_gene_model = TRUE,
  min_hap_samples = 5L,
  show_reference_row = TRUE,
  variant_label = c("variant_id", "pos"),
  show_gene_pos_axis = TRUE,
  gene_pos_axis_n = 5L,
  gene_pos_axis_label = NULL,
  gene_pos_x_angle = 0,
  gene_track_legend_position = c("right", "top", "none"),
  direction_mode = c("transcript", "gene", "end", "none"),
  text_size = 14,
  table_x_angle = 90,
  genotype_text_size = 3.2,
  gene_track_height = 1.25,
  connector_height = 0.35,
  table_height = NULL,
  exon_height = 0.22,
  cds_height = 0.44,
  gene_palette = "Paired",
  gene_colors = NULL,
  gene_border_color = NA,
  table_palette = "Paired",
  table_colors = NULL,
  table_alpha = 0.6,
  reference_fill = "white",
  variant_palette = "Paired",
  variant_colors = NULL,
  variant_alpha = 0.6,
  show_variant_marker = TRUE,
  variant_marker_size = 2.8
)
```

## Arguments

- hap:

  A HapVariant object from
  [`hap_variant()`](https://renscq.github.io/GeneTrackR/reference/hap_variant.md).

- annotation:

  Optional Feature/GenePred annotation object. When supplied, a compact
  gene track is drawn above the haplotype table.

- show_gene_model:

  Logical. Whether to draw the gene track when `annotation` is supplied.

- min_hap_samples:

  Minimum sample number required for a haplotype group to be displayed.

- show_reference_row:

  Logical. Whether to add two reference rows showing REF and ALT alleles
  for each variant.

- variant_label:

  Column used for variant labels. One of `variant_id`, `pos`, or an
  existing column in `hap$variants`.

- show_gene_pos_axis:

  Logical. Whether to show genomic coordinate labels above the gene
  track.

- gene_pos_axis_n:

  Approximate number of genomic coordinate ticks above the gene track.

- gene_pos_axis_label:

  Optional x-axis title for genomic coordinate labels.

- gene_pos_x_angle:

  Angle of gene-position x-axis labels. Default is 0.

- gene_track_legend_position:

  Legend position for variant-type markers in the gene track. One of
  `right`, `top`, or `none`.

- direction_mode:

  Gene-strand arrow style for the compact gene track. `transcript` draws
  one arrow per transcript, `gene` draws one arrow per gene, `end` draws
  a short arrow at the directional end of each gene, and `none` hides
  direction arrows.

- text_size:

  Base text size for gene-track labels, axes, legends, and table axes.

- table_x_angle:

  Angle of haplotype table x-axis labels. Default is 90.

- genotype_text_size:

  Genotype cell text size only. Default is 3.2.

- gene_track_height:

  Relative height of the gene track panel.

- connector_height:

  Relative height of the connector panel.

- table_height:

  Relative height of the haplotype table panel.

- exon_height:

  Height of exon boxes in the compact gene track.

- cds_height:

  Height of CDS boxes in the compact gene track.

- gene_palette:

  RColorBrewer palette name used for gene model feature fills.

- gene_colors:

  Optional custom fill colors for gene model features.

- gene_border_color:

  Border color for gene model feature boxes.

- table_palette:

  RColorBrewer palette name used for haplotype table fills. For
  allele-string haplotypes, the first five colors are always assigned in
  the fixed order A, T, C, G, and indel.

- table_colors:

  Optional custom fill colors for haplotype table genotypes. For
  allele-string haplotypes, use a vector named `A`, `T`, `C`, `G`, and
  `indel`, or supply five unnamed colors in that order.

- table_alpha:

  Alpha value for haplotype table fill colors.

- reference_fill:

  Background fill color for the REF and ALT reference rows.

- variant_palette:

  RColorBrewer palette name used for solid variant-type triangle marker
  colors.

- variant_colors:

  Optional custom colors for solid variant-type triangle markers.

- variant_alpha:

  Alpha value for solid variant-type triangle marker colors.

- show_variant_marker:

  Logical. Whether to draw natural-variant triangle markers.

- variant_marker_size:

  Size of natural-variant triangle markers. Set to 0 to hide markers.

## Value

A patchwork object with attributes `plot_data`, `variant_data`, and
`gene_data`.

## Examples

``` r
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
anno_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string", min_variant_number = 1)
#> [GeneTrackR] Retrieved variants: 11.
plot_hap_variant(hap, annotation = anno, min_hap_samples = 1)

plot_hap_variant(
  hap,
  annotation = anno,
  min_hap_samples = 1,
  table_x_angle = 90,
  table_palette = "Paired"
)
```
