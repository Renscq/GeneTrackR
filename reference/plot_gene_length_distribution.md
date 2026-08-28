# Plot gene annotation feature length distributions

Plot the length distribution of genes, transcripts, exons, CDS, UTRs, 5'
UTRs, and 3' UTRs from any Feature/GenePred-compatible annotation
object. The function is designed for annotation quality control and
comparison of coding versus non-coding feature lengths.

## Usage

``` r
plot_gene_length_distribution(
  object,
  feature = c("all", "gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr"),
  unit = c("auto", "transcript", "segment"),
  transcript_length = c("spliced", "genomic"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  mode = c("overlap", "within", "trim"),
  group_by = c("gene_type", "feature", "strand", "chrom", "none"),
  plot_type = c("density", "histogram", "boxplot", "violin"),
  scale = c("log10", "linear"),
  bins = 60L,
  facet = TRUE,
  keep_zero = FALSE,
  fill_palette = "Paired",
  fill_colors = NULL,
  border_color = NA,
  return_data = FALSE
)

plot_genepred_length_distribution(...)
```

## Arguments

- object:

  A Feature or GenePred-compatible annotation object.

- feature:

  Feature type to plot. Use `all` to facet multiple feature types.

- unit:

  Output unit used for length extraction. See
  [`get_gene_length_distribution_table()`](https://renscq.github.io/GeneTrackR/reference/get_gene_length_distribution_table.md).

- transcript_length:

  Transcript length definition for transcript records.

- chrom:

  Optional chromosome filter.

- start:

  Optional region start in 1-based closed coordinates.

- end:

  Optional region end in 1-based closed coordinates.

- mode:

  Region selection mode passed to
  [`retrieve_feature()`](https://renscq.github.io/GeneTrackR/reference/retrieve_feature.md).

- group_by:

  Grouping variable. Common choices are `gene_type`, `feature`,
  `strand`, and `chrom`.

- plot_type:

  Plot type. One of `density`, `histogram`, `boxplot`, or `violin`.

- scale:

  Length scale. `log10` is recommended for genomic features because
  feature lengths are usually right-skewed.

- bins:

  Number of bins for histogram.

- facet:

  Logical. Whether to facet by feature when multiple features are
  requested.

- keep_zero:

  Logical. Whether to keep zero-length CDS/UTR records.

- fill_palette:

  RColorBrewer palette name used for grouped fills. Default is `Paired`.

- fill_colors:

  Optional custom fill colors. Named vectors are matched to the values
  of `group_by`, for example
  `c(coding = "#1f78b4", `non-coding` = "#a6cee3")`. Unnamed colors are
  matched in plotting order and are automatically extended when needed.

- border_color:

  Optional border color for histograms, boxplots, and violin plots. Use
  `NA` to hide borders.

- return_data:

  Logical. If `TRUE`, return a list containing the plot and the
  underlying length table.

## Value

A ggplot object, or a list with `plot` and `data` when
`return_data = TRUE`.

## Details

`fill_colors` follows the values of `group_by`. For example, if
`group_by = "gene_type"`, names should match `coding` and `non-coding`;
if `group_by = "feature"`, names can include `gene`, `transcript`,
`exon`, `cds`, `utr`, `five_utr`, and `three_utr`. Use
`return_data = TRUE` to inspect the exact groups before assigning named
colors. See also
[GeneTrackR-advanced-parameters](https://renscq.github.io/GeneTrackR/reference/GeneTrackR-advanced-parameters.md).

## Examples

``` r
if (FALSE) { # \dontrun{
gp <- read_genepred("annotation.genePredExt", format = "genePredExt")
plot_gene_length_distribution(gp, feature = "gene", group_by = "gene_type")
plot_gene_length_distribution(
  gp,
  feature = "gene",
  group_by = "gene_type",
  fill_colors = c(coding = "#1f78b4", `non-coding` = "#a6cee3"),
  border_color = NA
)
plot_gene_length_distribution(
  gp,
  feature = "all",
  group_by = "feature",
  plot_type = "density",
  facet = TRUE
)
res <- plot_gene_length_distribution(
  gp,
  feature = "all",
  group_by = "feature",
  return_data = TRUE
)
unique(res$data$feature)
length_table <- get_gene_length_distribution_table(gp, feature = "cds")
head(length_table)
} # }
```
