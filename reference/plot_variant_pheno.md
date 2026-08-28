# Plot phenotype values grouped by a single variant genotype

Draws phenotype distributions for genotype groups at one variant site.
The function is designed as the single-variant counterpart of
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md).
Groups are ordered by sample number from left to right, and fill colors
are mapped to the median phenotype value of each genotype group.

## Usage

``` r
plot_variant_pheno(
  variant,
  phenotype,
  traits = NULL,
  sample_col = "sample_id",
  variant_id = NULL,
  chrom = NULL,
  pos = NULL,
  start = NULL,
  end = NULL,
  samples = NULL,
  genotype_mode = c("code", "string"),
  missing_genotype = NA_character_,
  min_group_samples = 2L,
  plot_type = c("violin_boxplot", "violin", "boxplot"),
  test_method = c("t.test", "wilcox.test", "ks.test"),
  p_adjust = "BH",
  p_label = c("stars", "number", "both"),
  p_cutoff = 0.05,
  p_value_type = c("raw", "adjusted"),
  show_signif_only = TRUE,
  show_points = FALSE,
  show_outliers = FALSE,
  fill_palette = "Paired",
  fill_colors = NULL,
  fill_alpha = 0.75,
  violin_width = 0.9,
  box_width = 0.18,
  bracket_step = 0.08,
  bracket_tip_fraction = 0.12,
  x_text_angle = 90,
  strip_label_width = 24,
  strip_fill = "white",
  strip_border_color = NULL,
  strip_text_lineheight = 0.9,
  text_size = 14
)
```

## Arguments

- variant:

  A VariantTrack object, a VCF file path, a HapVariant object, or a
  VCF-like data.frame/data.table containing one or more variants with
  genotype sample columns.

- phenotype:

  A phenotype table returned by
  [`read_pheno()`](https://renscq.github.io/GeneTrackR/reference/read_pheno.md)
  or a compatible data.frame.

- traits:

  Phenotype trait names. If NULL, all numeric traits are used.

- sample_col:

  Sample column name in phenotype table.

- variant_id:

  Optional variant ID to select.

- chrom:

  Optional chromosome name used to select a variant.

- pos:

  Optional 1-based variant position. This is equivalent to setting both
  `start` and `end` to the same value.

- start:

  Optional 1-based start position for region-based selection.

- end:

  Optional 1-based end position for region-based selection.

- samples:

  Optional sample names to keep.

- genotype_mode:

  Genotype representation. `code` converts genotypes to 0/1 states,
  where 0 means reference genotype and 1 means any alternate allele is
  present. `string` converts genotypes to compact allele labels; long
  InDel alleles are compressed as `iN`, where `N` is allele length.

- missing_genotype:

  Missing genotype label. Default is `NA_character_`.

- min_group_samples:

  Minimum sample number required for a genotype group.

- plot_type:

  Plot type. One of `violin`, `boxplot`, or `violin_boxplot`.

- test_method:

  Pairwise test method. One of `t.test`, `wilcox.test`, or `ks.test`.

- p_adjust:

  P-value adjustment method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- p_label:

  P-value label style. One of `stars`, `number`, or `both`.

- p_cutoff:

  Significance cutoff used for displaying pairwise comparisons.

- p_value_type:

  Which p-value is used for filtering and labeling. One of `raw` or
  `adjusted`.

- show_signif_only:

  Logical. Whether to only display significant comparisons.

- show_points:

  Logical. Whether to show sample points.

- show_outliers:

  Logical. Whether to show boxplot outliers.

- fill_palette:

  RColorBrewer palette name used for median-based genotype fills.

- fill_colors:

  Optional custom fill colors for genotype groups.

- fill_alpha:

  Alpha value for violin/boxplot fill colors.

- violin_width:

  Violin plot width.

- box_width:

  Boxplot width.

- bracket_step:

  Fraction of y-range used to separate significance brackets.

- bracket_tip_fraction:

  Fraction of bracket vertical spacing used for the short downward
  bracket tips.

- x_text_angle:

  Rotation angle for genotype labels on the x-axis.

- strip_label_width:

  Maximum character width for wrapping long facet strip labels.

- strip_fill:

  Strip background fill color. Default is white.

- strip_border_color:

  Strip border color. Default NULL removes the strip border.

- text_size:

  Text size.

## Value

A list with `figure` and `pvalue` elements. Additional elements include
`summary`, `bracket`, `plot_data`, and `variant_data`.

## Examples

``` r
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
pheno_file <- system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR")
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
pheno <- read_pheno(pheno_file, verbose = FALSE)
plot_variant_pheno(vcf, phenotype = pheno, variant_id = "varA03",
                   traits = "protein_content", min_group_samples = 3)

plot_variant_pheno(vcf, phenotype = pheno, chrom = "chr1", pos = 12342550,
                   traits = "protein_content", genotype_mode = "string",
                   min_group_samples = 3, test_method = "wilcox.test")
```
