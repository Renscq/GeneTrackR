# Create an LDTrack object

`LDTrack()` stores pairwise linkage disequilibrium results, retained
variant metadata, the selected genomic region, and optionally the dosage
matrix and the generated figure. It follows the same lightweight S3
list-class style used by `VariantTrack` and `HapVariant`.

## Usage

``` r
LDTrack(
  data = NULL,
  matrix = NULL,
  variants = NULL,
  region = list(),
  genotype = NULL,
  figure = NULL,
  meta = list(),
  plot = NULL
)
```

## Arguments

- data:

  A data.frame/data.table of pairwise LD statistics.

- matrix:

  A symmetric LD matrix. Diagonal values are usually 1.

- variants:

  Variant metadata used in the calculation.

- region:

  A list with `chrom`, `start`, and `end`.

- genotype:

  Optional numeric genotype dosage matrix.

- figure:

  Optional figure object produced by
  [`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md).

- meta:

  A metadata list.

- plot:

  Deprecated alias for `figure`, kept only for compatibility with older
  development versions.

## Value

An LDTrack object.
