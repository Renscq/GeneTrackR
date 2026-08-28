# Summarize a VariantTrack object

`summary_vcf()` is the unified summary API for VCF-derived VariantTrack
objects. When no genomic range is supplied, the full in-memory object is
summarized. For lazy tracks, an omitted range causes the source VCF to
be read before summarization.

## Usage

``` r
summary_vcf(
  object,
  chrom = NULL,
  start = NULL,
  end = NULL,
  by = c("chrom", "variant_type")
)
```

## Arguments

- object:

  A VariantTrack object.

- chrom:

  Optional chromosome filter. May be used without `start`/`end`.

- start:

  Optional 1-based start coordinate. Must be paired with `end`.

- end:

  Optional 1-based end coordinate. Must be paired with `start`.

- by:

  Grouping columns. Default `c("chrom", "variant_type")`.

## Value

A data.table summary.

## Examples

``` r
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
summary_vcf(vcf)
#>     chrom variant_type n_variants
#>    <char>       <char>      <int>
#> 1:   chr1          DEL          1
#> 2:   chr1          INS          2
#> 3:   chr1          SNP         33
#> 4:   chr2          SNP         20
summary_vcf(vcf, chrom = "chr1")
#>     chrom variant_type n_variants
#>    <char>       <char>      <int>
#> 1:   chr1          DEL          1
#> 2:   chr1          INS          2
#> 3:   chr1          SNP         33
```
