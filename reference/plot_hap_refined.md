# Plot phenotype distributions for refined haplotypes

Deprecated compatibility alias of
[`plot_refined_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_refined_hap_pheno.md).

## Usage

``` r
plot_hap_refined(refined_hap, phenotype, ...)
```

## Arguments

- refined_hap:

  A `HapRefined` object returned by
  [`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md)
  or a refined `HapVariant` object.

- phenotype:

  A phenotype table returned by
  [`read_pheno()`](https://renscq.github.io/GeneTrackR/reference/read_pheno.md)
  or a compatible data.frame.

- ...:

  Additional parameters passed to
  [`plot_refined_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_refined_hap_pheno.md).

## Value

A `GeneTrackRPhenoPlot` object returned by
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md).
