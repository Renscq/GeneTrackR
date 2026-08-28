# Read a VCF file as a VariantTrack

Backward-compatible alias of
[`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md).

## Usage

``` r
read_vcf_track(
  file,
  keep_genotype = TRUE,
  mode = c("auto", "memory", "lazy"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  verbose = TRUE,
  progress = interactive() && isTRUE(verbose)
)
```

## Arguments

- file:

  VCF file path. Plain VCF, gzip-compressed VCF, and bgzip VCF are
  supported.

- keep_genotype:

  Logical. Whether to keep FORMAT and sample genotype columns.

- mode:

  Reading mode. `memory` parses records immediately, `lazy` stores
  indexed VCF metadata only, and `auto` uses lazy mode for indexed VCF
  files when no region is supplied.

- chrom:

  Optional chromosome name for indexed regional reading.

- start:

  Optional 1-based region start for indexed regional reading.

- end:

  Optional 1-based region end for indexed regional reading.

- verbose:

  Logical. Whether to print progress messages.

- progress:

  Logical. Whether to print a compact stage-level progress indicator.

## Value

A VariantTrack object.
