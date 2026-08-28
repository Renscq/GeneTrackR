# Write a VariantTrack object

Backward-compatible wrapper around
[`write_vcf()`](https://renscq.github.io/GeneTrackR/reference/write_vcf.md).

## Usage

``` r
write_variant_track(object, file, overwrite = FALSE, include_header = TRUE)
```

## Arguments

- object:

  A `VariantTrack` object created by
  [`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md).

- file:

  Output VCF file path.

- overwrite:

  Whether to overwrite an existing file.

- include_header:

  Whether to include minimal VCF header lines.
