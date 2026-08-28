# Write a VariantTrack object to VCF

Write a VariantTrack object to VCF

## Usage

``` r
write_vcf(object, file, overwrite = FALSE, include_header = TRUE)
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

## Value

Invisibly returns the output file path.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- read_vcf("variants.vcf")
write_vcf(vcf, "variants.out.vcf")
} # }
```
