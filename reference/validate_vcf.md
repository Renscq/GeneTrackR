# Validate a VariantTrack object

`validate_vcf()` validates genome-level variant records stored in a
`VariantTrack` object. It checks chromosome and position fields, allele
fields, duplicated variant IDs, and consistency between `pos`, `start`,
and `end` when those columns are available.

## Usage

``` r
validate_vcf(object)
```

## Arguments

- object:

  A VariantTrack object.

## Value

A validation list with `invalid_records`, `invalid_summary`, and
`warnings`.
