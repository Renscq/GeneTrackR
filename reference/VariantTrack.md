# Create a VariantTrack object

`VariantTrack()` stores VCF-like variant records as genomic point
features. Coordinates are stored internally as 1-based positions with
`start == end`.

## Usage

``` r
VariantTrack(data = NULL, meta = list())
```

## Arguments

- data:

  A data.frame or data.table containing at least `chrom` and `pos`.

- meta:

  A list of metadata.

## Value

A VariantTrack object.
