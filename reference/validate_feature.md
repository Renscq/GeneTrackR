# Validate a Feature annotation object

`validate_feature()` is the unified validation entry for annotation
data. It validates BED/GFF/GTF/GenePred-derived `Feature` objects. If
gene-model tables are present, transcript/exon consistency is checked in
addition to the flat standardized feature table.

## Usage

``` r
validate_feature(object, check_gene_model = TRUE)
```

## Arguments

- object:

  A Feature-compatible annotation object.

- check_gene_model:

  Logical. Whether to check derived gene-model tables when transcripts
  and exons are available.

## Value

A validation list with `invalid_records`, `invalid_summary`, and
`warnings`.
