# Summarize a Feature-compatible annotation object

`summary_feature()` is the unified summary API for GenePred, GTF, GFF,
and BED-derived annotation objects.

## Usage

``` r
summary_feature(
  object,
  chrom = NULL,
  start = NULL,
  end = NULL,
  level = c("feature", "gene", "transcript", "exon"),
  by = c("chrom", "type")
)
```

## Arguments

- object:

  A Feature-compatible annotation object.

- chrom:

  Optional chromosome filter.

- start:

  Optional region start.

- end:

  Optional region end.

- level:

  Summary level. Use `feature`, `gene`, `transcript`, or `exon`.

- by:

  Grouping columns used when `level = "feature"`.

## Value

A data.table summary.
