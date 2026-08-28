# Test whether a bedGraph tabix query returns records without full-file fallback

Test whether a bedGraph tabix query returns records without full-file
fallback

## Usage

``` r
benchmark_tabix_query(object, chrom, start, end, samples = NULL)
```

## Arguments

- object:

  A BwgTrack object.

- chrom:

  Chromosome name.

- start:

  Region start.

- end:

  Region end.

- samples:

  Optional sample IDs.

## Value

A data.table with per-sample query timing and number of records.
