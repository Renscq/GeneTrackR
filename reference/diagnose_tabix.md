# Diagnose tabix availability for bedGraph tracks

Diagnose tabix availability for bedGraph tracks

## Usage

``` r
diagnose_tabix(object)
```

## Arguments

- object:

  A BwgTrack object created by
  [`read_bwg()`](https://renscq.github.io/GeneTrackR/reference/read_bwg.md).

## Value

A data.table reporting index presence, backend availability, and whether
indexed querying is enabled.

## Examples

``` r
if (FALSE) { # \dontrun{
rnaseq <- read_bwg(
  system.file(
    "extdata",
    c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
    package = "GeneTrackR"
  ),
  sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
  strand = c("+", "-"),
  mode = "lazy"
)
diagnose_tabix(rnaseq)
} # }
```
