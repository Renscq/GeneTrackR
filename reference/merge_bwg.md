# Merge signal track objects

Merge signal track objects

## Usage

``` r
merge_bwg(
  ...,
  sample_conflict = c("error", "rename", "sum", "mean", "keep_first"),
  require_same_norm = TRUE
)
```

## Arguments

- ...:

  BwgTrack objects.

- sample_conflict:

  Conflict strategy for duplicated sample IDs. `error` stops, `rename`
  makes sample IDs unique, `sum`/`mean` combines records with identical
  genomic intervals, and `keep_first` keeps the first occurrence.

- require_same_norm:

  Require identical normalization labels before merging.

## Value

A merged BwgTrack object.

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
  mode = "memory"
)
riboseq <- read_bwg(
  system.file(
    "extdata",
    c("gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"),
    package = "GeneTrackR"
  ),
  sample_names = c("Ribo_seq_plus", "Ribo_seq_minus"),
  strand = c("+", "-"),
  mode = "memory"
)
signal_all <- merge_bwg(rnaseq, riboseq)
} # }
```
