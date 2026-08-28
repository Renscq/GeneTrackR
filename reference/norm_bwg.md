# Normalize signal tracks

Normalize signal tracks

## Usage

``` r
norm_bwg(
  object,
  method = c("none", "RPM", "CPM", "scale", "custom"),
  library_size = NULL,
  scale_factor = NULL,
  custom_factor = NULL
)
```

## Arguments

- object:

  A BwgTrack object.

- method:

  Normalization method. Use `RPM`/`CPM` with `library_size`, `scale`
  with `scale_factor`, or `custom` with `custom_factor`.

- library_size:

  Optional named library sizes.

- scale_factor:

  Optional named scale factors.

- custom_factor:

  Optional named custom factors.

## Value

A normalized BwgTrack object.

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
norm_bwg(
  rnaseq,
  method = "scale",
  scale_factor = c(RNA_seq_plus = 0.5, RNA_seq_minus = 0.5)
)
} # }
```
