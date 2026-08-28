# Bin signal tracks into fixed-width windows

Bin signal tracks into fixed-width windows

## Usage

``` r
bin_bwg(data, bin_size = 50L)
```

## Arguments

- data:

  Signal table returned by retrieve_bwg. The table must contain
  `sample_id`, `chrom`, `start`, `end`, and `value` columns.

- bin_size:

  Bin size in bases.

## Value

A binned signal table.

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
dt <- retrieve_bwg(
  rnaseq,
  chrom = "chr1",
  start = 12339001,
  end = 12352000,
  strand = "+"
)
bin_bwg(dt, bin_size = 50)
} # }
```
