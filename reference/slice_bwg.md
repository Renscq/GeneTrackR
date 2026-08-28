# Slice signal tracks by genomic region

Slice signal tracks by genomic region

## Usage

``` r
slice_bwg(
  object,
  chrom,
  start,
  end,
  samples = NULL,
  strand = "ignore",
  as = c("BwgTrack", "data.frame", "GRanges")
)
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

- strand:

  Strand selector. For unstranded bigWig/wig tracks, '+' and '-' do not
  filter records.

- as:

  Return type. Use `BwgTrack` to keep the object structure, `data.frame`
  for a plain table, or `GRanges` for genomic interval operations.

## Value

A BwgTrack, data.table, or GRanges object.

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
slice_bwg(
  rnaseq,
  chrom = "chr1",
  start = 12339001,
  end = 12352000,
  strand = "+",
  as = "BwgTrack"
)
} # }
```
