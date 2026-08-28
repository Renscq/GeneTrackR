# Write a FeatureTrack object

Backward-compatible wrapper around
[`write_feature()`](https://renscq.github.io/GeneTrackR/reference/write_feature.md).

## Usage

``` r
write_feature_track(
  object,
  file,
  format = c("bed12", "bed6", "gff", "gtf", "genepred", "genepredext"),
  sort_output = TRUE,
  chrom_order = NULL
)
```

## Arguments

- object:

  A `Feature`, `FeatureTrack`, or GenePred-compatible object.

- file:

  Output file path.

- format:

  Output format. Use `auto`, `genepred`, `genepredext`, `gff`, `gtf`,
  `bed6`, or `bed12`. `bed6` writes six-column BED intervals. For
  gene-model objects, BED6 is transcript-level. For generic interval
  objects, BED6 is feature-level. `bed12` writes transcript-level BED12
  gene models with exon blocks.

- sort_output:

  Whether to sort records before writing. Default TRUE. GenePred output
  is sorted by chromosome and transcript start. GFF/GTF output is sorted
  by chromosome, start, and feature hierarchy.

- chrom_order:

  Optional chromosome order used for sorting. It can be NULL, a
  character vector of chromosome names, a data frame whose first column
  contains chromosome names, or a FASTA index `.fai` file path. If NULL,
  a natural chromosome order is used.
