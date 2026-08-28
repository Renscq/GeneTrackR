# Write Feature annotation data

Writes a unified `Feature` annotation object to GenePred, GenePredExt,
GFF3, GTF, or BED format. This is the main output function for
annotation tracks produced by
[`read_genepred()`](https://renscq.github.io/GeneTrackR/reference/read_genepred.md),
[`read_gff()`](https://renscq.github.io/GeneTrackR/reference/read_gff.md),
[`read_gtf()`](https://renscq.github.io/GeneTrackR/reference/read_gtf.md),
or
[`read_bed()`](https://renscq.github.io/GeneTrackR/reference/read_bed.md).

## Usage

``` r
write_feature(
  object,
  file,
  format = c("auto", "genepred", "genepredext", "gff", "gtf", "bed6", "bed12"),
  coordinate = c("ucsc", "granges"),
  overwrite = FALSE,
  keep_attributes = TRUE,
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

- coordinate:

  Output coordinate system for GenePred/BED-like formats. `ucsc` means
  0-based half-open where applicable; `granges` means 1-based closed.

- overwrite:

  Whether to overwrite an existing file.

- keep_attributes:

  Whether to reuse existing `attribute` strings for GFF/GTF output when
  available.

- sort_output:

  Whether to sort records before writing. Default TRUE. GenePred output
  is sorted by chromosome and transcript start. GFF/GTF output is sorted
  by chromosome, start, and feature hierarchy.

- chrom_order:

  Optional chromosome order used for sorting. It can be NULL, a
  character vector of chromosome names, a data frame whose first column
  contains chromosome names, or a FASTA index `.fai` file path. If NULL,
  a natural chromosome order is used.

## Value

Invisibly returns the output file path.

## Examples

``` r
if (FALSE) { # \dontrun{
write_feature(gp, "annotation.genePredExt", format = "genepredext")
write_feature(gtf, "annotation.gtf", format = "gtf")
write_feature(gff, "annotation.gff3", format = "gff")
write_feature(gp, "annotation.bed6", format = "bed6")
write_feature(gp, "annotation.bed12", format = "bed12")
} # }
```
