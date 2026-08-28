# Write BwgTrack signal data

Writes a `BwgTrack` object to bedGraph, wig, or bigWig files. bedGraph
and wig output are implemented directly in R. bigWig output is written
directly with the bundled third-party libBigWig library and requires
chromosome sizes.

For in-memory `BwgTrack` objects, signal records are written from
`object$data`. For lazy objects, direct conversion is not possible
because signal records are not loaded; if the requested format matches
the original file format, the original files are copied to the output
directory.

## Usage

``` r
write_bwg(
  object,
  outdir,
  format = c("bedgraph", "wig", "bigwig"),
  samples = NULL,
  chrom_sizes = NULL,
  overwrite = FALSE,
  compress = FALSE
)
```

## Arguments

- object:

  A `BwgTrack` object.

- outdir:

  Output directory.

- format:

  Output format. One of `bedgraph`, `wig`, or `bigwig`.

- samples:

  Optional sample IDs to write. Default writes all samples.

- chrom_sizes:

  Chromosome sizes for bigWig output. Can be a file path or a
  data.frame/data.table with two columns: chromosome and size.

- overwrite:

  Whether to overwrite existing files.

- compress:

  Whether to gzip-compress bedGraph or wig output.

## Value

Invisibly returns a data.table containing sample IDs and output files.

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
write_bwg(rnaseq, outdir = tempdir(), format = "bedgraph", overwrite = TRUE)
} # }
```
