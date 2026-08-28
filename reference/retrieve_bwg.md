# Retrieve signal records from a BwgTrack object

`retrieve_bwg()` is the unified signal retrieval API in GeneTrackR. It
retrieves signal intervals from bedGraph, wig, or bigWig tracks by
genomic region.

## Usage

``` r
retrieve_bwg(
  object,
  chrom = NULL,
  start = NULL,
  end = NULL,
  annotation = NULL,
  gene_id = NULL,
  transcript_id = NULL,
  upstream = 0L,
  downstream = 0L,
  strand_aware = TRUE,
  samples = NULL,
  strand = c("ignore", "+", "-", "both", "auto"),
  strand_policy = c("ignore_unstranded", "strict"),
  as = c("data.table", "BwgTrack", "GRanges"),
  verbose = FALSE,
  progress = interactive() && isTRUE(verbose),
  keep_empty_samples = FALSE,
  tabix_empty_fallback = NULL
)
```

## Arguments

- object:

  A BwgTrack object returned by
  [`read_bwg()`](https://renscq.github.io/GeneTrackR/reference/read_bwg.md).

- chrom:

  Chromosome name. Required for direct region queries.

- start:

  Region start in 1-based closed coordinates. Required for direct region
  queries.

- end:

  Region end in 1-based closed coordinates. Required for direct region
  queries.

- annotation:

  Optional GenePred/Feature annotation object used for
  gene/transcript-aware retrieval.

- gene_id:

  Optional gene ID. When supplied, `annotation` is used to resolve the
  gene range.

- transcript_id:

  Optional transcript ID. When supplied, `annotation` is used to resolve
  the transcript range.

- upstream:

  Upstream flanking length in bp for gene/transcript queries.

- downstream:

  Downstream flanking length in bp for gene/transcript queries.

- strand_aware:

  Logical. Whether upstream/downstream should follow gene/transcript
  strand direction.

- samples:

  Optional character vector of sample IDs to retrieve. If NULL, all
  samples are queried.

- strand:

  Strand selector. Use `"ignore"` for unstranded retrieval, `"+"` or
  `"-"` for strand-specific tracks, `"both"` for both strands, or
  `"auto"` when called by higher-level gene-aware plotting functions.

- strand_policy:

  How to handle strand filtering for unstranded signal files.
  `"ignore_unstranded"` returns unstranded tracks for any strand
  request; `"strict"` only returns explicitly matching strand records.

- as:

  Output type. `"data.table"` returns a signal table, `"BwgTrack"`
  returns an in-memory BwgTrack sub-object, and `"GRanges"` returns a
  GRanges object.

- verbose:

  Logical. Whether to print progress messages.

- progress:

  Logical. Whether to show a text progress bar for multi-sample queries.

- keep_empty_samples:

  Logical. If TRUE, samples without signal in the requested region are
  returned as zero-valued placeholder intervals. This is useful for
  preserving empty facets in plots.

- tabix_empty_fallback:

  Optional logical. If TRUE, empty tabix results are verified by
  full-file fread. Keep FALSE for speed unless debugging tabix
  coordinate/index issues.

## Value

A data.table, BwgTrack object, or GRanges object.

## Examples

``` r
if (FALSE) { # \dontrun{
rnaseq_files <- system.file(
  "extdata",
  c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
  package = "GeneTrackR"
)
rnaseq <- read_bwg(
  rnaseq_files,
  format = "bedgraph",
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

anno <- read_genepred(
  system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
  format = "genePredExt",
  verbose = FALSE
)
gene_signal <- retrieve_bwg(
  rnaseq,
  annotation = anno,
  gene_id = "GeneA",
  upstream = 500,
  downstream = 500,
  strand = "auto"
)
} # }
```
