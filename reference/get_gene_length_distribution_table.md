# Extract length records from a gene annotation object

Build a tidy feature-length table from a Feature/GenePred-compatible
annotation object. The function supports gene spans, transcript lengths,
exon lengths, CDS lengths, total UTR lengths, 5' UTR lengths, and 3' UTR
lengths. CDS and UTR lengths can be summarized per transcript or
returned as individual genomic segments.

## Usage

``` r
get_gene_length_distribution_table(
  object,
  feature = c("all", "gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr"),
  unit = c("auto", "transcript", "segment"),
  transcript_length = c("spliced", "genomic"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  mode = c("overlap", "within", "trim"),
  keep_zero = FALSE
)

get_genepred_length_table(...)
```

## Arguments

- object:

  A Feature or GenePred-compatible annotation object.

- feature:

  Feature type to extract. Use `all` to return gene, transcript, exon,
  CDS, UTR, 5' UTR, and 3' UTR records.

- unit:

  Output unit. `auto` uses gene-level records for genes,
  transcript-level records for transcripts, transcript-level totals for
  CDS/UTR, and segment-level records for exons. `transcript` summarizes
  exonic features per transcript where applicable. `segment` returns
  individual exon/CDS/UTR genomic segments where applicable.

- transcript_length:

  Transcript length definition. `spliced` uses the sum of exon lengths.
  `genomic` uses transcript span length.

- chrom:

  Optional chromosome filter.

- start:

  Optional region start in 1-based closed coordinates.

- end:

  Optional region end in 1-based closed coordinates.

- mode:

  Region selection mode passed to
  [`retrieve_feature()`](https://renscq.github.io/GeneTrackR/reference/retrieve_feature.md).

- keep_zero:

  Logical. Whether to keep zero-length transcript-level CDS/UTR records.
  Default is `FALSE`.

## Value

A data.table with feature-length records.

## Details

Use `unit = "segment"` to inspect each exon/CDS/UTR segment separately,
and `unit = "transcript"` to summarize feature lengths per transcript.
`transcript_length = "spliced"` uses exon lengths, whereas `"genomic"`
uses transcript span length including introns.

## Examples

``` r
if (FALSE) { # \dontrun{
gp <- read_genepred("annotation.genePredExt", format = "genePredExt")
get_gene_length_distribution_table(gp, feature = "cds", unit = "transcript")
get_gene_length_distribution_table(gp, feature = "five_utr", unit = "segment")
get_gene_length_distribution_table(
  gp,
  feature = "transcript",
  transcript_length = "spliced",
  chrom = "chr1",
  start = 1,
  end = 1000
)
} # }
```
