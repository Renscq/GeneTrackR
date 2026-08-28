# Create a standardized Feature object

`Feature()` is the unified annotation class used by GeneTrackR. BED,
GFF, GTF, and GenePred annotations are all normalized to a
`Feature`-compatible object with a flat interval table in `object$data`.
GenePred-compatible inputs additionally store derived `genes`,
`transcripts`, and `exons` tables, allowing the same object to be used
by gene model plotting functions.

Required coordinate convention is 1-based closed intervals.

## Usage

``` r
Feature(
  data,
  meta = list(),
  genes = NULL,
  transcripts = NULL,
  exons = NULL,
  validation = make_empty_validation()
)
```

## Arguments

- data:

  A data.frame or data.table containing standardized interval columns.
  Required columns are `chrom`, `start`, and `end`. Recommended columns
  include `feature_id`, `name`, `type`, `score`, `strand`, `source`,
  `gene_id`, `transcript_id`, `parent_id`, `level`, and `gene_type`.

- meta:

  Metadata list such as source file, original format, and coordinate
  system.

- genes:

  Optional gene-level table.

- transcripts:

  Optional transcript-level table.

- exons:

  Optional exon-level table.

- validation:

  Optional validation list.

## Value

A Feature object, also inheriting from FeatureTrack.
