# Create a FeatureTrack object

`FeatureTrack()` is kept as a user-facing alias of
[`Feature()`](https://renscq.github.io/GeneTrackR/reference/Feature.md)
for interval tracks. It returns an object inheriting from both `Feature`
and `FeatureTrack`.

## Usage

``` r
FeatureTrack(
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

A Feature object inheriting from FeatureTrack.
