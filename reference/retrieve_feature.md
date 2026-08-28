# Retrieve annotation features or sub-objects

`retrieve_feature()` extracts a true sub-object from a
Feature/GenePred-compatible annotation object. By default it returns a
Feature/GenePred-compatible object. Use `as = "data.table"` only when a
plain table is required.

## Usage

``` r
retrieve_feature(
  object,
  pattern = NULL,
  level = c("feature", "gene", "transcript", "exon"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  mode = c("overlap", "within", "trim"),
  gene_id = NULL,
  transcript_id = NULL,
  type = NULL,
  ignore_case = TRUE,
  fixed = FALSE,
  as = c("Feature", "data.table")
)
```

## Arguments

- object:

  A Feature-compatible annotation object.

- pattern:

  Optional pattern used for ID/name matching.

- level:

  Output table level when `as = "data.table"`. One of feature, gene,
  transcript, or exon.

- chrom:

  Optional chromosome.

- start:

  Optional start coordinate.

- end:

  Optional end coordinate.

- mode:

  Region matching mode. `overlap` keeps overlapping
  transcripts/features, `within` keeps fully contained records, and
  `trim` clips gene-model records to the query interval.

- gene_id:

  Optional gene ID filter. Exact matching is used.

- transcript_id:

  Optional transcript ID filter. Exact matching is used.

- type:

  Optional feature type filter for feature-table retrieval.

- ignore_case:

  Logical. Whether pattern matching ignores case.

- fixed:

  Logical. Whether pattern is fixed text.

- as:

  Output type. Default `Feature` returns a sub-object. Use `data.table`
  for a plain table.

## Value

A Feature/GenePred-compatible object or a data.table.
