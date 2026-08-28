# Merge annotation feature objects

`merge_feature()` combines annotations read from GenePred, GTF, GFF, or
BED files into one unified `Feature` object for downstream retrieval and
plotting. Inputs may represent different genes, genomic regions,
annotation formats, or independently retrieved subsets.

Duplicate identifiers are detected across input objects. A warning
reports duplicated gene, transcript, and feature IDs before the selected
conflict strategy is applied.

With `conflict = "deduplicate"`, input order defines precedence. The
first complete copy of a duplicated transcript is retained, duplicated
feature records are removed, and gene ranges are recalculated from all
retained non-duplicated transcripts. With `conflict = "rename"`, only
conflicting identifiers in later inputs are renamed, and all
hierarchical references are updated together.

## Usage

``` r
merge_feature(
  ...,
  source_names = NULL,
  conflict = c("deduplicate", "rename", "keep_all", "error", "keep_first"),
  rename_prefix = NULL,
  sort = TRUE
)
```

## Arguments

- ...:

  One or more `Feature`, `FeatureTrack`, or `GenePred` objects. A single
  list containing such objects is also accepted.

- source_names:

  Optional unique source labels. If omitted, names from a supplied list
  or named arguments are used where available; remaining labels are
  generated as `track1`, `track2`, and so on. The labels are stored in
  the `track_source` column.

- conflict:

  Duplicate-ID strategy. `"deduplicate"` keeps the first compatible
  annotation model according to input order. `"rename"` renames
  conflicting IDs in later inputs. `"keep_all"` retains conflicts
  unchanged, and `"error"` stops after duplicate detection.
  `"keep_first"` is accepted as a backward-compatible alias of
  `"deduplicate"`.

- rename_prefix:

  Optional character vector with one prefix per input. Used only when
  `conflict = "rename"`. By default, `source_names` followed by an
  underscore are used.

- sort:

  Logical. Whether to sort the merged feature and hierarchy tables.

## Value

A merged `Feature` object. If transcript and exon tables are present,
the result also inherits from `GenePred`.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_a <- retrieve_feature(annotation, gene_id = "GeneA")
gene_b <- retrieve_feature(annotation, gene_id = "GeneB")
merged_genes <- merge_feature(gene_a, gene_b)

merged_formats <- merge_feature(
  list(genepred_annotation, gtf_annotation),
  source_names = c("genePred", "GTF"),
  conflict = "rename"
)
} # }
```
