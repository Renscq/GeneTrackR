# Create a GenePred object

Create a GenePred object

## Usage

``` r
GenePred(
  transcripts,
  exons,
  genes = NULL,
  meta = list(),
  validation = make_empty_validation()
)
```

## Arguments

- transcripts:

  Transcript-level annotation table.

- exons:

  Exon-level annotation table.

- genes:

  Gene-level annotation table.

- meta:

  Metadata list.

- validation:

  Validation result list.

## Value

A GenePred object.
