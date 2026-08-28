# Summarize phenotype traits

Summarize phenotype traits

## Usage

``` r
summary_pheno(pheno, sample_col = "sample_id")
```

## Arguments

- pheno:

  A phenotype data.frame/data.table returned by
  [`read_pheno()`](https://renscq.github.io/GeneTrackR/reference/read_pheno.md).

- sample_col:

  Sample column name.

## Value

A data.table summarizing trait type, missing values, and basic
statistics.
