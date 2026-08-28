# Plot phenotype distributions

Plot phenotype distributions

## Usage

``` r
plot_pheno(
  pheno,
  traits = NULL,
  sample_col = "sample_id",
  bins = 30L,
  text_size = 14
)
```

## Arguments

- pheno:

  A phenotype data.frame/data.table.

- traits:

  Optional trait names. If NULL, all traits are plotted.

- sample_col:

  Sample column name.

- bins:

  Histogram bins for numeric traits.

- text_size:

  Text size.

## Value

A ggplot object.
