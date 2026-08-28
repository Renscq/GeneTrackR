# Read a GFF3 file as a FeatureTrack

Read a GFF3 file as a FeatureTrack

## Usage

``` r
read_gff(
  file,
  feature_types = NULL,
  verbose = TRUE,
  progress = interactive() && isTRUE(verbose)
)
```

## Arguments

- file:

  GFF3 file path.

- feature_types:

  Optional feature types to keep, such as `gene`, `mRNA`, `exon`, or
  `CDS`.

- verbose:

  Whether to print progress messages.

- progress:

  Whether to show a stage progress bar. The parser reports major stages
  using the same style as
  [`read_genepred()`](https://renscq.github.io/GeneTrackR/reference/read_genepred.md).

## Value

A FeatureTrack object.

## Examples

``` r
gff_file <- system.file("extdata", "gtr_demo.gff3", package = "GeneTrackR")
gff <- read_gff(gff_file, verbose = FALSE, progress = FALSE)
gff
#> <Feature>
#>   records    : 234
#>   genes      : 20
#>   transcripts: 24
#>   exons      : 76
#>   format     : GFF
#>   coordinate : 1-based closed
head(gff$genes)
#>    gene_id  chrom strand gene_start gene_end n_transcripts gene_type
#>     <char> <char> <char>      <int>    <int>         <int>    <char>
#> 1:   GeneD   chr1      +    1000001  1008000             1    coding
#> 2:   GeneE   chr1      +    3000001  3009000             1    coding
#> 3:   GeneF   chr1      -    5000001  5010000             1    coding
#> 4:   GeneG   chr1      +    7000001  7007000             1    coding
#> 5:   GeneH   chr1      -    9000001  9009000             1    coding
#> 6:   GeneA   chr1      +   12340001 12352000             2    coding
```
