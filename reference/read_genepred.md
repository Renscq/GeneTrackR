# Read a GenePred annotation file

Read a GenePred annotation file

## Usage

``` r
read_genepred(
  file,
  format = c("auto", "genePred", "genePredExt"),
  coordinate = c("ucsc", "granges"),
  remove_invalid = TRUE,
  report_invalid = TRUE,
  gene_col = c("name2", "name"),
  transcript_col = "name",
  verbose = TRUE,
  progress = interactive() && isTRUE(verbose)
)
```

## Arguments

- file:

  GenePred or GenePredExt file path.

- format:

  Input format. Use auto, genePred, or genePredExt.

- coordinate:

  Input coordinate system. ucsc means 0-based half-open.

- remove_invalid:

  Remove invalid records after validation.

- report_invalid:

  Store invalid records in the object validation slot.

- gene_col:

  Gene identifier source. name2 is recommended for GenePredExt.

- transcript_col:

  Transcript identifier source.

- verbose:

  Logical. Whether to print step-level progress messages. Default TRUE.

- progress:

  Logical. Whether to show a stage-level progress bar. By default, a
  progress bar is shown only in interactive sessions when
  `verbose = TRUE`.

## Value

A GenePred object.

## Details

Standard GenePred uses `name` as transcript ID. GenePredExt commonly
uses `name` as transcript ID and `name2` as gene ID, so the default
`gene_col` tries `name2` first and falls back to `name` when `name2` is
unavailable. `coordinate = "ucsc"` converts 0-based half-open GenePred
coordinates to the package's internal 1-based closed coordinate system.

For large files, `verbose = TRUE` prints major processing stages, while
`progress = TRUE` additionally shows a stage-level text progress bar.

## Examples

``` r
gp_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
gp <- read_genepred(gp_file, format = "genePredExt", verbose = FALSE)
gp
#> <Feature>
#>   records    : 284
#>   genes      : 20
#>   transcripts: 24
#>   exons      : 76
#>   format     : genePredExt
#>   coordinate : 1-based closed
head(gp$genes)
#>    gene_id  chrom strand gene_start gene_end n_transcripts gene_type
#>     <char> <char> <char>      <int>    <int>         <int>    <char>
#> 1:   GeneD   chr1      +    1000001  1008000             1    coding
#> 2:   GeneE   chr1      +    3000001  3009000             1    coding
#> 3:   GeneF   chr1      -    5000001  5010000             1    coding
#> 4:   GeneG   chr1      +    7000001  7007000             1    coding
#> 5:   GeneH   chr1      -    9000001  9009000             1    coding
#> 6:   GeneA   chr1      +   12340001 12352000             2    coding
head(gp$transcripts)
#>    transcript_id gene_id  chrom strand tx_start   tx_end cds_start  cds_end
#>           <char>  <char> <char> <char>    <int>    <int>     <int>    <int>
#> 1:          TxD1   GeneD   chr1      +  1000001  1008000   1000201  1007800
#> 2:          TxE1   GeneE   chr1      +  3000001  3009000   3000201  3008800
#> 3:          TxF1   GeneF   chr1      -  5000001  5010000   5000201  5009800
#> 4:          TxG1   GeneG   chr1      +  7000001  7007000   7000201  7006800
#> 5:          TxH1   GeneH   chr1      -  9000001  9009000   9000201  9008800
#> 6:          TxA1   GeneA   chr1      + 12340001 12352000  12340501 12351500
#>    exon_count gene_type score cds_start_stat cds_end_stat  exon_frames
#>         <int>    <char> <num>         <char>       <char>       <char>
#> 1:          3    coding     0       complete     complete    -1,-1,-1,
#> 2:          3    coding     0       complete     complete    -1,-1,-1,
#> 3:          3    coding     0       complete     complete    -1,-1,-1,
#> 4:          3    coding     0       complete     complete    -1,-1,-1,
#> 5:          3    coding     0       complete     complete    -1,-1,-1,
#> 6:          4    coding     0       complete     complete -1,-1,-1,-1,
```
