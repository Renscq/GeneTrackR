# Read a BED file as a FeatureTrack

Reads BED3-BED12 style files and converts BED coordinates from 0-based
half-open to 1-based closed intervals by default.

## Usage

``` r
read_bed(
  file,
  coordinate = c("bed", "granges"),
  name_col = 4L,
  type = "BED",
  verbose = TRUE,
  progress = interactive() && isTRUE(verbose)
)
```

## Arguments

- file:

  BED file path. Gzip-compressed files are supported by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).

- coordinate:

  Input coordinate system. `bed` means 0-based half-open; `granges`
  means 1-based closed.

- name_col:

  BED column used as feature name. Default 4.

- type:

  Feature type assigned to BED records.

- verbose:

  Whether to print progress messages.

- progress:

  Whether to show a stage progress bar.

## Value

A FeatureTrack object.

## Examples

``` r
bed_file <- system.file("extdata", "gtr_demo_features.bed", package = "GeneTrackR")
features <- read_bed(bed_file, verbose = FALSE, progress = FALSE)
features
#> <Feature>
#>   records    : 15
#>   genes      : 0
#>   transcripts: 0
#>   exons      : 0
#>   format     : BED
#>   coordinate : 1-based closed
plot_feature_track(features, chrom = "chr1", start = 12338201, end = 12374500)
```
