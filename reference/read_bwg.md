# Read bedGraph, wig, or bigWig signal tracks

Read bedGraph, wig, or bigWig signal tracks

## Usage

``` r
read_bwg(
  files,
  format = c("auto", "bedgraph", "bigwig", "wig"),
  sample_names = NULL,
  strand = NULL,
  mode = c("lazy", "memory"),
  genome = NULL,
  check_chrom = TRUE,
  use_tabix = c("auto", "yes", "no"),
  tabix_empty_fallback = FALSE,
  verbose = TRUE
)
```

## Arguments

- files:

  Signal file paths.

- format:

  Input format. Use auto, bedgraph, bigwig, or wig.

- sample_names:

  Optional sample names. If NULL, names are inferred from file basenames
  after removing compression suffixes and signal-track suffixes such as
  `.bedgraph.gz`, `.bigwig`, `.bw`, or `.wig`.

- strand:

  Optional strand labels for files.

- mode:

  lazy stores file paths, memory reads data into R.

- genome:

  Optional genome label.

- check_chrom:

  Whether to check chromosome names for in-memory files.

- use_tabix:

  Whether to use tabix/Rsamtools for indexed bedGraph-like files in lazy
  mode. Accepts `"auto"`, `"yes"`, `"no"`, TRUE, or FALSE. `TRUE` is
  equivalent to `"yes"`, and `FALSE` is equivalent to `"no"`.

- tabix_empty_fallback:

  Logical. Whether an empty tabix result should be verified by a
  full-file fread query. Default FALSE for performance. Set TRUE only
  when you suspect the tabix index coordinate convention is incompatible
  with the queried region.

- verbose:

  Logical. Whether to print read/setup progress messages.

## Value

A BwgTrack object.

## Details

`format = "auto"` infers the input type from file extensions after
removing compression suffixes such as `.gz` or `.bgz`.
`sample_names = NULL` infers names from file basenames, removing
suffixes such as `.bedgraph.gz`, `.bw`, `.bigwig`, and `.wig`.

For bigWig and wig files, strand information is not stored in the file,
so the sample is treated as unstranded. For paired plus/minus bedGraph
files, provide `strand = c("+", "-")` and matching `sample_names` if
strand-specific filtering is required.

`use_tabix = "auto"` uses indexed querying only when a `.tbi` index and
an available backend are detected. GeneTrackR first checks the system
`tabix` command and then the R package `Rsamtools`. Otherwise, lazy
bedGraph queries fall back to full-file reading.

## Examples

``` r
if (FALSE) { # \dontrun{
rnaseq_files <- system.file(
  "extdata",
  c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph"),
  package = "GeneTrackR"
)
riboseq_files <- system.file(
  "extdata",
  c("gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph"),
  package = "GeneTrackR"
)

rnaseq <- read_bwg(
  rnaseq_files,
  format = "bedgraph",
  sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
  strand = c("+", "-"),
  mode = "memory"
)
riboseq <- read_bwg(
  riboseq_files,
  format = "bedgraph",
  sample_names = c("Ribo_seq_plus", "Ribo_seq_minus"),
  strand = c("+", "-"),
  mode = "memory"
)
} # }
```
