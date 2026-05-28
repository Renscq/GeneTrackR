# GeneTrackR

GeneTrackR is a lightweight R package for reading GenePred annotations and genomic signal tracks, then drawing gene models, signal tracks, and genome-browser-like combined figures.

The package is designed for workflows that need to inspect local gene structures together with RNA-seq, Ribo-seq, bedGraph, wig, or bigWig signal tracks.

## Installation

```r
# From a local source directory
devtools::document()
devtools::install()
```

## Example data

GeneTrackR ships with small example files for testing and learning:

```r
library(GeneTrackR)

gp_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
bg_files <- system.file(
  "extdata",
  c("example_signal_A.bedgraph", "example_signal_B.bedgraph"),
  package = "GeneTrackR"
)
```

The example annotation contains four genes across two chromosomes, including coding and non-coding genes, positive and negative strand genes, and a multi-transcript gene.

## Read GenePred annotation

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = TRUE
)

gp
```

For large annotation files, `read_genepred()` reports progress by default. To suppress progress messages:

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE
)
```

## Summarize annotation

```r
summary_genepred(gp, level = "gene")
summary_genepred(gp, level = "transcript")
summary_genepred(gp, level = "exon")
```

You can also summarize a genomic interval:

```r
summary_genepred(
  gp,
  chrom = "chr1",
  start = 1,
  end = 2000,
  level = "gene"
)
```

## Plot gene models

### Plot a gene

```r
plot_gene(
  gp,
  gene_id = "GeneA",
  label_position = "axis",
  direction_mode = "end"
)
```

### Plot a transcript

```r
plot_transcript(
  gp,
  transcript_id = "TxA1",
  coordinate = "transcript"
)
```

### Plot a genomic region

```r
plot_region(
  gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  mode = "overlap",
  label_by = "gene"
)
```

### Highlight a region

`highlight` must be a data frame with at least `start` and `end` columns.

```r
plot_region(
  gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  highlight = data.frame(start = 250, end = 500)
)
```

## Read signal tracks

GeneTrackR supports bedGraph, wig, and bigWig inputs. For bigWig files, the package uses its local libBigWig backend.

```r
bg <- read_bwg(
  bg_files,
  format = "bedgraph",
  mode = "lazy"
)

bg
```

Query signal over a region:

```r
dt <- query_bwg(
  bg,
  chrom = "chr1",
  start = 101,
  end = 900
)

head(dt)
```

Select specific samples:

```r
dt_a <- query_bwg(
  bg,
  chrom = "chr1",
  start = 101,
  end = 900,
  samples = "example_signal_A"
)
```

## Plot signal tracks

### Signal over a gene

```r
plot_signal_gene(
  bg,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "bar",
  signal_palette = "Blues",
  signal_y_scale = "free",
  signal_y_ticks = "range"
)
```

### Signal over a transcript

```r
plot_signal_transcript(
  bg,
  annotation = gp,
  transcript_id = "TxA1",
  plot_type = "bar",
  show_gene_model = TRUE
)
```

### Signal over a region

```r
plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 101,
  end = 900,
  plot_type = "bar",
  signal_transform = "sqrt"
)
```

## Combined genome-browser-like tracks

`plot_tracks()` is the recommended high-level entry point. It supports three locator modes.

### By gene ID

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  signal_type = "bar",
  signal_palette = "Blues"
)
```

### By transcript ID

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  transcript_id = "TxA1",
  signal_type = "bar"
)
```

### By genomic region

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  chrom = "chr1",
  start = 1,
  end = 1200,
  signal_type = "bar"
)
```

## Testing

The package includes a small `testthat` suite based on the example data:

```r
devtools::test()
devtools::check()
```

## Notes on coordinates

GenePred input is usually UCSC-style 0-based half-open. GeneTrackR converts it internally to 1-based closed coordinates. bedGraph input is also converted from 0-based half-open to 1-based closed coordinates.

## Search annotation IDs

Use `search_gene()` and `search_transcript()` when you only remember part of an ID.

```r
search_gene(gp, "GeneA")
search_transcript(gp, "TxA")
```

## Signal sample groups and replicate summaries

For KO/WT or treatment/control designs, you can provide a sample-to-group mapping.
This can be used for group-level colors:

```r
groups <- c(
  example_signal_A = "TreatmentA",
  example_signal_B = "TreatmentB"
)

plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 101,
  end = 900,
  sample_groups = groups,
  signal_color_by = "group",
  plot_type = "bar"
)
```

Replicate summaries are also supported. When raw interval boundaries differ among samples, set `bin_size` first so that samples are summarized on comparable windows.

```r
plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 101,
  end = 900,
  sample_groups = groups,
  signal_summary = "mean",
  signal_color_by = "group",
  bin_size = 100,
  plot_type = "line"
)
```

## bedGraph.gz tabix support

For large bgzip-compressed bedGraph files, GeneTrackR can use tabix when a `.tbi` index is present and the `tabix` command is available on `PATH`.

```r
bg <- read_bwg(
  files = "sample.bedgraph.gz",
  format = "bedgraph",
  mode = "lazy",
  use_tabix = "auto"
)

query_bwg(
  bg,
  chrom = "chr1",
  start = 1000,
  end = 2000,
  verbose = TRUE
)
```

If no usable tabix index is detected, GeneTrackR falls back to full-file `fread()` for that file.

## Advanced parameters and common customization

This section summarizes parameters that are easy to miss when using GeneTrackR.
Most of them are also documented in `?GeneTrackR-advanced-parameters`.

### Custom gene-model colors

Gene structure plots use `color_palette` by default. You can override feature
colors with `fill_colors`. The safest approach is to use a named vector:

```r
plot_gene(
  gp,
  gene_id = "GeneA",
  fill_colors = c(
    CDS = "#33a02c",
    UTR = "#b2df8a",
    exon = "#fb9a99"
  ),
  border_color = "black"
)
```

For transcript and region plots the same rule applies:

```r
plot_transcript(
  gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  fill_colors = c(CDS = "#33a02c", UTR = "#b2df8a", exon = "#fb9a99")
)

plot_region(
  gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  fill_colors = c(CDS = "#33a02c", UTR = "#b2df8a", exon = "#fb9a99"),
  label_position = "feature",
  label_side = "above"
)
```

### Highlight intervals

`highlight` is a data frame with at least `start` and `end` columns:

```r
highlight_region <- data.frame(start = 250, end = 500)

plot_region(
  gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  highlight = highlight_region
)
```

For genomic plots, `start` and `end` are genomic coordinates. For
`coordinate = "transcript"`, they are spliced transcript coordinates.

### Sample selection and signal colors

Use `samples` to select a subset of samples. Use `signal_palette` for automatic
RColorBrewer colors, or `signal_colors` for manual colors:

```r
plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  samples = c("example_signal_A", "example_signal_B"),
  signal_colors = c(
    example_signal_A = "#2166AC",
    example_signal_B = "#B2182B"
  )
)
```

All RColorBrewer palettes are supported. When there are more samples than the
palette maximum, colors are interpolated automatically.

### Group-level coloring and replicate summaries

Use `sample_groups` together with `signal_color_by = "group"` to color samples by
experimental group. Use `signal_summary` to summarize replicates by group.
When raw interval boundaries differ between samples, set `bin_size` before
summarizing.

```r
groups <- c(
  example_signal_A = "WT",
  example_signal_B = "KO"
)

plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 1,
  end = 1200,
  sample_groups = groups,
  signal_color_by = "group",
  signal_summary = "mean",
  bin_size = 50,
  signal_palette = "Set1"
)
```

### Signal y-axis transformation and scale

Signal plots support y-axis transformations and per-sample or shared y-axis
ranges:

```r
plot_signal_gene(
  bg,
  annotation = gp,
  gene_id = "GeneA",
  signal_transform = "sqrt",
  signal_y_scale = "fixed",
  signal_y_ticks = "range"
)
```

Available transformations are `none`, `log2`, `log10`, and `sqrt`. `range` ticks
show integer axis limits only, whereas `pretty` uses default ggplot-style breaks.

### Strand policy in signal queries

Standard bigWig and wig files are unstranded. With the default
`strand_policy = "ignore_unstranded"`, querying `strand = "+"` or `"-"` still
returns unstranded signal. For explicitly stranded bedGraph files, use:

```r
query_bwg(
  bg,
  chrom = "chr1",
  start = 1,
  end = 1200,
  strand = "+",
  strand_policy = "strict"
)
```

### Tabix-indexed bedGraph files

For large `.bedgraph.gz` files, `read_bwg(use_tabix = "auto")` will use tabix if
both a `.tbi` index and the `tabix` command are available. Otherwise it falls
back to full-file reading.

```r
bg_tabix <- read_bwg(
  "sample.bedgraph.gz",
  format = "bedgraph",
  mode = "lazy",
  use_tabix = "auto"
)
```


### Strand-specific bedGraph pairs

When plus and minus strand signals are stored in separate bedGraph files, provide
one `strand` value for each file. Querying `strand = "+"` returns only plus-strand
tracks, while querying `strand = "-"` returns only minus-strand tracks.

```r
bg_strand <- read_bwg(
  files = c("sample_minus.bedgraph.gz", "sample_plus.bedgraph.gz"),
  format = "bedgraph",
  strand = c("-", "+"),
  mode = "lazy"
)

plus_signal <- query_bwg(
  bg_strand,
  chrom = "chr1",
  start = 1,
  end = 10000,
  strand = "+"
)
```

Use `strand = "ignore"` when you want to combine or display all strand files.
