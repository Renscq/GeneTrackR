<!-- README.md is generated from README.qmd and docs/*.qmd. Edit the QMD sources, not README.md. -->

**GeneTrackR** is a lightweight R package for reading, querying, writing, and visualizing gene annotations, genomic signal tracks, variant tracks, haplotypes, and phenotype associations. It is designed for programmable IGV-like visualization in R while keeping direct access to the underlying annotation, coverage, variant, haplotype, and phenotype tables.

## Main features

- Read and standardize annotations from **GenePred**, **GenePredExt**, **GTF**, **GFF3**, and **BED**.
- Read and query **bedGraph**, **wig**, and **bigWig** signal tracks.
- Read, lazily query, retrieve, merge, plot, and write **VCF** variant tracks.
- Draw gene models by gene, transcript, or genomic region.
- Draw signal tracks as bar, line, area, or heatmap tracks.
- Combine gene models, signal tracks, feature tracks, and variant tracks with `plot_tracks()`.
- Extract haplotypes from genes, transcripts, or arbitrary genomic regions.
- Draw haplotype-variant diagrams with gene models, variant markers, connector lines, and genotype tables.
- Read phenotype tables, summarize missingness/type, and draw phenotype distributions.
- Draw haplotype-phenotype and single-variant phenotype association plots with pairwise statistical tests.
- Export annotations, signal tracks, and variants to standard formats.

## Installation

GeneTrackR depends on CRAN and Bioconductor packages. Install the dependencies first, then install GeneTrackR from GitHub.

```{r}
## CRAN dependencies
install.packages(c(
  "devtools",
  "data.table",
  "ggplot2",
  "patchwork",
  "rlang",
  "Rcpp",
  "RColorBrewer"
))

## Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomicRanges",
  "IRanges",
  "Rsamtools"
), ask = FALSE, update = FALSE)

## Install GeneTrackR from GitHub
devtools::install_github(
  "Renscq/GeneTrackR",
  dependencies = TRUE,
  build_vignettes = FALSE
)
```

Load the package:

```{r}
library(GeneTrackR)
```

For local development:

```bash
git clone https://github.com/Renscq/GeneTrackR.git
cd GeneTrackR
```

```{r}
devtools::document()
devtools::install()
```

## Built-in example files

GeneTrackR ships with one deterministic demo genome in `inst/extdata`. The dataset is compact but intentionally designed to exercise the major modules without external files.

```{r}
gp_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
gtf_file <- system.file("extdata", "gtr_demo.gtf", package = "GeneTrackR")
gff_file <- system.file("extdata", "gtr_demo.gff3", package = "GeneTrackR")
bed_file <- system.file("extdata", "gtr_demo_features.bed", package = "GeneTrackR")

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

vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
pheno_file <- system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR")
```

The deterministic demo genome contains 2 chromosomes, 20 genes, 24 transcripts, 36 samples, 50 designed variants, strand-specific RNA-seq/Ribo-seq tracks, four balanced GeneA haplotype groups, and phenotype traits with known positive and negative-control associations. All protein-coding demo transcripts use CDS lengths divisible by three. The Ribo-seq tracks contain moderately dense heterogeneous integer P-site-like counts with designed frame0/frame1/frame2 total-count proportions of approximately 80%/10%/10%. Frame 0 is broadly occupied and variable, while frame 1 and frame 2 use different subsets of codons with lower irregular counts. Initiation and termination frame-0 counts are approximately two times the internal frame-0 mean.

## Annotation files: read, subset, merge, and write

GeneTrackR supports five annotation inputs: **GenePred**, **GenePredExt**, **GTF**, **GFF3**, and **BED**. All readers standardize coordinates internally to **1-based closed** coordinates, but the source formats are not equivalent: GenePred/GenePredExt, GTF, and GFF3 can provide hierarchical gene models, whereas BED is primarily treated as an interval annotation track.

### Supported formats and object types

| Format | Typical content | Reader | Main return | Writer format |
|---|---|---|---|---|
| GenePred | 10-column transcript model | `read_genepred(..., format = "genePred")` | `GenePred` | `genepred` |
| GenePredExt | 15-column GenePred with gene ID/status/frame fields | `read_genepred(..., format = "genePredExt")` | `GenePred` | `genepredext` |
| GTF | 9-column hierarchical annotation with attributes | `read_gtf()` | `FeatureTrack` with derived gene model | `gtf` |
| GFF3 | 9-column hierarchical annotation with `ID`/`Parent` relationships | `read_gff()` | `FeatureTrack` with derived gene model | `gff` |
| BED | BED3-BED12 interval annotation | `read_bed()` | interval-style `FeatureTrack` | `bed6`; `bed12` for gene-model objects |

Coordinate conventions are handled at the file boundary:

- GenePred/GenePredExt default to UCSC-style **0-based half-open** input with `coordinate = "ucsc"`.
- BED defaults to **0-based half-open** input with `coordinate = "bed"`.
- GTF and GFF3 use **1-based closed** coordinates.
- GeneTrackR stores all annotation coordinates internally as **1-based closed** coordinates.
- `write_feature()` converts GenePred/BED-like outputs back to UCSC coordinates by default with `coordinate = "ucsc"`.

The built-in annotation files describe the same demo genome:

```{r}
gp_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
gtf_file <- system.file("extdata", "gtr_demo.gtf", package = "GeneTrackR")
gff_file <- system.file("extdata", "gtr_demo.gff3", package = "GeneTrackR")
bed_file <- system.file("extdata", "gtr_demo_features.bed", package = "GeneTrackR")
```

### GenePred and GenePredExt

GenePred stores one transcript per row with transcript/CDS boundaries and exon blocks. Standard GenePred contains 10 core columns. GenePredExt adds `score`, `name2`, CDS status, and exon-frame information; GeneTrackR normally uses `name` as `transcript_id` and `name2` as `gene_id` for GenePredExt.

#### Read

The built-in file is GenePredExt:

```{r}
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE,
  progress = FALSE
)

gp
head(gp$genes)
head(gp$transcripts)
head(gp$exons)
```

Standard GenePred uses the same reader. Here a standard GenePred file is created from the demo object and read back:

```{r}
gp_standard_file <- file.path(tempdir(), "gtr_demo.genePred")

write_feature(
  gp,
  gp_standard_file,
  format = "genepred",
  overwrite = TRUE
)

gp_standard <- read_genepred(
  gp_standard_file,
  format = "genePred",
  verbose = FALSE,
  progress = FALSE
)
```

#### Subset

GenePred-compatible objects can be subset by gene, transcript, or genomic region:

```{r}
gp_gene_a <- retrieve_feature(
  gp,
  gene_id = "GeneA"
)

gp_tx_a1 <- retrieve_feature(
  gp,
  transcript_id = "TxA1"
)

gp_region <- retrieve_feature(
  gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  mode = "overlap"
)
```

For a plain table instead of a sub-object:

```{r}
gp_gene_table <- retrieve_feature(
  gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  level = "gene",
  as = "data.table"
)
```

#### Merge

Subsets from the same or different GenePred objects can be recombined with `merge_feature()`:

```{r}
gp_gene_b <- retrieve_feature(gp, gene_id = "GeneB")

gp_merged <- merge_feature(
  gp_gene_a,
  gp_gene_b,
  source_names = c("GeneA", "GeneB"),
  conflict = "deduplicate"
)
```

#### Write

A GenePred-compatible object can be written as either GenePred or GenePredExt:

```{r}
write_feature(
  gp_gene_a,
  file.path(tempdir(), "GeneA.genePred"),
  format = "genepred",
  overwrite = TRUE
)

write_feature(
  gp_gene_a,
  file.path(tempdir(), "GeneA.genePredExt"),
  format = "genepredext",
  overwrite = TRUE
)
```

`write_genepred()` remains available as a backward-compatible wrapper around `write_feature()`.

### GTF

GTF stores gene-model records in a 9-column table and uses the final attribute column for identifiers such as `gene_id` and `transcript_id`. GeneTrackR parses the hierarchical records into a `FeatureTrack` with standardized `$genes`, `$transcripts`, `$exons`, and `$data` tables when the required records are present.

#### Read

```{r}
gtf <- read_gtf(
  gtf_file,
  verbose = FALSE,
  progress = FALSE
)

gtf
head(gtf$genes)
head(gtf$transcripts)
```

Selected feature types can be loaded when the full annotation is not required:

```{r}
gtf_core <- read_gtf(
  gtf_file,
  feature_types = c("gene", "transcript", "exon", "CDS"),
  verbose = FALSE,
  progress = FALSE
)
```

#### Subset

```{r}
gtf_gene_a <- retrieve_feature(gtf, gene_id = "GeneA")

gtf_gene_a_exons <- retrieve_feature(
  gtf,
  gene_id = "GeneA",
  level = "exon",
  as = "data.table"
)

gtf_region <- retrieve_feature(
  gtf,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  mode = "overlap"
)
```

#### Merge

```{r}
gtf_gene_b <- retrieve_feature(gtf, gene_id = "GeneB")

gtf_merged <- merge_feature(
  gtf_gene_a,
  gtf_gene_b,
  source_names = c("GeneA", "GeneB")
)
```

#### Write

```{r}
write_feature(
  gtf_gene_a,
  file.path(tempdir(), "GeneA.gtf"),
  format = "gtf",
  overwrite = TRUE
)
```

Because a gene-model GTF contains transcript and exon hierarchy, it can also be converted directly to GenePred/GenePredExt or BED12 with `write_feature()`.

### GFF3

GFF3 also uses 9 columns, but hierarchy is represented primarily through `ID` and `Parent` attributes. `read_gff()` reconstructs gene/transcript/exon relationships and returns the same unified `FeatureTrack` interface used for GTF.

#### Read

```{r}
gff <- read_gff(
  gff_file,
  verbose = FALSE,
  progress = FALSE
)

gff
head(gff$genes)
head(gff$transcripts)
```

Feature-type filtering is also supported:

```{r}
gff_core <- read_gff(
  gff_file,
  feature_types = c("gene", "mRNA", "exon", "CDS"),
  verbose = FALSE,
  progress = FALSE
)
```

#### Subset

```{r}
gff_gene_a <- retrieve_feature(gff, gene_id = "GeneA")

gff_region <- retrieve_feature(
  gff,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  mode = "overlap"
)
```

#### Merge

```{r}
gff_gene_b <- retrieve_feature(gff, gene_id = "GeneB")

gff_merged <- merge_feature(
  gff_gene_a,
  gff_gene_b,
  source_names = c("GeneA", "GeneB")
)
```

#### Write

```{r}
write_feature(
  gff_gene_a,
  file.path(tempdir(), "GeneA.gff3"),
  format = "gff",
  overwrite = TRUE
)
```

### BED

BED is treated differently from GenePred/GTF/GFF3. `read_bed()` reads BED3-BED12-style files as **interval annotations** and standardizes the interval coordinates, name, score, and strand fields into a `FeatureTrack`. It does not reconstruct BED12 block columns into a transcript/exon hierarchy; use GenePred, GTF, or GFF3 when exon-level gene models are required.

#### Read

```{r}
features <- read_bed(
  bed_file,
  verbose = FALSE,
  progress = FALSE
)

features
head(features$data)
```

The demo BED contains promoters, enhancers, candidate regions, repeats, conserved regions, and QTL intervals. Their semantic type is encoded in the BED name field, for example `GeneA_promoter|promoter`.

#### Subset

BED intervals can be selected by region or text pattern:

```{r}
bed_gene_a <- retrieve_feature(
  features,
  chrom = "chr1",
  start = 12338201,
  end = 12352000,
  mode = "overlap"
)

bed_promoters <- retrieve_feature(
  features,
  pattern = "promoter"
)
```

#### Merge

```{r}
bed_gene_b <- retrieve_feature(
  features,
  chrom = "chr1",
  start = 12355501,
  end = 12368500,
  mode = "overlap"
)

bed_merged <- merge_feature(
  bed_gene_a,
  bed_gene_b,
  source_names = c("GeneA_region", "GeneB_region"),
  conflict = "keep_all"
)
```

#### Write

Generic BED-derived intervals should normally be written as BED6:

```{r}
write_feature(
  bed_merged,
  file.path(tempdir(), "merged_features.bed6"),
  format = "bed6",
  overwrite = TRUE
)
```

BED12 output requires transcript/exon gene-model information. It is therefore appropriate for `GenePred` or gene-model GTF/GFF3 objects rather than a generic BED interval track:

```{r}
write_feature(
  gp_gene_a,
  file.path(tempdir(), "GeneA.bed12"),
  format = "bed12",
  overwrite = TRUE
)
```

### Region-selection modes

`retrieve_feature()` uses the same region-selection interface across annotation formats:

| `mode` | Meaning |
|---|---|
| `"overlap"` | Keep records overlapping the requested region. This is the default. |
| `"within"` | Keep records fully contained in the requested region. |
| `"trim"` | For hierarchical gene models, keep overlapping records and clip transcript/exon boundaries to the requested interval. |

For generic BED-style interval tracks, `"trim"` behaves as overlap selection because there is no transcript/exon hierarchy to rebuild.

### Merge conflict handling

`merge_feature()` accepts `Feature`, `FeatureTrack`, and `GenePred` objects, including a list of objects. Duplicate identifiers are handled explicitly:

| `conflict` | Behavior |
|---|---|
| `"deduplicate"` | Keep the first compatible copy according to input order. Default. |
| `"rename"` | Rename conflicting IDs in later inputs and update hierarchy references. |
| `"keep_all"` | Retain conflicting records unchanged. |
| `"error"` | Stop when duplicate IDs are detected. |
| `"keep_first"` | Backward-compatible alias of `"deduplicate"`. |

When the same locus is loaded from different annotation formats, `rename` is useful if both representations should be retained:

```{r}
cross_format <- merge_feature(
  list(gp_gene_a, gtf_gene_a),
  source_names = c("GenePred", "GTF"),
  conflict = "rename"
)
```

If the two inputs represent the same annotation and only one copy is needed, use `conflict = "deduplicate"` instead.

### Cross-format conversion and standardized tables

The readers share a common internal representation, so annotation formats can be converted without reparsing text manually:

```{r}
# GTF/GFF3 gene models -> GenePred-compatible object
gtf_as_gp <- as_genepred(gtf)
gff_as_gp <- as_genepred(gff)

# Any annotation object -> standardized feature table
feature_dt <- as_feature_table(gp)

# Gene-model tables
gene_dt <- as_gene_table(gp)
tx_dt <- as_transcript_table(gp)
exon_dt <- as_exon_table(gp)

# Bioconductor interoperability
gene_gr <- as_granges(gp, level = "gene")
```

`write_feature()` is also the main cross-format writer. For example, a GTF-derived gene model can be exported as GenePredExt or BED12:

```{r}
write_feature(
  gtf_gene_a,
  file.path(tempdir(), "GeneA.from_gtf.genePredExt"),
  format = "genepredext",
  overwrite = TRUE
)

write_feature(
  gtf_gene_a,
  file.path(tempdir(), "GeneA.from_gtf.bed12"),
  format = "bed12",
  overwrite = TRUE
)
```

### Annotation object contract

After reading, GeneTrackR keeps annotation data in a small set of compatible objects:

| Input | Main object | Important contents |
|---|---|---|
| GenePred / GenePredExt | `GenePred` | `$genes`, `$transcripts`, `$exons`, `$data` |
| Gene-model GTF / GFF3 | `FeatureTrack` | `$data` plus derived `$genes`, `$transcripts`, `$exons` |
| BED | interval `FeatureTrack` | primarily `$data` |

`retrieve_feature()` returns a compatible sub-object by default (`as = "Feature"`) and a `data.table` when `as = "data.table"`. `merge_feature()` always returns a unified `Feature` object and additionally inherits from `GenePred` when transcript/exon tables are available.

### How annotation objects connect to the rest of GeneTrackR

Annotation objects are the first input to most downstream workflows. `GenePred` or gene-model `FeatureTrack` objects can be passed directly to gene-model plotting, signal plotting, haplotype extraction, LD visualization, and integrated browser-style tracks.

### Core object flow and return contracts

GeneTrackR uses a small set of S3 objects throughout the workflow. The examples below keep these objects intact instead of repeatedly converting them to plain tables.

| Step | Function | Main return | Important contents / next consumer |
|---|---|---|---|
| Annotation | `read_genepred()` | `GenePred` | `$genes`, `$transcripts`, `$exons`, `$data`; used by gene/signal/haplotype plotting |
| GTF/GFF/BED | `read_gtf()`, `read_gff()`, `read_bed()` | `FeatureTrack` | `$data` plus derived hierarchy when available; gene-model FeatureTracks are GenePred-convertible |
| Signal | `read_bwg()` / `merge_bwg()` | `BwgTrack` | `$samples`, optional in-memory `$data`, `$meta`; consumed by `plot_signal_*()` and `plot_tracks()` |
| Variant | `read_vcf()` | `VariantTrack` | `$data`, `$meta`; consumed by haplotype, LD, phenotype, and track functions |
| Variant subset | `retrieve_vcf()` | `data.table` by default; `VariantTrack` with `as = "VariantTrack"` | use the table for inspection or request `VariantTrack` for downstream track/object workflows |
| Haplotype | `hap_gene_variant()` / `hap_region_variant()` | `HapVariant` | `$variants`, `$haplotypes`, `$sample_haplotypes`, `$genotype_wide` |
| Phenotype association | `plot_hap_pheno()` / `plot_variant_pheno()` | `GeneTrackRPhenoPlot` | `$figure`, `$pvalue`, `$summary`, `$plot_data` |
| LD | `compute_ld_block()` | `LDTrack` | `$data`, `$matrix`, `$variants`, `$region`; `plot_ld_block()` stores the figure in `$figure` |
| Refinement | `refine_haplotype()` | `HapRefined` | `$refined_hap`, `$refined_haplotypes`, `$haplotype_map`, `$pairwise_test` |
| Variant effect | `plot_variant_effect()` | `GeneTrackRVariantEffectPlot` | `$figure`, `$effect`, `$plot_data` |

Plotting functions such as `plot_gene()`, `plot_signal_gene()`, `plot_variant()`, `plot_tracks()`, `plot_hap_variant()`, and `plot_refined_hap_variant()` return a ggplot/patchwork figure directly. In contrast, phenotype/effect plotting functions return result objects containing `$figure` plus analysis tables. `plot_ld_block()` is intentionally stateful: by default it returns an updated `LDTrack` with the plot stored in `$figure`; use `return_object = FALSE` only when a figure-only return is required.

## Gene model plotting

### Plot a gene

```{r}
plot_gene(
  gp,
  gene_id = "GeneA",
  label_position = "axis",
  direction_mode = "end"
)
```

### Plot a transcript

```{r}
plot_transcript(
  gp,
  transcript_id = "TxA1",
  coordinate = "genomic"
)
```

Spliced transcript coordinate mode removes introns:

```{r}
plot_transcript(
  gp,
  transcript_id = "TxA1",
  coordinate = "transcript"
)
```

### Plot a genomic region

```{r}
plot_region(
  gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  mode = "overlap",
  label_by = "gene"
)
```

### Customize gene model colors

```{r}
plot_gene(
  gp,
  gene_id = "GeneA",
  gene_palette = "Set2",
  gene_border_color = "black"
)

plot_gene(
  gp,
  gene_id = "GeneA",
  gene_colors = c(
    CDS = "#1b9e77",
    UTR = "#a6d854",
    exon = "#7570b3"
  )
)
```

## RNA-seq and Ribo-seq signal tracks

GeneTrackR uses `BwgTrack` objects for continuous genomic signals. This module follows one complete workflow for the strand-specific RNA-seq and Ribo-seq demo data:

1. read bedGraph files;
2. write bigWig files;
3. plot gene-level signal tracks;
4. plot transcript-level signal tracks;
5. combine RNA-seq and Ribo-seq with `plot_tracks()`;
6. plot signals across a genomic region.

The examples use the same `GeneA`/`TxA1` positive-strand locus whenever possible so that the outputs can be compared directly between steps.

### Step 1. Read RNA-seq and Ribo-seq bedGraph data

Load the annotation and locate the four strand-specific signal files:

```{r}
gp_file <- system.file(
  "extdata",
  "gtr_demo.genePredExt",
  package = "GeneTrackR"
)

rnaseq_files <- system.file(
  "extdata",
  c(
    "gtr_demo_rnaseq_plus.bedgraph",
    "gtr_demo_rnaseq_minus.bedgraph"
  ),
  package = "GeneTrackR"
)

riboseq_files <- system.file(
  "extdata",
  c(
    "gtr_demo_riboseq_plus.bedgraph",
    "gtr_demo_riboseq_minus.bedgraph"
  ),
  package = "GeneTrackR"
)

gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE,
  progress = FALSE
)
```

Read the RNA-seq plus/minus bedGraph files into one `BwgTrack` object:

```{r}
rnaseq <- read_bwg(
  rnaseq_files,
  format = "bedgraph",
  sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
  strand = c("+", "-"),
  mode = "memory"
)

rnaseq
summary_bwg(rnaseq)
```

Read the Ribo-seq plus/minus bedGraph files in the same way:

```{r}
riboseq <- read_bwg(
  riboseq_files,
  format = "bedgraph",
  sample_names = c("Ribo_seq_plus", "Ribo_seq_minus"),
  strand = c("+", "-"),
  mode = "memory"
)

riboseq
summary_bwg(riboseq)
```

The two assays intentionally have different signal structures:

- **RNA-seq** signal is enriched across exons, including UTRs, while intronic and intergenic regions have little or no coverage in the demo data.
- **Ribo-seq** signal uses moderately dense 1-bp integer P-site-like counts within CDS regions. Frame 0 is broadly occupied with variable heights; frame 1 and frame 2 occur at different subsets of codons with lower irregular counts. Zero-count bases are omitted from bedGraph. Total counts remain approximately 80%/10%/10% for frame 0/frame 1/frame 2, and the initiation/termination frame-0 counts are approximately two times the internal frame-0 mean.

Because bedGraph does not store strand metadata, `strand = c("+", "-")` explicitly records the strand associated with each input file in the `BwgTrack` sample table.

### Step 2. Write RNA-seq and Ribo-seq bigWig files

`write_bwg()` can convert an in-memory `BwgTrack` to bigWig directly with the bundled third-party libBigWig library in `src/`. No external conversion program is required; only chromosome sizes are needed for the bigWig header.

```{r}
chrom_sizes_file <- system.file(
  "extdata",
  "gtr_demo.chrom.sizes",
  package = "GeneTrackR"
)

bigwig_dir <- file.path(tempdir(), "GeneTrackR_demo_bigwig")
dir.create(bigwig_dir, recursive = TRUE, showWarnings = FALSE)
```

Write the RNA-seq tracks:

```{r}
rnaseq_bigwig <- write_bwg(
  rnaseq,
  outdir = bigwig_dir,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)

rnaseq_bigwig
```

Write the Ribo-seq tracks:

```{r}
riboseq_bigwig <- write_bwg(
  riboseq,
  outdir = bigwig_dir,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)

riboseq_bigwig
```

Each call invisibly returns a table with `sample_id`, output `file`, and `format`. With the sample names used above, the output directory contains:

```text
RNA_seq_plus.bigwig
RNA_seq_minus.bigwig
Ribo_seq_plus.bigwig
Ribo_seq_minus.bigwig
```

BigWig export has a single backend in GeneTrackR: the bundled libBigWig implementation. This keeps the write path deterministic across platforms and avoids an external executable dependency.

### Step 3. Plot RNA-seq and Ribo-seq tracks for a gene

`plot_signal_gene()` retrieves the genomic span of a gene and optionally adds the gene model below the signal panel. With `strand = "auto"`, the gene strand is used to select the matching signal track. `GeneA` is on the positive strand, so the following examples use the plus RNA-seq and Ribo-seq samples automatically.

RNA-seq gene track:

```{r}
p_rnaseq_gene <- plot_signal_gene(
  signal = rnaseq,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "bar",
  strand = "auto",
  signal_palette = "Blues",
  signal_palette_direction = -1,
  signal_y_scale = "fixed",
  signal_y_ticks = "pretty",
  signal_track_height = 3,
  gene_track_height = 1
)

p_rnaseq_gene
```

Ribo-seq gene track:

```{r}
p_riboseq_gene <- plot_signal_gene(
  signal = riboseq,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "bar",
  strand = "auto",
  signal_palette = "Reds",
  signal_palette_direction = -1,
  signal_y_scale = "fixed",
  signal_y_ticks = "pretty",
  signal_track_height = 3,
  gene_track_height = 1
)

p_riboseq_gene
```

At a whole-gene genomic scale, thousands of bases are compressed into the available plot width, so a three-nucleotide pattern cannot be resolved visually even when the underlying RPF counts are frame-biased. The gene-level bar track is intended to show where translation signal occurs across CDS exons. Use the transcript-level `frame` view below, or a short genomic window, to inspect three-nucleotide periodicity.

The relative vertical space occupied by the signal and gene-model panels is controlled directly by `signal_track_height` and `gene_track_height`.

### Step 4. Plot RNA-seq and Ribo-seq tracks for a transcript

`plot_signal_transcript()` focuses on a single transcript. The RNA-seq example uses the standard bar representation in transcript coordinates:

```{r}
p_rnaseq_transcript <- plot_signal_transcript(
  signal = rnaseq,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  plot_type = "bar",
  strand = "auto",
  signal_palette = "Blues",
  signal_palette_direction = -1,
  signal_y_scale = "fixed",
  signal_track_height = 3,
  gene_track_height = 1
)

p_rnaseq_transcript
```

For Ribo-seq, `plot_type = "frame"` maps P-site counts back onto the transcript reading frame. The demo is designed so that frame 0 contributes about 80% of total RPF counts, while frame 1 and frame 2 each contribute about 10%. Counts within each frame are deliberately heterogeneous rather than repeated at nearly identical heights. Discrete palettes such as `Set1` are assigned strictly in palette order: `frame0`, `frame1`, and `frame2` receive the first, second, and third palette colors, respectively:

```{r}
p_riboseq_transcript <- plot_signal_transcript(
  signal = riboseq,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  plot_type = "frame",
  strand = "auto",
  frame_palette = "Set1",
  signal_track_height = 3,
  gene_track_height = 1
)

p_riboseq_transcript
```

Use `plot_type = "bar"` instead when the goal is to compare the genomic distribution of Ribo-seq P-site counts with the RNA-seq coverage representation:

```{r}
plot_signal_transcript(
  signal = riboseq,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  plot_type = "bar",
  strand = "auto",
  signal_palette = "Reds",
  signal_palette_direction = -1
)
```

### Step 5. Plot RNA-seq and Ribo-seq together with `plot_tracks()`

Merge the two `BwgTrack` objects so RNA-seq and Ribo-seq can be displayed in one integrated track figure:

```{r}
signal_all <- merge_bwg(rnaseq, riboseq)

signal_all
summary_bwg(signal_all)
```

For the positive-strand `GeneA` locus, explicitly select the two plus-strand samples:

```{r}
p_signal_tracks <- plot_tracks(
  annotation = gp,
  signal = signal_all,
  gene_id = "GeneA",
  samples = c("RNA_seq_plus", "Ribo_seq_plus"),
  strand = "+",
  signal_type = "bar",
  signal_palette = "Dark2",
  signal_palette_direction = 1,
  signal_y_scale = "free",
  signal_y_ticks = "pretty",
  heights = c(
    signal = 4,
    gene = 1,
    feature = 0.8,
    variant = 0.7
  )
)

p_signal_tracks
```

`plot_tracks()` uses its named `heights` vector to control panel proportions. The `signal` and `gene` entries play the same role as `signal_track_height` and `gene_track_height` in the dedicated signal plotting functions.

### Step 6. Plot RNA-seq and Ribo-seq tracks across a genomic region

`plot_signal_region()` is useful when the region contains several genes or when the target interval is not defined by a single gene/transcript.

For `bar` and `line` region plots, sample/group colors follow the sample/group level order and the standard RColorBrewer class order. Discrete colors are not selected by interpolating between the first and last palette colors. Heatmaps remain continuous gradients.

The following region contains the positive-strand `GeneA`, negative-strand `GeneB`, and positive-strand `GeneC`, making it useful for displaying both strand-specific RNA-seq tracks:

```{r}
p_rnaseq_region <- plot_signal_region(
  signal = rnaseq,
  annotation = gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  strand = "both",
  plot_type = "bar",
  signal_palette = "Blues",
  signal_palette_direction = -1,
  signal_y_scale = "free",
  signal_y_ticks = "pretty",
  signal_track_height = 4,
  gene_track_height = 1
)

p_rnaseq_region
```

Plot the Ribo-seq tracks over the same interval:

```{r}
p_riboseq_region <- plot_signal_region(
  signal = riboseq,
  annotation = gp,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  strand = "both",
  plot_type = "bar",
  signal_palette = "Reds",
  signal_palette_direction = -1,
  signal_y_scale = "free",
  signal_y_ticks = "pretty",
  signal_track_height = 4,
  gene_track_height = 1
)

p_riboseq_region
```

For large real-world signal files, use lazy `BwgTrack` access where appropriate so that only the requested genomic interval is retrieved during plotting.

## VCF variant tracks

GeneTrackR stores VCF records in a `VariantTrack` object. This module follows one continuous workflow for the demo VCF:

1. read VCF data and inspect the `VariantTrack` object;
2. validate and summarize variants;
3. retrieve variants from a genomic region;
4. retrieve variants by gene or transcript;
5. filter variants by ID, type, or text pattern;
6. plot variant tracks.

The examples use the same `GeneA` / `TxA1` locus used elsewhere in the documentation so that variant, haplotype, LD, phenotype, and signal examples can be compared directly.

### Step 1. Read VCF data

Locate the demo VCF and annotation files:

```{r}
vcf_file <- system.file(
  "extdata",
  "gtr_demo_variants.vcf",
  package = "GeneTrackR"
)

gp_file <- system.file(
  "extdata",
  "gtr_demo.genePredExt",
  package = "GeneTrackR"
)

gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE,
  progress = FALSE
)
```

Read the complete demo VCF into memory:

```{r}
vcf <- read_vcf(
  vcf_file,
  mode = "memory",
  keep_genotype = TRUE,
  verbose = FALSE,
  progress = FALSE
)

vcf
```

`read_vcf()` returns a `VariantTrack`. Standard VCF fields are normalized internally to columns such as `chrom`, `pos`, `variant_id`, `ref`, `alt`, `qual`, `filter`, `info`, and `variant_type`. Genotype columns are retained when `keep_genotype = TRUE`.

```{r}
head(vcf$data[, c(
  "chrom",
  "pos",
  "variant_id",
  "ref",
  "alt",
  "variant_type"
)])

vcf$meta$sample_names
```

Coordinates inside `VariantTrack` use **1-based genomic positions**, consistent with the VCF specification.

#### Large indexed VCF files

For a bgzip-compressed VCF with a `.tbi` or `.csi` index, `mode = "auto"` or `mode = "lazy"` avoids loading all records into memory. A lazy `VariantTrack` stores the source file and sample metadata first; variants are loaded only when a genomic region is requested.

```{r}
# large_vcf <- "/path/to/large.vcf.gz"
#
# vcf_lazy <- read_vcf(
#   large_vcf,
#   mode = "lazy"
# )
#
# vcf_lazy
#
# vcf_lazy_region <- retrieve_vcf(
#   vcf_lazy,
#   chrom = "chr1",
#   start = 12339001,
#   end = 12374500,
#   as = "VariantTrack"
# )
```

Lazy regional queries require an indexed VCF and `Rsamtools`. For small files such as the bundled demo VCF, memory mode is simpler.

### Step 2. Validate and summarize variants

Use `validate_vcf()` to check required fields, coordinates, allele fields, duplicated variant IDs, and consistency between `pos`, `start`, and `end`:

```{r}
vcf_validation <- validate_vcf(vcf)

vcf_validation$invalid_summary
vcf_validation$warnings
```

The validation result contains:

- `invalid_records`: individual records that failed a check;
- `invalid_summary`: number of invalid records by reason;
- `warnings`: non-record-specific validation messages.

Use `summary_vcf()` to summarize the current variant collection. A genomic region is **not required** for an in-memory `VariantTrack`; with no range supplied, the complete object is summarized. By default, variants are counted by chromosome and inferred variant type:

```{r}
vcf_summary <- summary_vcf(vcf)
vcf_summary
```

A chromosome or complete genomic range can still be supplied when only a subset should be summarized:

```{r}
chr1_summary <- summary_vcf(vcf, chrom = "chr1")
chr1_summary
```

The demo VCF intentionally contains SNPs and indels so the same file can be reused for variant visualization, haplotype analysis, and LD analysis.

### Step 3. Retrieve variants from a genomic region

`retrieve_vcf()` can return either a plain `data.table` or a new `VariantTrack`.

Retrieve the broader `GeneA` / `GeneB` / `GeneC` demonstration region as a table:

```{r}
region_table <- retrieve_vcf(
  vcf,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  as = "data.table",
  verbose = FALSE
)

region_table[, .(
  chrom,
  pos,
  variant_id,
  ref,
  alt,
  variant_type
)]
```

Use `as = "VariantTrack"` when the subset will be passed to another GeneTrackR function:

```{r}
region_variants <- retrieve_vcf(
  vcf,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  as = "VariantTrack",
  verbose = FALSE
)

region_variants
```

The distinction is important:

- `as = "data.table"` is convenient for inspection, filtering, and tabular analysis;
- `as = "VariantTrack"` preserves the GeneTrackR object interface for downstream plotting and analysis.

### Step 4. Retrieve variants by gene or transcript

When annotation is available, `retrieve_vcf()` can resolve genomic coordinates directly from a gene or transcript ID.

Retrieve variants inside the `GeneA` gene body:

```{r}
genea_variants <- retrieve_vcf(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  as = "VariantTrack",
  verbose = FALSE
)

genea_variants
```

The canonical demo design contains **11 variants within the GeneA gene body**.

Upstream and downstream flanking regions can be added explicitly:

```{r}
genea_extended <- retrieve_vcf(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 500,
  strand_aware = TRUE,
  as = "VariantTrack",
  verbose = FALSE
)

genea_extended
```

For `GeneA`, the upstream extension additionally includes `varAup01`, which is intentionally designed with missing and heterozygous genotypes for edge-case testing.

Transcript-aware retrieval uses the same interface:

```{r}
txa1_variants <- retrieve_vcf(
  vcf,
  annotation = gp,
  transcript_id = "TxA1",
  as = "VariantTrack",
  verbose = FALSE
)

txa1_variants
```

For negative-strand genes or transcripts, `strand_aware = TRUE` makes `upstream` and `downstream` follow transcriptional direction rather than simply subtracting or adding genomic coordinates.

### Step 5. Filter variants by ID, type, or pattern

These filters do not require a genomic region for an in-memory `VariantTrack`; they can be used independently or combined with region/gene/transcript filters.

Retrieve one or more known variant IDs:

```{r}
selected_variants <- retrieve_vcf(
  vcf,
  variant_id = c("varA03", "varA04", "varA05"),
  as = "data.table",
  verbose = FALSE
)

selected_variants[, .(
  chrom,
  pos,
  variant_id,
  ref,
  alt,
  variant_type
)]
```

Filter by inferred variant type:

```{r}
genea_indels <- retrieve_vcf(
  genea_variants,
  variant_type = c("INS", "DEL"),
  as = "data.table",
  verbose = FALSE
)

genea_indels[, .(
  pos,
  variant_id,
  ref,
  alt,
  variant_type
)]
```

`pattern` searches `variant_id`, `REF`, `ALT`, `INFO`, and `variant_type`. For example, the demo high-LD variants contain `high_ld` in the INFO role field:

```{r}
high_ld_variants <- retrieve_vcf(
  vcf,
  pattern = "high_ld",
  fixed = TRUE,
  as = "VariantTrack",
  verbose = FALSE
)

high_ld_variants
```

These filters can be combined with genomic, gene, or transcript queries when a more specific subset is needed.

### Step 6. Plot variant tracks

`plot_variant()` returns a ggplot object directly. The GeneA region contains SNP, insertion, and deletion examples, so it is useful for demonstrating variant-type colors.

```{r}
p_genea_variants <- plot_variant(
  genea_variants,
  color_by = "variant_type",
  label_by = "variant_id",
  variant_shape = "lollipop",
  variant_palette = "Set1",
  point_size = 2.5,
  text_size = 12
)

p_genea_variants
```

For a broader genomic region, use the previously retrieved `VariantTrack`:

```{r}
p_region_variants <- plot_variant(
  region_variants,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  color_by = "variant_type",
  label_by = "none",
  variant_shape = "rug",
  variant_palette = "Set1",
  text_size = 12
)

p_region_variants
```

The main plotting controls are:

- `color_by`: `"variant_type"`, `"filter"`, or `"none"`;
- `label_by`: variant ID, type, REF, ALT, or no labels;
- `variant_shape`: `"lollipop"`, `"point"`, or `"rug"`;
- `variant_palette` / `variant_colors`: automatic palette or explicit colors.

Variant subsets created here are reused directly by the later haplotype, phenotype, LD, and integrated-track modules.

VCF export is covered separately in the **Export variants and analysis results** module. The current `write_vcf()` interface writes standardized site-level VCF fields; the original genotype-containing VCF should be retained as the archival source when genotype round-trip preservation is required.

## Browser-like combined tracks

`plot_tracks()` combines gene models, signal tracks, BED feature tracks, and VCF variant tracks.
It returns the assembled ggplot/patchwork figure directly rather than wrapping it in a result list.

```{r}
plot_tracks(
  annotation = gp,
  signal = rnaseq,
  gene_id = "GeneA",
  signal_type = "bar"
)
```

Customize gene model colors in combined tracks:

```{r}
plot_tracks(
  annotation = gp,
  signal = rnaseq,
  gene_id = "GeneA",
  signal_type = "bar",
  gene_palette = "Set2",
  gene_border_color = "black"
)

plot_tracks(
  annotation = gp,
  signal = rnaseq,
  gene_id = "GeneA",
  gene_colors = c(
    CDS = "#1b9e77",
    UTR = "#a6d854",
    exon = "#7570b3"
  )
)
```

Control the complete track style with standardized parameters:

```{r}
plot_tracks(
  annotation = gp,
  signal = rnaseq,
  gene_id = "GeneA",
  signal_type = "bar",
  signal_transform = "sqrt",
  signal_y_scale = "fixed",
  signal_y_limits = c(0, 20),
  signal_alpha = 0.80,
  signal_bar_width = 0.85,
  plot_theme = "classic",
  show_panel_border = FALSE
)
```

Add feature and variant tracks:

```{r}
plot_tracks(
  annotation = gp,
  signal = signal_all,
  features = features,
  variants = vcf,
  chrom = "chr1",
  start = 12339001,
  end = 12374500,
  signal_type = "bar"
)
```

## Haplotype extraction

GeneTrackR now separates gene/transcript-based and region-based haplotype extraction.

### Gene or transcript haplotypes

```{r}
hap_gene <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "string",
  min_variant_number = 1
)

hap_tx <- hap_gene_variant(
  vcf,
  annotation = gp,
  transcript_id = "TxA1",
  genotype_mode = "code",
  min_variant_number = 1
)
```

Include upstream/downstream variants around a gene or transcript:

```{r}
hap_gene_ext <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  strand_aware = TRUE,
  genotype_mode = "string",
  min_variant_number = 1
)
```

### Region haplotypes

```{r}
hap_region <- hap_region_variant(
  vcf,
  chrom = "chr1",
  start = 12339700,
  end = 12352000,
  genotype_mode = "code",
  min_variant_number = 1
)
```

`hap_variant()` is retained as a compatibility wrapper, but new code should prefer `hap_gene_variant()` and `hap_region_variant()`.

### Inspect haplotype tables

```{r}
hap_gene$region
hap_gene$variants
hap_gene$haplotypes
hap_gene$sample_haplotypes
```

`hap_gene` is a `HapVariant`. Downstream haplotype plots and phenotype/refinement functions consume this object directly; there is no need to reconstruct genotype tables manually.

## Haplotype-variant plot

`plot_hap_variant()` draws a gene model, variant markers, connector lines, and a genotype table.

```{r}
hap_variant_figure <- plot_hap_variant(
  hap_gene,
  annotation = gp,
  min_hap_samples = 3,
  show_reference_row = TRUE,
  table_x_angle = 90
)

hap_variant_figure
```

Customize table and variant colors:

```{r}
plot_hap_variant(
  hap_gene,
  annotation = gp,
  min_hap_samples = 3,
  table_palette = "RdBu",
  table_alpha = 0.6,
  variant_palette = "Set1",
  genotype_text_size = 3
)
```

## Phenotype input and summary

The phenotype table should contain sample/taxa IDs in the first column or a named sample column. Each additional column is treated as one trait.

```{r}
pheno <- read_pheno(pheno_file)

summary_pheno(pheno)

plot_pheno(
  pheno,
  traits = c("plant_height", "seed_weight")
)
```

`read_pheno()` returns a `data.table`; `plot_pheno()` returns a figure directly.

## Haplotype-phenotype association

`plot_hap_pheno()` compares phenotype distributions among haplotype groups. It returns both the figure and the p-value table.

```{r}
hap_res <- plot_hap_pheno(
  hap = hap_gene,
  phenotype = pheno,
  traits = "plant_height",
  min_hap_samples = 3,
  test_method = "t.test",
  p_value_type = "raw",
  p_label = "stars"
)

hap_res$figure
hap_res$pvalue
hap_res$summary
```

Use Wilcoxon or KS tests:

```{r}
plot_hap_pheno(
  hap = hap_gene,
  phenotype = pheno,
  traits = "seed_weight",
  test_method = "wilcox.test",
  min_hap_samples = 3
)$figure
```

Multiple traits preserve the requested order and are arranged over multiple facet columns:

```{r}
plot_hap_pheno(
  hap = hap_gene,
  phenotype = pheno,
  traits = c(
    "seed_weight",
    "plant_height",
    "protein_content",
    "flowering_time"
  ),
  facet_ncol = 2,
  min_hap_samples = 3
)$figure
```

## Single-variant phenotype association

`plot_variant_pheno()` is the single-variant version of `plot_hap_pheno()`. It groups samples by genotype or allele state at one variant.

```{r}
variant_res <- plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  variant_id = "varA03",
  traits = "protein_content",
  genotype_mode = "string",
  min_group_samples = 3
)

variant_res$figure
variant_res$pvalue
variant_res$variant_data
```

You can also select a variant by genomic position:

```{r}
plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  chrom = "chr1",
  pos = 12342550,
  traits = "protein_content",
  genotype_mode = "code",
  min_group_samples = 3
)$figure
```

## Linkage disequilibrium analysis

The demo VCF contains a designed six-variant high-LD block (`varLD01`-`varLD06`) in `GeneA`. These variants share the same genotype pattern and therefore provide a deterministic `r2 = 1` example.

```{r}
ld <- compute_ld_block(
  vcf,
  chrom = "chr1",
  start = 12342620,
  end = 12343180,
  method = "r2",
  verbose = FALSE
)

ld$data
```

Add the compact GeneA structure above the triangular heatmap:

```{r}
ld <- plot_ld_block(
  ld,
  show_region = TRUE,
  annotation = gp,
  show_variant_labels = FALSE
)

ld$figure
```

The assignment is intentional: `plot_ld_block()` returns the updated `LDTrack`, not the figure, when `return_object = TRUE` (the default).

The `GeneT` region contains exactly two variants and is retained as the deterministic two-variant LD plotting case:

```{r}
ld_pair <- compute_ld_block(
  vcf,
  chrom = "chr2",
  start = 16995001,
  end = 17006000,
  method = "r2",
  verbose = FALSE
)

plot_ld_block(ld_pair, return_object = FALSE)
```

## Haplotype refinement

`seed_weight` has a designed GeneA haplotype effect and can be used to demonstrate phenotype-guided haplotype refinement.

```{r}
refined <- refine_haplotype(
  hap_gene,
  phenotype = pheno,
  traits = "seed_weight",
  min_hap_samples = 3
)

refined$refined_haplotypes
```

`refined` is a `HapRefined`; `refined$refined_hap` is the corresponding refined `HapVariant` used internally by the refined plotting wrappers.

The refined object can be plotted with the same phenotype and variant interfaces:

```{r}
refined_pheno <- plot_refined_hap_pheno(
  refined,
  phenotype = pheno,
  traits = "seed_weight",
  min_hap_samples = 3
)
refined_pheno$figure
refined_pheno$pvalue

refined_variant_figure <- plot_refined_hap_variant(
  refined,
  annotation = gp,
  min_hap_samples = 3
)
refined_variant_figure
```

## Variant effect prioritization

`protein_content` was designed around `varA03`, providing a deterministic positive-effect example for `plot_variant_effect()`.

```{r}
effect_res <- plot_variant_effect(
  hap_gene,
  phenotype = pheno,
  traits = "protein_content",
  min_group_samples = 3,
  x_axis = "position"
)

effect_res$figure
effect_res$effect
```

## Export variants and analysis results

Annotation writing and cross-format conversion are covered in the **Annotation files** module for GenePred, GenePredExt, GTF, GFF3, BED6, and BED12.

### Export variants

```{r}
write_vcf(
  vcf,
  file = file.path(tempdir(), "gtr_demo.output.vcf"),
  overwrite = TRUE
)
```

### Save figures and tables

```{r}
ggplot2::ggsave(
  filename = file.path(tempdir(), "hap_pheno.pdf"),
  plot = hap_res$figure,
  width = 6,
  height = 5
)

data.table::fwrite(
  hap_res$pvalue,
  file = file.path(tempdir(), "hap_pheno.pvalue.tsv"),
  sep = "\t"
)
```

## Recommended workflow

```{r}
library(GeneTrackR)

## Input files
gp_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
pheno_file <- system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR")
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

## Read data
gp <- read_genepred(gp_file, format = "genePredExt", verbose = FALSE)
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
pheno <- read_pheno(pheno_file, verbose = FALSE)
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
signal_all <- merge_bwg(rnaseq, riboseq)

## Browser-like view (direct figure return)
browser_figure <- plot_tracks(
  annotation = gp,
  signal = signal_all,
  variants = vcf,
  gene_id = "GeneA",
  signal_type = "bar"
)
browser_figure

## Haplotype extraction
hap <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "string",
  min_variant_number = 1
)

## Haplotype-variant figure (direct figure return)
hap_plot <- plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3
)
hap_plot

## Haplotype-phenotype association
res <- plot_hap_pheno(
  hap,
  phenotype = pheno,
  traits = "seed_weight",
  min_hap_samples = 3
)

res$figure
res$pvalue

## LD block
ld <- compute_ld_block(
  vcf,
  chrom = "chr1",
  start = 12342620,
  end = 12343180,
  verbose = FALSE
)
ld <- plot_ld_block(ld, show_region = TRUE, annotation = gp, show_variant_labels = FALSE)
ld$figure

## Phenotype-guided refinement
refined <- refine_haplotype(
  hap,
  phenotype = pheno,
  traits = "seed_weight",
  min_hap_samples = 3
)
refined$refined_haplotypes

## Variant effects
effect_res <- plot_variant_effect(
  hap,
  phenotype = pheno,
  traits = "protein_content",
  min_group_samples = 3,
  x_axis = "position"
)
effect_res$figure
```

## Notes

- For large indexed VCF files, use `read_vcf(file, mode = "lazy")` and query regions with `retrieve_vcf()` or directly through `hap_gene_variant()` / `hap_region_variant()`.
- For bigWig files, region-based access avoids loading the whole signal file into memory.
- For long signal regions, use `bin_size`, `heatmap_bin_size`, or `heatmap_max_bins` to keep plots readable.
- `hap_variant()` remains available for compatibility, but new code should use `hap_gene_variant()` or `hap_region_variant()`.
- `plot_hap_pheno()` and `plot_variant_pheno()` return structured lists. Use `$figure` for plotting and `$pvalue` for statistical results.
