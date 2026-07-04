# GeneTrackR

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

```r
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

```r
library(GeneTrackR)
```

For local development:

```r
git clone https://github.com/Renscq/GeneTrackR.git
setwd("GeneTrackR")

devtools::document()
devtools::install()
```

## Built-in example files

GeneTrackR ships with example files in `inst/extdata`. The example data are large enough to test the major modules without external files.

```r
gp_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
gtf_file <- system.file("extdata", "example_annotation.gtf", package = "GeneTrackR")
gff_file <- system.file("extdata", "example_annotation.gff3", package = "GeneTrackR")
bed_file <- system.file("extdata", "example_features.bed", package = "GeneTrackR")

bg_files <- system.file(
  "extdata",
  c("example_signal_A.bedgraph", "example_signal_B.bedgraph"),
  package = "GeneTrackR"
)

vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
variant_file <- system.file("extdata", "example_variants.vcf", package = "GeneTrackR")
pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
```

The example set includes approximately 100 genes, multiple chromosomes, coding/non-coding transcripts, positive/negative strands, random signal coverage, around 500 VCF records, 60 samples, and multiple phenotype traits.

## 1. Annotation input

### GenePred / GenePredExt

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = TRUE
)

gp
head(gp$genes)
head(gp$transcripts)
head(gp$exons)
```

Silent loading:

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE,
  progress = FALSE
)
```

### GTF / GFF3 / BED

```r
gtf <- read_gtf(gtf_file)
gff <- read_gff(gff_file)
bed <- read_bed(bed_file)
```

`read_gtf()` and `read_gff()` standardize gene, transcript, exon, CDS, and UTR records when possible. `read_bed()` reads interval-style feature tracks.

### Summarize annotation objects

```r
summary_feature(gp, level = "gene")
summary_feature(gp, level = "transcript")
summary_feature(gp, level = "exon")
```

## 2. Annotation retrieval and conversion

### Retrieve by gene or transcript

```r
gene_feature <- retrieve_feature(
  gp,
  gene_id = "GeneA"
)

tx_feature <- retrieve_feature(
  gp,
  transcript_id = "TxA1"
)
```

### Retrieve by genomic region

```r
region_feature <- retrieve_feature(
  gp,
  chrom = "chr1",
  start = 1,
  end = 20000
)

unique(region_feature$genes$chrom)
```

### Return standardized tables

```r
gene_table <- retrieve_feature(
  gp,
  chrom = "chr1",
  start = 1,
  end = 20000,
  level = "gene",
  as = "data.table"
)

gene_dt <- as_gene_table(gp)
tx_dt <- as_transcript_table(gp)
exon_dt <- as_exon_table(gp)
feature_dt <- as_feature_table(gp)
```

## 3. Gene model plotting

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
  coordinate = "genomic"
)
```

Spliced transcript coordinate mode removes introns:

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
  end = 20000,
  mode = "overlap",
  label_by = "gene"
)
```

### Customize gene model colors

```r
plot_gene(
  gp,
  gene_id = "GeneA",
  color_palette = "Set2",
  border_color = "black"
)

plot_gene(
  gp,
  gene_id = "GeneA",
  fill_colors = c(
    CDS = "#1b9e77",
    UTR = "#a6d854",
    exon = "#7570b3"
  )
)
```

## 4. Signal track module

### Read bedGraph / wig / bigWig

```r
bg <- read_bwg(
  bg_files,
  format = "bedgraph"
)

bg
summary_bwg(bg)
```

For bigWig files, GeneTrackR can use lazy region-based access, so large signal files do not need to be fully loaded before plotting.

### Plot signal over a gene

```r
plot_signal_gene(
  signal = bg,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "bar",
  signal_y_scale = "fixed",
  signal_y_ticks = "pretty"
)
```

### Plot signal over a transcript

```r
plot_signal_transcript(
  signal = bg,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  plot_type = "bar"
)
```

### Heatmap signal with binning

For long regions, `plot_type = "heatmap"` may produce very narrow tiles. Use binning to make the heatmap readable:

```r
plot_signal_gene(
  signal = bg,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "heatmap",
  heatmap_bin_size = 50,
  heatmap_summary = "mean"
)
```

## 5. Variant track module

### Read VCF into memory

```r
vcf <- read_vcf(vcf_file, mode = "memory")
vcf
summary_vcf(vcf)
```

### Lazy indexed VCF access

For large `vcf.gz` files with `.tbi` or `.csi` indexes, use lazy mode. This reads the header and sample names first, then reads variants only for requested regions.

```r
vcf_lazy <- read_vcf(
  "large.vcf.gz",
  mode = "lazy"
)

vcf_region <- retrieve_vcf(
  vcf_lazy,
  chrom = "chr1",
  start = 1,
  end = 20000,
  as = "VariantTrack"
)
```

If no index is available, use memory mode or create a tabix index before lazy retrieval.

### Retrieve and plot variants

```r
vcf_region <- retrieve_vcf(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000
)

plot_variant(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000,
  color_by = "variant_type"
)
```

## 6. Browser-like combined tracks

`plot_tracks()` combines gene models, signal tracks, BED feature tracks, and VCF variant tracks.

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  signal_type = "bar"
)
```

Customize gene model colors in combined tracks:

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  signal_type = "bar",
  gene_color_palette = "Set2",
  gene_border_color = "black"
)

plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  gene_fill_colors = c(
    CDS = "#1b9e77",
    UTR = "#a6d854",
    exon = "#7570b3"
  )
)
```

Add feature and variant tracks:

```r
features <- read_bed(bed_file)
variants <- read_vcf(variant_file)

plot_tracks(
  annotation = gp,
  signal = bg,
  features = features,
  variants = variants,
  chrom = "chr1",
  start = 1,
  end = 20000,
  signal_type = "bar"
)
```

## 7. Haplotype extraction

GeneTrackR now separates gene/transcript-based and region-based haplotype extraction.

### Gene or transcript haplotypes

```r
hap_gene <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "string"
)

hap_tx <- hap_gene_variant(
  vcf,
  annotation = gp,
  transcript_id = "TxA1",
  genotype_mode = "code"
)
```

Include upstream/downstream variants around a gene or transcript:

```r
hap_gene_ext <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  strand_aware = TRUE,
  genotype_mode = "string"
)
```

### Region haplotypes

```r
hap_region <- hap_region_variant(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000,
  genotype_mode = "code"
)
```

`hap_variant()` is retained as a compatibility wrapper, but new code should prefer `hap_gene_variant()` and `hap_region_variant()`.

### Inspect haplotype tables

```r
hap_gene$region
hap_gene$variants
hap_gene$haplotypes
hap_gene$sample_haplotypes
```

## 8. Haplotype-variant plot

`plot_hap_variant()` draws a gene model, variant markers, connector lines, and a genotype table.

```r
plot_hap_variant(
  hap_gene,
  annotation = gp,
  min_hap_samples = 3,
  show_reference_row = TRUE,
  table_x_angle = 90
)
```

Customize table and variant colors:

```r
plot_hap_variant(
  hap_gene,
  annotation = gp,
  min_hap_samples = 3,
  table_fill_palette = "RdBu",
  table_fill_alpha = 0.6,
  variant_palette = "Set1",
  genotype_text_size = 3
)
```

## 9. Phenotype input and summary

The phenotype table should contain sample/taxa IDs in the first column or a named sample column. Each additional column is treated as one trait.

```r
pheno <- read_pheno(pheno_file)

summary_pheno(pheno)

plot_pheno(
  pheno,
  traits = c("plant_height", "seed_weight")
)
```

## 10. Haplotype-phenotype association

`plot_hap_pheno()` compares phenotype distributions among haplotype groups. It returns both the figure and the p-value table.

```r
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

```r
plot_hap_pheno(
  hap = hap_gene,
  phenotype = pheno,
  traits = "seed_weight",
  test_method = "wilcox.test",
  min_hap_samples = 3
)$figure
```

Long trait names are wrapped in facet strips:

```r
plot_hap_pheno(
  hap = hap_gene,
  phenotype = pheno,
  traits = c(
    "plant_height",
    "very_long_trait_name_related_to_seed_weight_under_stress"
  ),
  strip_label_width = 12,
  strip_border_color = NULL,
  min_hap_samples = 3
)$figure
```

## 11. Single-variant phenotype association

`plot_variant_pheno()` is the single-variant version of `plot_hap_pheno()`. It groups samples by genotype or allele state at one variant.

```r
variant_res <- plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  variant_id = "rsA1",
  traits = "plant_height",
  genotype_mode = "string",
  min_group_samples = 3
)

variant_res$figure
variant_res$pvalue
variant_res$variant_data
```

You can also select a variant by genomic position:

```r
plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  chrom = "chr1",
  pos = 1000,
  traits = "plant_height",
  genotype_mode = "code",
  min_group_samples = 3
)$figure
```

## 12. Export module

### Export annotations

```r
write_feature(
  gp,
  file = "example.output.gtf",
  format = "gtf",
  overwrite = TRUE
)

write_feature(
  gp,
  file = "example.output.bed6",
  format = "bed6",
  overwrite = TRUE
)

write_feature(
  gp,
  file = "example.output.bed12",
  format = "bed12",
  overwrite = TRUE
)
```

### Export variants

```r
write_vcf(
  vcf,
  file = "example.output.vcf",
  overwrite = TRUE
)
```

### Save figures and tables

```r
ggplot2::ggsave(
  filename = "hap_pheno.pdf",
  plot = hap_res$figure,
  width = 6,
  height = 5
)

data.table::fwrite(
  hap_res$pvalue,
  file = "hap_pheno.pvalue.tsv",
  sep = "\t"
)
```

## 13. Recommended workflow

```r
library(GeneTrackR)

## Input files
gp_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
bg_files <- system.file(
  "extdata",
  c("example_signal_A.bedgraph", "example_signal_B.bedgraph"),
  package = "GeneTrackR"
)

## Read data
gp <- read_genepred(gp_file, format = "genePredExt", verbose = FALSE)
vcf <- read_vcf(vcf_file, mode = "memory")
pheno <- read_pheno(pheno_file)
bg <- read_bwg(bg_files, format = "bedgraph")

## Browser-like view
plot_tracks(
  annotation = gp,
  signal = bg,
  variants = vcf,
  gene_id = "GeneA",
  signal_type = "bar"
)

## Haplotype extraction
hap <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  genotype_mode = "string"
)

## Haplotype-variant figure
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3
)

## Haplotype-phenotype association
res <- plot_hap_pheno(
  hap,
  phenotype = pheno,
  traits = "plant_height",
  min_hap_samples = 3
)

res$figure
res$pvalue
```

## Notes

- For large indexed VCF files, use `read_vcf(file, mode = "lazy")` and query regions with `retrieve_vcf()` or directly through `hap_gene_variant()` / `hap_region_variant()`.
- For bigWig files, region-based access avoids loading the whole signal file into memory.
- For long signal regions, use `bin_size`, `heatmap_bin_size`, or `heatmap_max_bins` to keep plots readable.
- `hap_variant()` remains available for compatibility, but new code should use `hap_gene_variant()` or `hap_region_variant()`.
- `plot_hap_pheno()` and `plot_variant_pheno()` return structured lists. Use `$figure` for plotting and `$pvalue` for statistical results.
