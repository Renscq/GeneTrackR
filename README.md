# GeneTrackR

GeneTrackR is a lightweight R package for reading, slicing, writing, and visualizing genome annotations, genomic signal tracks, variants, haplotypes, and phenotype associations. It is designed for local genome-browser-like inspection of gene structures together with bedGraph / wig / bigWig coverage, BED/GFF/GTF feature tracks, VCF variant tracks, and haplotype-phenotype results.

The package is useful when you need an IGV-like figure in R, but also want programmable access to the underlying annotation, coverage, variant, haplotype, and phenotype tables.

## Main features

- Read and standardize gene annotations from GenePred, GenePredExt, GTF, GFF3, and BED.
- Read and query bedGraph, wig, and bigWig signal tracks.
- Read, retrieve, merge, plot, and write VCF variant tracks.
- Draw gene models by gene, transcript, or genomic region.
- Draw signal coverage tracks as bar, line, area, or heatmap tracks.
- Combine gene models, signal tracks, feature tracks, and variant tracks using `plot_tracks()`.
- Extract haplotypes from a gene, transcript, or genomic region.
- Draw gene-level haplotype-variant diagrams with variant markers and genotype tables.
- Read phenotype tables and summarize missing values / trait types.
- Draw haplotype-phenotype and single-variant phenotype association plots with pairwise statistical tests.
- Export annotation, signal, and variant tracks to standard formats.

## Installation

```r
# From a local source directory
devtools::document()
devtools::install()

# Or install from a local package archive
remotes::install_local("GeneTrackR.tar.gz", force = TRUE)
```

Load the package:

```r
library(GeneTrackR)
```

## Built-in example files

GeneTrackR ships with example files under `inst/extdata`. These files are large enough to test the main modules without needing external data.

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

The example data include roughly 100 genes, multiple chromosomes, positive/negative strands, coding/non-coding transcripts, random signal coverage across genes, 500 VCF records, 60 samples, and multiple numeric/categorical phenotypes.

## 1. Annotation input module

### Read GenePred / GenePredExt

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = TRUE
)

gp
```

For silent loading:

```r
gp <- read_genepred(
  gp_file,
  format = "genePredExt",
  verbose = FALSE,
  progress = FALSE
)
```

### Read GTF / GFF3 / BED

```r
gtf <- read_gtf(gtf_file)
gff <- read_gff(gff_file)
bed <- read_bed(bed_file)
```

GenePred-like inputs return gene-model-compatible objects with `genes`, `transcripts`, and `exons` tables. Generic feature files return Feature-compatible objects that can be plotted as feature tracks or converted to standardized tables.

### Inspect annotation tables

```r
head(gp$genes)
head(gp$transcripts)
head(gp$exons)

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

### Return a data.table instead of a Feature object

```r
gene_table <- retrieve_feature(
  gp,
  chrom = "chr1",
  start = 1,
  end = 20000,
  level = "gene",
  as = "data.table"
)
```

### Convert to standard tables

```r
gene_dt <- as_gene_table(gp)
tx_dt <- as_transcript_table(gp)
exon_dt <- as_exon_table(gp)
feature_dt <- as_feature_table(gp)
```

## 3. Gene model plotting module

### Plot one gene

```r
plot_gene(
  gp,
  gene_id = "GeneA",
  label_position = "axis",
  direction_mode = "end"
)
```

### Plot one transcript

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

### Highlight a region

```r
plot_region(
  gp,
  chrom = "chr1",
  start = 1,
  end = 20000,
  highlight = data.frame(start = 5000, end = 8000)
)
```

## 4. Signal track module

### Read bedGraph / wig / bigWig signal tracks

```r
bg <- read_bwg(
  bg_files,
  format = "bedgraph",
  mode = "lazy"
)

bg
```

### Retrieve signal in a region

```r
sig_region <- retrieve_bwg(
  bg,
  chrom = "chr1",
  start = 1,
  end = 20000
)

head(sig_region)
```

### Plot signal over a gene

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

### Plot signal over a transcript

```r
plot_signal_transcript(
  bg,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "genomic",
  plot_type = "bar",
  show_gene_model = TRUE
)
```

### Plot signal in transcript coordinate

```r
plot_signal_transcript(
  bg,
  annotation = gp,
  transcript_id = "TxA1",
  coordinate = "transcript",
  plot_type = "line"
)
```

### Plot signal over a genomic region

```r
plot_signal_region(
  bg,
  annotation = gp,
  chrom = "chr1",
  start = 1,
  end = 20000,
  plot_type = "area",
  signal_transform = "sqrt"
)
```

### Heatmap signal with binning

For long regions, heatmap tiles may become too narrow. Use automatic or explicit binning:

```r
plot_signal_gene(
  bg,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "heatmap",
  heatmap_max_bins = 500,
  heatmap_summary = "mean"
)

plot_signal_gene(
  bg,
  annotation = gp,
  gene_id = "GeneA",
  plot_type = "heatmap",
  heatmap_bin_size = 50,
  heatmap_summary = "max"
)
```

### Group and summarize signal samples

```r
groups <- c(
  example_signal_A = "sample_A",
  example_signal_B = "sample_B"
)

plot_signal_gene(
  bg,
  annotation = gp,
  gene_id = "GeneA",
  sample_groups = groups,
  signal_color_by = "group",
  signal_summary = "mean",
  bin_size = 50,
  signal_palette = "Set1"
)
```

## 5. Variant track module

### Read VCF

```r
vcf <- read_vcf(vcf_file)
vcf
```

### Retrieve variants by region

```r
vcf_region <- retrieve_vcf(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000
)

head(vcf_region$data)
```

For indexed large VCF files, `retrieve_vcf()` can query a bgzip-compressed VCF with a `.tbi` index:

```r
vcf_region <- retrieve_vcf(
  "large.vcf.gz",
  chrom = "chr1",
  start = 1,
  end = 20000,
  as = "VariantTrack"
)
```

### Plot variants

```r
plot_variant(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000,
  color_by = "variant_type"
)
```

Show variant IDs:

```r
plot_variant(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000,
  label_by = "variant_id"
)
```

## 6. Combined browser-like plotting with `plot_tracks()`

`plot_tracks()` is the high-level function for combining gene model, signal, feature, and variant tracks.

### Gene + signal by gene ID

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  signal_type = "bar"
)
```

### Gene + signal by transcript ID

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  transcript_id = "TxA1",
  signal_type = "line"
)
```

### Gene + signal by genomic region

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  chrom = "chr1",
  start = 1,
  end = 20000,
  signal_type = "area"
)
```

### Add feature and variant tracks

```r
features <- read_bed(bed_file)
vars_basic <- read_vcf(variant_file)

plot_tracks(
  annotation = gp,
  signal = bg,
  features = features,
  variants = vars_basic,
  chrom = "chr1",
  start = 1,
  end = 20000,
  signal_type = "bar"
)
```

### Customize gene track colors in combined plots

```r
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
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

## 7. Haplotype extraction module

GeneTrackR supports two explicit haplotype extraction interfaces:

- `hap_gene_variant()` for a gene or transcript.
- `hap_region_variant()` for a genomic interval.

`hap_variant()` is retained as a compatibility wrapper.

### Extract haplotypes for a gene

```r
hap <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "string"
)

hap
hap$variants
hap$haplotypes
hap$sample_haplotypes
```

### Extract haplotypes for a transcript

```r
hap_tx <- hap_gene_variant(
  vcf,
  annotation = gp,
  transcript_id = "TxA1",
  genotype_mode = "code"
)
```

### Include upstream and downstream variants

By default, upstream/downstream extension is strand-aware.

```r
hap_updown <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  strand_aware = TRUE,
  genotype_mode = "string"
)
```

For coordinate-based extension independent of strand:

```r
hap_updown2 <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  strand_aware = FALSE,
  genotype_mode = "string"
)
```

### Extract haplotypes from a region

```r
hap_region <- hap_region_variant(
  vcf,
  chrom = "chr1",
  start = 1,
  end = 20000,
  genotype_mode = "code"
)
```

### Genotype display modes

```r
hap_code <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "code"
)

hap_string <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  genotype_mode = "string"
)
```

- `code`: shows `0`, `1`, or `NA`.
- `string`: shows compact allele strings such as `A`, `G`, `i2`, `i6`, or `NA`.

Long insertion/deletion alleles are compressed as `iN`, where `N` is the allele length. For example, `AT/A` becomes `i2/A`, and `G/GTTACA` becomes `G/i6`.

## 8. Haplotype-variant plotting module

`plot_hap_variant()` draws a gene model, variant markers, connector lines, and a haplotype genotype table.

```r
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3
)
```

### Show reference row and compact InDel labels

```r
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3,
  show_reference_row = TRUE,
  genotype_text_size = 3,
  table_x_angle = 90
)
```

### Customize table colors

```r
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3,
  table_fill_palette = "RdBu",
  table_fill_alpha = 0.6,
  reference_fill = "white"
)
```

### Customize variant marker colors

```r
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3,
  variant_palette = "Set2"
)

plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3,
  variant_colors = c(
    SNP = "#1b9e77",
    Insertion = "#d95f02",
    Deletion = "#7570b3"
  )
)
```

### Show genomic coordinate axis above the gene track

```r
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3,
  show_gene_position_axis = TRUE,
  gene_position_axis_n = 5
)
```

The coordinate title includes chromosome information, for example `Chromosome chr1 position (bp)`.

## 9. Phenotype module

### Read phenotype table

The phenotype table should contain a sample column. By default, GeneTrackR expects `sample_id`.

```r
pheno <- read_pheno(pheno_file)
head(pheno)
```

### Summarize phenotype types and missing values

```r
summary_pheno(pheno)
```

### Plot phenotype distributions

```r
plot_pheno(pheno)

plot_pheno(
  pheno,
  traits = c("plant_height", "seed_weight", "protein_content")
)
```

## 10. Haplotype-phenotype association plotting

`plot_hap_pheno()` uses a HapVariant object and a phenotype table. It returns a structured list:

```r
list(
  figure = p,
  pvalue = test_table,
  summary = summary_table,
  bracket = bracket_table,
  plot_data = plot_data
)
```

### Basic haplotype phenotype plot

```r
res <- plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = "plant_height",
  min_hap_samples = 3
)

res$figure
res$pvalue
```

### Multiple traits

```r
res_multi <- plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = c("plant_height", "seed_weight", "protein_content"),
  min_hap_samples = 3,
  strip_label_width = 18
)

res_multi$figure
res_multi$pvalue
```

### Statistical test options

```r
plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = "plant_height",
  test_method = "t.test",
  p_label = "stars",
  p_value_type = "raw"
)$figure

plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = "plant_height",
  test_method = "wilcox.test",
  p_label = "number"
)$figure
```

Supported tests:

- `t.test`
- `wilcox.test`
- `ks.test`

Supported p-value labels:

- `stars`: `*`, `**`, `***`
- `number`: numeric p-value, using scientific notation for very small values
- `both`: stars and numeric p-value

### Customize phenotype plot style

```r
res <- plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = "plant_height",
  plot_type = "violin_boxplot",
  fill_palette = "RdBu",
  show_points = FALSE,
  show_outliers = FALSE,
  x_text_angle = 90
)

res$figure
```

The default plot is a violin plot with a narrow boxplot in the middle. Haplotype groups are ordered by sample number from left to right, and x-axis labels are shown as `Hap (n)`.

## 11. Single-variant phenotype association plotting

`plot_variant_pheno()` is the single-variant version of `plot_hap_pheno()`. It groups samples by genotype/allele state at one variant site.

### Plot by variant ID

```r
res_var <- plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  variant_id = "rsA1",
  traits = "plant_height",
  min_group_samples = 3
)

res_var$figure
res_var$pvalue
res_var$variant_data
```

### Plot by chromosome position

```r
res_var_pos <- plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  chrom = "chr1",
  pos = 5000,
  traits = "plant_height",
  genotype_mode = "string",
  min_group_samples = 3
)

res_var_pos$figure
res_var_pos$pvalue
```

### Use a VCF file path directly

```r
res_var_file <- plot_variant_pheno(
  variant = vcf_file,
  phenotype = pheno,
  variant_id = "rsA1",
  traits = "plant_height",
  min_group_samples = 3
)
```

For indexed large VCF files, pass `chrom` and `pos` to retrieve only the target site.

## 12. Writing output files

### Write annotation files

```r
write_feature(
  gp,
  file = "example_output.gtf",
  format = "gtf",
  overwrite = TRUE
)

write_feature(
  gp,
  file = "example_output.bed12",
  format = "bed12",
  overwrite = TRUE
)

write_feature(
  gp,
  file = "example_output.bed6",
  format = "bed6",
  overwrite = TRUE
)
```

### Write VCF

```r
write_vcf(
  vcf,
  file = "example_output.vcf",
  overwrite = TRUE
)
```

### Save figures

```r
p <- plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA"
)

ggplot2::ggsave(
  filename = "GeneA_tracks.pdf",
  plot = p,
  width = 8,
  height = 5
)

res <- plot_hap_pheno(
  hap = hap,
  phenotype = pheno,
  traits = "plant_height",
  min_hap_samples = 3
)

ggplot2::ggsave(
  filename = "GeneA_hap_pheno.pdf",
  plot = res$figure,
  width = 6,
  height = 5
)

data.table::fwrite(
  res$pvalue,
  file = "GeneA_hap_pheno_pvalue.tsv",
  sep = "\t"
)
```

## 13. Suggested analysis workflow

A typical workflow is:

```r
library(GeneTrackR)

# 1. Read annotation, signal, VCF, and phenotype files
gp <- read_genepred(gp_file, format = "genePredExt", verbose = FALSE)
bg <- read_bwg(bg_files, format = "bedgraph", mode = "lazy")
vcf <- read_vcf(vcf_file)
pheno <- read_pheno(pheno_file)

# 2. Inspect a locus with gene model and signal coverage
plot_tracks(
  annotation = gp,
  signal = bg,
  gene_id = "GeneA",
  signal_type = "bar"
)

# 3. Extract gene-level haplotypes
hap <- hap_gene_variant(
  vcf,
  annotation = gp,
  gene_id = "GeneA",
  upstream = 1000,
  downstream = 1000,
  genotype_mode = "string"
)

# 4. Visualize variants and haplotype table
plot_hap_variant(
  hap,
  annotation = gp,
  min_hap_samples = 3
)

# 5. Test haplotype-phenotype association
hap_pheno <- plot_hap_pheno(
  hap,
  phenotype = pheno,
  traits = "plant_height",
  min_hap_samples = 3
)

hap_pheno$figure
hap_pheno$pvalue

# 6. Test a single variant-phenotype association
variant_pheno <- plot_variant_pheno(
  variant = vcf,
  phenotype = pheno,
  variant_id = "rsA1",
  traits = "plant_height",
  min_group_samples = 3
)

variant_pheno$figure
variant_pheno$pvalue
```

## Notes

- GenePred / GTF / GFF3 annotations are represented internally as standardized Feature/GenePred-compatible objects.
- Coordinates in gene model objects are 1-based closed unless otherwise specified.
- BED output follows standard BED conventions: 0-based half-open intervals.
- `plot_tracks()` requires exactly one locator: `gene_id`, `transcript_id`, or `chrom + start + end`.
- `plot_hap_pheno()` and `plot_variant_pheno()` return both the figure and the statistical result table.
- For large VCF files, bgzip compression and tabix indexing are recommended.
