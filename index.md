# Any annotation object -\> standardized feature table

**GeneTrackR** is a lightweight R package for reading, querying,
writing, and visualizing gene annotations, genomic signal tracks,
variant tracks, haplotypes, and phenotype associations. It is designed
for programmable IGV-like visualization in R while keeping direct access
to the underlying annotation, coverage, variant, haplotype, and
phenotype tables.

## Main features

- Read and standardize annotations from **GenePred**, **GenePredExt**,
  **GTF**, **GFF3**, and **BED**.
- Read and query **bedGraph**, **wig**, and **bigWig** signal tracks.
- Read, lazily query, retrieve, merge, plot, and write **VCF** variant
  tracks.
- Draw gene models by gene, transcript, or genomic region.
- Draw signal tracks as bar, line, area, or heatmap tracks.
- Combine gene models, signal tracks, feature tracks, and variant tracks
  with
  [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md).
- Extract haplotypes from genes, transcripts, or arbitrary genomic
  regions.
- Draw haplotype-variant diagrams with gene models, variant markers,
  connector lines, and genotype tables.
- Read phenotype tables, summarize missingness/type, and draw phenotype
  distributions.
- Draw haplotype-phenotype and single-variant phenotype association
  plots with pairwise statistical tests.
- Export annotations, signal tracks, and variants to standard formats.

## Installation

GeneTrackR depends on CRAN and Bioconductor packages. Install the
dependencies first, then install GeneTrackR from GitHub.

\`\`\`{r} \## CRAN dependencies install.packages(c( “devtools”,
“data.table”, “ggplot2”, “patchwork”, “rlang”, “Rcpp”, “RColorBrewer” ))

## Bioconductor dependencies

if (!requireNamespace(“BiocManager”, quietly = TRUE)) {
install.packages(“BiocManager”) }

BiocManager::install(c( “GenomicRanges”, “IRanges”, “Rsamtools” ), ask =
FALSE, update = FALSE)

## Install GeneTrackR from GitHub

devtools::install_github( “Renscq/GeneTrackR”, dependencies = TRUE,
build_vignettes = FALSE )


    Load the package:

    ```{r}
    library(GeneTrackR)

For local development:

``` bash
git clone https://github.com/Renscq/GeneTrackR.git
cd GeneTrackR
```

`{r} devtools::document() devtools::install()`

## Built-in example files

GeneTrackR ships with one deterministic demo genome in `inst/extdata`.
The dataset is compact but intentionally designed to exercise the major
modules without external files.

\`\`\`{r} gp_file \<- system.file(“extdata”, “gtr_demo.genePredExt”,
package = “GeneTrackR”) gtf_file \<- system.file(“extdata”,
“gtr_demo.gtf”, package = “GeneTrackR”) gff_file \<-
system.file(“extdata”, “gtr_demo.gff3”, package = “GeneTrackR”) bed_file
\<- system.file(“extdata”, “gtr_demo_features.bed”, package =
“GeneTrackR”)

rnaseq_files \<- system.file( “extdata”,
c(“gtr_demo_rnaseq_plus.bedgraph”, “gtr_demo_rnaseq_minus.bedgraph”),
package = “GeneTrackR” )

riboseq_files \<- system.file( “extdata”,
c(“gtr_demo_riboseq_plus.bedgraph”, “gtr_demo_riboseq_minus.bedgraph”),
package = “GeneTrackR” )

vcf_file \<- system.file(“extdata”, “gtr_demo_variants.vcf”, package =
“GeneTrackR”) pheno_file \<- system.file(“extdata”,
“gtr_demo_pheno.tsv”, package = “GeneTrackR”)


    The deterministic demo genome contains 2 chromosomes, 20 genes, 24 transcripts, 36 samples, 56 designed variants, strand-specific RNA-seq/Ribo-seq tracks, four balanced GeneA haplotype groups, and phenotype traits with known positive and negative-control associations. The GeneA genotype design is hierarchical: `Hap1/Hap2` are the closest pair and `Hap3/Hap4` are the second closest pair. `protein_content` follows the same two-cluster structure through `varA03`, while `seed_weight` and `plant_height` retain four distinct but hierarchically ordered haplotype means. All protein-coding demo transcripts use CDS lengths divisible by three. The Ribo-seq tracks contain moderately dense heterogeneous integer P-site-like counts with designed frame0/frame1/frame2 total-count proportions of approximately 80%/10%/10%. Frame 0 is broadly occupied and variable, while frame 1 and frame 2 use different subsets of codons with lower irregular counts. Initiation and termination frame-0 counts are approximately two times the internal frame-0 mean.

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

### GenePred and GenePredExt

GenePred stores one transcript per row with transcript/CDS boundaries
and exon blocks. Standard GenePred contains 10 core columns. GenePredExt
adds `score`, `name2`, CDS status, and exon-frame information;
GeneTrackR normally uses `name` as `transcript_id` and `name2` as
`gene_id` for GenePredExt.

#### Read

The built-in file is GenePredExt:

\`\`\`{r} gp \<- read_genepred( gp_file, format = “genePredExt”, verbose
= FALSE, progress = FALSE )

gp head(gp$`genes)
head(gp`$transcripts) head(gp\$exons)


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

#### Subset

GenePred-compatible objects can be subset by gene, transcript, or
genomic region:

\`\`\`{r} gp_gene_a \<- retrieve_feature( gp, gene_id = “GeneA” )

gp_tx_a1 \<- retrieve_feature( gp, transcript_id = “TxA1” )

gp_region \<- retrieve_feature( gp, chrom = “chr1”, start = 12339001,
end = 12374500, mode = “overlap” )


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

#### Merge

Subsets from the same or different GenePred objects can be recombined
with
[`merge_feature()`](https://renscq.github.io/GeneTrackR/reference/merge_feature.md):

\`\`\`{r} gp_gene_b \<- retrieve_feature(gp, gene_id = “GeneB”)

gp_merged \<- merge_feature( gp_gene_a, gp_gene_b, source_names =
c(“GeneA”, “GeneB”), conflict = “deduplicate” )


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

[`write_genepred()`](https://renscq.github.io/GeneTrackR/reference/write_genepred.md)
remains available as a backward-compatible wrapper around
[`write_feature()`](https://renscq.github.io/GeneTrackR/reference/write_feature.md).

### GTF

GTF stores gene-model records in a 9-column table and uses the final
attribute column for identifiers such as `gene_id` and `transcript_id`.
GeneTrackR parses the hierarchical records into a `FeatureTrack` with
standardized `$genes`, `$transcripts`, `$exons`, and `$data` tables when
the required records are present.

#### Read

\`\`\`{r} gtf \<- read_gtf( gtf_file, verbose = FALSE, progress = FALSE
)

gtf head(gtf$`genes)
head(gtf`$transcripts)


    Selected feature types can be loaded when the full annotation is not required:

    ```{r}
    gtf_core <- read_gtf(
      gtf_file,
      feature_types = c("gene", "transcript", "exon", "CDS"),
      verbose = FALSE,
      progress = FALSE
    )

#### Subset

\`\`\`{r} gtf_gene_a \<- retrieve_feature(gtf, gene_id = “GeneA”)

gtf_gene_a_exons \<- retrieve_feature( gtf, gene_id = “GeneA”, level =
“exon”, as = “data.table” )

gtf_region \<- retrieve_feature( gtf, chrom = “chr1”, start = 12339001,
end = 12374500, mode = “overlap” )


    #### Merge

    ```{r}
    gtf_gene_b <- retrieve_feature(gtf, gene_id = "GeneB")

    gtf_merged <- merge_feature(
      gtf_gene_a,
      gtf_gene_b,
      source_names = c("GeneA", "GeneB")
    )

#### Write

`{r} write_feature( gtf_gene_a, file.path(tempdir(), "GeneA.gtf"), format = "gtf", overwrite = TRUE )`

Because a gene-model GTF contains transcript and exon hierarchy, it can
also be converted directly to GenePred/GenePredExt or BED12 with
[`write_feature()`](https://renscq.github.io/GeneTrackR/reference/write_feature.md).

### GFF3

GFF3 also uses 9 columns, but hierarchy is represented primarily through
`ID` and `Parent` attributes.
[`read_gff()`](https://renscq.github.io/GeneTrackR/reference/read_gff.md)
reconstructs gene/transcript/exon relationships and returns the same
unified `FeatureTrack` interface used for GTF.

#### Read

\`\`\`{r} gff \<- read_gff( gff_file, verbose = FALSE, progress = FALSE
)

gff head(gff$`genes)
head(gff`$transcripts)


    Feature-type filtering is also supported:

    ```{r}
    gff_core <- read_gff(
      gff_file,
      feature_types = c("gene", "mRNA", "exon", "CDS"),
      verbose = FALSE,
      progress = FALSE
    )

#### Subset

\`\`\`{r} gff_gene_a \<- retrieve_feature(gff, gene_id = “GeneA”)

gff_region \<- retrieve_feature( gff, chrom = “chr1”, start = 12339001,
end = 12374500, mode = “overlap” )


    #### Merge

    ```{r}
    gff_gene_b <- retrieve_feature(gff, gene_id = "GeneB")

    gff_merged <- merge_feature(
      gff_gene_a,
      gff_gene_b,
      source_names = c("GeneA", "GeneB")
    )

#### Write

`{r} write_feature( gff_gene_a, file.path(tempdir(), "GeneA.gff3"), format = "gff", overwrite = TRUE )`

### BED

BED is treated differently from GenePred/GTF/GFF3.
[`read_bed()`](https://renscq.github.io/GeneTrackR/reference/read_bed.md)
reads BED3-BED12-style files as **interval annotations** and
standardizes the interval coordinates, name, score, and strand fields
into a `FeatureTrack`. It does not reconstruct BED12 block columns into
a transcript/exon hierarchy; use GenePred, GTF, or GFF3 when exon-level
gene models are required.

#### Read

\`\`\`{r} features \<- read_bed( bed_file, verbose = FALSE, progress =
FALSE )

features head(features\$data)


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

#### Merge

\`\`\`{r} bed_gene_b \<- retrieve_feature( features, chrom = “chr1”,
start = 12355501, end = 12368500, mode = “overlap” )

bed_merged \<- merge_feature( bed_gene_a, bed_gene_b, source_names =
c(“GeneA_region”, “GeneB_region”), conflict = “keep_all” )


    #### Write

    Generic BED-derived intervals should normally be written as BED6:

    ```{r}
    write_feature(
      bed_merged,
      file.path(tempdir(), "merged_features.bed6"),
      format = "bed6",
      overwrite = TRUE
    )

BED12 output requires transcript/exon gene-model information. It is
therefore appropriate for `GenePred` or gene-model GTF/GFF3 objects
rather than a generic BED interval track:

`{r} write_feature( gp_gene_a, file.path(tempdir(), "GeneA.bed12"), format = "bed12", overwrite = TRUE )`

### Region-selection modes

[`retrieve_feature()`](https://renscq.github.io/GeneTrackR/reference/retrieve_feature.md)
uses the same region-selection interface across annotation formats:

| `mode` | Meaning |
|----|----|
| `"overlap"` | Keep records overlapping the requested region. This is the default. |
| `"within"` | Keep records fully contained in the requested region. |
| `"trim"` | For hierarchical gene models, keep overlapping records and clip transcript/exon boundaries to the requested interval. |

For generic BED-style interval tracks, `"trim"` behaves as overlap
selection because there is no transcript/exon hierarchy to rebuild.

### Merge conflict handling

[`merge_feature()`](https://renscq.github.io/GeneTrackR/reference/merge_feature.md)
accepts `Feature`, `FeatureTrack`, and `GenePred` objects, including a
list of objects. Duplicate identifiers are handled explicitly:

| `conflict` | Behavior |
|----|----|
| `"deduplicate"` | Keep the first compatible copy according to input order. Default. |
| `"rename"` | Rename conflicting IDs in later inputs and update hierarchy references. |
| `"keep_all"` | Retain conflicting records unchanged. |
| `"error"` | Stop when duplicate IDs are detected. |
| `"keep_first"` | Backward-compatible alias of `"deduplicate"`. |

When the same locus is loaded from different annotation formats,
`rename` is useful if both representations should be retained:

`{r} cross_format <- merge_feature( list(gp_gene_a, gtf_gene_a), source_names = c("GenePred", "GTF"), conflict = "rename" )`

If the two inputs represent the same annotation and only one copy is
needed, use `conflict = "deduplicate"` instead.

### Cross-format conversion and standardized tables

The readers share a common internal representation, so annotation
formats can be converted without reparsing text manually:

\`\`\`{r} \# GTF/GFF3 gene models -\> GenePred-compatible object
gtf_as_gp \<- as_genepred(gtf) gff_as_gp \<- as_genepred(gff)

feature_dt \<- as_feature_table(gp)

# Gene-model tables

gene_dt \<- as_gene_table(gp) tx_dt \<- as_transcript_table(gp) exon_dt
\<- as_exon_table(gp)

# Bioconductor interoperability

gene_gr \<- as_granges(gp, level = “gene”)


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

### Annotation object contract

After reading, GeneTrackR keeps annotation data in a small set of
compatible objects:

| Input | Main object | Important contents |
|----|----|----|
| GenePred / GenePredExt | `GenePred` | `$genes`, `$transcripts`, `$exons`, `$data` |
| Gene-model GTF / GFF3 | `FeatureTrack` | `$data` plus derived `$genes`, `$transcripts`, `$exons` |
| BED | interval `FeatureTrack` | primarily `$data` |

[`retrieve_feature()`](https://renscq.github.io/GeneTrackR/reference/retrieve_feature.md)
returns a compatible sub-object by default (`as = "Feature"`) and a
`data.table` when `as = "data.table"`.
[`merge_feature()`](https://renscq.github.io/GeneTrackR/reference/merge_feature.md)
always returns a unified `Feature` object and additionally inherits from
`GenePred` when transcript/exon tables are available.

### How annotation objects connect to the rest of GeneTrackR

Annotation objects are the first input to most downstream workflows.
`GenePred` or gene-model `FeatureTrack` objects can be passed directly
to gene-model plotting, signal plotting, haplotype extraction, LD
visualization, and integrated browser-style tracks.

### Core object flow and return contracts

GeneTrackR uses a small set of S3 objects throughout the workflow. The
examples below keep these objects intact instead of repeatedly
converting them to plain tables.

| Step | Function | Main return | Important contents / next consumer |
|----|----|----|----|
| Annotation | [`read_genepred()`](https://renscq.github.io/GeneTrackR/reference/read_genepred.md) | `GenePred` | `$genes`, `$transcripts`, `$exons`, `$data`; used by gene/signal/haplotype plotting |
| GTF/GFF/BED | [`read_gtf()`](https://renscq.github.io/GeneTrackR/reference/read_gtf.md), [`read_gff()`](https://renscq.github.io/GeneTrackR/reference/read_gff.md), [`read_bed()`](https://renscq.github.io/GeneTrackR/reference/read_bed.md) | `FeatureTrack` | `$data` plus derived hierarchy when available; gene-model FeatureTracks are GenePred-convertible |
| Signal | [`read_bwg()`](https://renscq.github.io/GeneTrackR/reference/read_bwg.md) / [`merge_bwg()`](https://renscq.github.io/GeneTrackR/reference/merge_bwg.md) | `BwgTrack` | `$samples`, optional in-memory `$data`, `$meta`; consumed by `plot_signal_*()` and [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md) |
| Variant | [`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md) | `VariantTrack` | `$data`, `$meta`; consumed by haplotype, LD, phenotype, and track functions |
| Variant subset | [`retrieve_vcf()`](https://renscq.github.io/GeneTrackR/reference/retrieve_vcf.md) | `data.table` by default; `VariantTrack` with `as = "VariantTrack"` | use the table for inspection or request `VariantTrack` for downstream track/object workflows |
| Haplotype | [`hap_gene_variant()`](https://renscq.github.io/GeneTrackR/reference/hap_gene_variant.md) / [`hap_region_variant()`](https://renscq.github.io/GeneTrackR/reference/hap_region_variant.md) | `HapVariant` | `$variants`, `$haplotypes`, `$sample_haplotypes`, `$genotype_wide` |
| Phenotype association | [`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md) / [`plot_variant_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_pheno.md) | `GeneTrackRPhenoPlot` | `$figure`, `$pvalue`, `$summary`, `$plot_data` |
| LD | [`compute_ld_block()`](https://renscq.github.io/GeneTrackR/reference/compute_ld_block.md) | `LDTrack` | `$data`, `$matrix`, `$variants`, `$region`; [`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md) stores the figure in `$figure` |
| Refinement | [`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md) | `HapRefined` | `$refined_hap`, `$refined_haplotypes`, `$haplotype_map`, `$pairwise_test` |
| Variant effect | [`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md) | `GeneTrackRVariantEffectPlot` | `$figure`, `$effect`, `$plot_data` |

Plotting functions such as
[`plot_gene()`](https://renscq.github.io/GeneTrackR/reference/plot_gene.md),
[`plot_signal_gene()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_gene.md),
[`plot_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_variant.md),
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md),
[`plot_hap_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_variant.md),
and
[`plot_refined_hap_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_refined_hap_variant.md)
return a ggplot/patchwork figure directly. In contrast, phenotype/effect
plotting functions return result objects containing `$figure` plus
analysis tables.
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
is intentionally stateful: by default it returns an updated `LDTrack`
with the plot stored in `$figure`; use `return_object = FALSE` only when
a figure-only return is required.

## Gene model plotting

### Plot a gene

`{r} plot_gene( gp, gene_id = "GeneA", label_position = "axis", direction_mode = "end" )`

### Plot a transcript

`{r} plot_transcript( gp, transcript_id = "TxA1", coordinate = "genomic" )`

Spliced transcript coordinate mode removes introns:

`{r} plot_transcript( gp, transcript_id = "TxA1", coordinate = "transcript" )`

### Plot a genomic region

`{r} plot_region( gp, chrom = "chr1", start = 12339001, end = 12374500, mode = "overlap", label_by = "gene" )`

### Customize gene model colors

\`\`\`{r} plot_gene( gp, gene_id = “GeneA”, gene_palette = “Set2”,
gene_border_color = “black” )

plot_gene( gp, gene_id = “GeneA”, gene_colors = c( CDS = “#1b9e77”, UTR
= “#a6d854”, exon = “#7570b3” ) )


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

Read the RNA-seq plus/minus bedGraph files into one `BwgTrack` object:

\`\`\`{r} rnaseq \<- read_bwg( rnaseq_files, format = “bedgraph”,
sample_names = c(“RNA_seq_plus”, “RNA_seq_minus”), strand = c(“+”, “-”),
mode = “memory” )

rnaseq summary_bwg(rnaseq)


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

The two assays intentionally have different signal structures:

- **RNA-seq** signal is enriched across exons, including UTRs, while
  intronic and intergenic regions have little or no coverage in the demo
  data.
- **Ribo-seq** signal uses moderately dense 1-bp integer P-site-like
  counts within CDS regions. Frame 0 is broadly occupied with variable
  heights; frame 1 and frame 2 occur at different subsets of codons with
  lower irregular counts. Zero-count bases are omitted from bedGraph.
  Total counts remain approximately 80%/10%/10% for frame 0/frame
  1/frame 2, and the initiation/termination frame-0 counts are
  approximately two times the internal frame-0 mean.

Because bedGraph does not store strand metadata, `strand = c("+", "-")`
explicitly records the strand associated with each input file in the
`BwgTrack` sample table.

### Step 2. Write RNA-seq and Ribo-seq bigWig files

[`write_bwg()`](https://renscq.github.io/GeneTrackR/reference/write_bwg.md)
can convert an in-memory `BwgTrack` to bigWig directly with the bundled
third-party libBigWig library in `src/`. No external conversion program
is required; only chromosome sizes are needed for the bigWig header.

\`\`\`{r} chrom_sizes_file \<- system.file( “extdata”,
“gtr_demo.chrom.sizes”, package = “GeneTrackR” )

bigwig_dir \<- file.path(tempdir(), “GeneTrackR_demo_bigwig”)
dir.create(bigwig_dir, recursive = TRUE, showWarnings = FALSE)


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

Write the Ribo-seq tracks:

\`\`\`{r} riboseq_bigwig \<- write_bwg( riboseq, outdir = bigwig_dir,
format = “bigwig”, chrom_sizes = chrom_sizes_file, overwrite = TRUE )

riboseq_bigwig


    Each call invisibly returns a table with `sample_id`, output `file`, and `format`. With the sample names used above, the output directory contains:

    ```text
    RNA_seq_plus.bigwig
    RNA_seq_minus.bigwig
    Ribo_seq_plus.bigwig
    Ribo_seq_minus.bigwig

BigWig export has a single backend in GeneTrackR: the bundled libBigWig
implementation. This keeps the write path deterministic across platforms
and avoids an external executable dependency.

### Step 3. Plot RNA-seq and Ribo-seq tracks for a gene

[`plot_signal_gene()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_gene.md)
retrieves the genomic span of a gene and optionally adds the gene model
below the signal panel. With `strand = "auto"`, the gene strand is used
to select the matching signal track. `GeneA` is on the positive strand,
so the following examples use the plus RNA-seq and Ribo-seq samples
automatically.

RNA-seq gene track:

\`\`\`{r} p_rnaseq_gene \<- plot_signal_gene( signal = rnaseq,
annotation = gp, gene_id = “GeneA”, plot_type = “bar”, strand = “auto”,
signal_palette = “Blues”, signal_palette_direction = -1, signal_y_scale
= “fixed”, signal_y_ticks = “pretty”, signal_track_height = 3,
gene_track_height = 1 )

p_rnaseq_gene


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

At a whole-gene genomic scale, thousands of bases are compressed into
the available plot width, so a three-nucleotide pattern cannot be
resolved visually even when the underlying RPF counts are frame-biased.
The gene-level bar track is intended to show where translation signal
occurs across CDS exons. Use the transcript-level `frame` view below, or
a short genomic window, to inspect three-nucleotide periodicity.

The relative vertical space occupied by the signal and gene-model panels
is controlled directly by `signal_track_height` and `gene_track_height`.

### Step 4. Plot RNA-seq and Ribo-seq tracks for a transcript

[`plot_signal_transcript()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_transcript.md)
focuses on a single transcript. The RNA-seq example uses the standard
bar representation in transcript coordinates:

\`\`\`{r} p_rnaseq_transcript \<- plot_signal_transcript( signal =
rnaseq, annotation = gp, transcript_id = “TxA1”, coordinate =
“transcript”, plot_type = “bar”, strand = “auto”, signal_palette =
“Blues”, signal_palette_direction = -1, signal_y_scale = “fixed”,
signal_track_height = 3, gene_track_height = 1 )

p_rnaseq_transcript


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

Use `plot_type = "bar"` instead when the goal is to compare the genomic
distribution of Ribo-seq P-site counts with the RNA-seq coverage
representation:

`{r} plot_signal_transcript( signal = riboseq, annotation = gp, transcript_id = "TxA1", coordinate = "transcript", plot_type = "bar", strand = "auto", signal_palette = "Reds", signal_palette_direction = -1 )`

### Step 5. Plot RNA-seq and Ribo-seq together with `plot_tracks()`

Merge the two `BwgTrack` objects so RNA-seq and Ribo-seq can be
displayed in one integrated track figure:

\`\`\`{r} signal_all \<- merge_bwg(rnaseq, riboseq)

signal_all summary_bwg(signal_all)


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

[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
uses its named `heights` vector to control panel proportions. The
`signal` and `gene` entries play the same role as `signal_track_height`
and `gene_track_height` in the dedicated signal plotting functions.

### Step 6. Plot RNA-seq and Ribo-seq tracks across a genomic region

[`plot_signal_region()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_region.md)
is useful when the region contains several genes or when the target
interval is not defined by a single gene/transcript.

For `bar` and `line` region plots, sample/group colors follow the
sample/group level order and the standard RColorBrewer class order.
Discrete colors are not selected by interpolating between the first and
last palette colors. Heatmaps remain continuous gradients.

The following region contains the positive-strand `GeneA`,
negative-strand `GeneB`, and positive-strand `GeneC`, making it useful
for displaying both strand-specific RNA-seq tracks:

\`\`\`{r} p_rnaseq_region \<- plot_signal_region( signal = rnaseq,
annotation = gp, chrom = “chr1”, start = 12339001, end = 12374500,
strand = “both”, plot_type = “bar”, signal_palette = “Blues”,
signal_palette_direction = -1, signal_y_scale = “free”, signal_y_ticks =
“pretty”, signal_track_height = 4, gene_track_height = 1 )

p_rnaseq_region


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

For large real-world signal files, use lazy `BwgTrack` access where
appropriate so that only the requested genomic interval is retrieved
during plotting.

## VCF variant tracks

GeneTrackR stores VCF records in a `VariantTrack` object. This module
follows one continuous workflow for the demo VCF:

1.  read VCF data and inspect the `VariantTrack` object;
2.  validate and summarize variants;
3.  retrieve variants from a genomic region;
4.  retrieve variants by gene or transcript;
5.  filter variants by ID, type, or text pattern;
6.  plot variant tracks.

The examples use the same `GeneA` / `TxA1` locus used elsewhere in the
documentation so that variant, haplotype, LD, phenotype, and signal
examples can be compared directly.

### Step 1. Read VCF data

Locate the demo VCF and annotation files:

\`\`\`{r} vcf_file \<- system.file( “extdata”, “gtr_demo_variants.vcf”,
package = “GeneTrackR” )

gp_file \<- system.file( “extdata”, “gtr_demo.genePredExt”, package =
“GeneTrackR” )

gp \<- read_genepred( gp_file, format = “genePredExt”, verbose = FALSE,
progress = FALSE )


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

[`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md)
returns a `VariantTrack`. Standard VCF fields are normalized internally
to columns such as `chrom`, `pos`, `variant_id`, `ref`, `alt`, `qual`,
`filter`, `info`, and `variant_type`. Genotype columns are retained when
`keep_genotype = TRUE`.

\`\`\`{r} head(vcf\$data\[, c( “chrom”, “pos”, “variant_id”, “ref”,
“alt”, “variant_type” )\])

vcf$`meta`$sample_names


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

Lazy regional queries require an indexed VCF and `Rsamtools`. For small
files such as the bundled demo VCF, memory mode is simpler.

### Step 2. Validate and summarize variants

Use
[`validate_vcf()`](https://renscq.github.io/GeneTrackR/reference/validate_vcf.md)
to check required fields, coordinates, allele fields, duplicated variant
IDs, and consistency between `pos`, `start`, and `end`:

\`\`\`{r} vcf_validation \<- validate_vcf(vcf)

vcf_validation$`invalid_summary
vcf_validation`$warnings


    The validation result contains:

    - `invalid_records`: individual records that failed a check;
    - `invalid_summary`: number of invalid records by reason;
    - `warnings`: non-record-specific validation messages.

    Use `summary_vcf()` to summarize the current variant collection. A genomic region is **not required** for an in-memory `VariantTrack`; with no range supplied, the complete object is summarized. By default, variants are counted by chromosome and inferred variant type:

    ```{r}
    vcf_summary <- summary_vcf(vcf)
    vcf_summary

A chromosome or complete genomic range can still be supplied when only a
subset should be summarized:

`{r} chr1_summary <- summary_vcf(vcf, chrom = "chr1") chr1_summary`

The demo VCF intentionally contains SNPs and indels so the same file can
be reused for variant visualization, haplotype analysis, and LD
analysis.

### Step 3. Retrieve variants from a genomic region

[`retrieve_vcf()`](https://renscq.github.io/GeneTrackR/reference/retrieve_vcf.md)
can return either a plain `data.table` or a new `VariantTrack`.

Retrieve the broader `GeneA` / `GeneB` / `GeneC` demonstration region as
a table:

\`\`\`{r} region_table \<- retrieve_vcf( vcf, chrom = “chr1”, start =
12339001, end = 12374500, as = “data.table”, verbose = FALSE )

region_table\[, .( chrom, pos, variant_id, ref, alt, variant_type )\]


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

The distinction is important:

- `as = "data.table"` is convenient for inspection, filtering, and
  tabular analysis;
- `as = "VariantTrack"` preserves the GeneTrackR object interface for
  downstream plotting and analysis.

### Step 4. Retrieve variants by gene or transcript

When annotation is available,
[`retrieve_vcf()`](https://renscq.github.io/GeneTrackR/reference/retrieve_vcf.md)
can resolve genomic coordinates directly from a gene or transcript ID.

Retrieve variants inside the `GeneA` gene body:

\`\`\`{r} genea_variants \<- retrieve_vcf( vcf, annotation = gp, gene_id
= “GeneA”, as = “VariantTrack”, verbose = FALSE )

genea_variants


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

For `GeneA`, the upstream extension additionally includes `varAup01`,
which is intentionally designed with missing and heterozygous genotypes
for edge-case testing.

Transcript-aware retrieval uses the same interface:

\`\`\`{r} txa1_variants \<- retrieve_vcf( vcf, annotation = gp,
transcript_id = “TxA1”, as = “VariantTrack”, verbose = FALSE )

txa1_variants


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

Filter by inferred variant type:

\`\`\`{r} genea_indels \<- retrieve_vcf( genea_variants, variant_type =
c(“INS”, “DEL”), as = “data.table”, verbose = FALSE )

genea_indels\[, .( pos, variant_id, ref, alt, variant_type )\]


    `pattern` searches `variant_id`, `REF`, `ALT`, `INFO`, and `variant_type`. Fixed-string matching supports `ignore_case = TRUE` without regular-expression interpretation. For example, the demo high-LD variants contain `high_ld` in the INFO role field:

    ```{r}
    high_ld_variants <- retrieve_vcf(
      vcf,
      pattern = "high_ld",
      fixed = TRUE,
      as = "VariantTrack",
      verbose = FALSE
    )

    high_ld_variants

These filters can be combined with genomic, gene, or transcript queries
when a more specific subset is needed.

### Step 6. Plot variant tracks

[`plot_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_variant.md)
returns a ggplot object directly. The GeneA region contains SNP,
insertion, and deletion examples, so it is useful for demonstrating
variant-type colors.

\`\`\`{r} p_genea_variants \<- plot_variant( genea_variants, color_by =
“variant_type”, label_by = “variant_id”, variant_shape = “lollipop”,
variant_palette = “Set1”, point_size = 2.5, text_size = 12 )

p_genea_variants


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

The main plotting controls are:

- `color_by`: `"variant_type"`, `"filter"`, or `"none"`;
- `label_by`: variant ID, type, REF, ALT, or no labels;
- `variant_shape`: `"lollipop"`, `"point"`, or `"rug"`;
- `variant_palette` / `variant_colors`: automatic palette or explicit
  colors.

Variant subsets created here are reused directly by the later haplotype,
phenotype, LD, and integrated-track modules.

VCF export is covered separately in the **Export variants and analysis
results** module. The current
[`write_vcf()`](https://renscq.github.io/GeneTrackR/reference/write_vcf.md)
interface writes standardized site-level VCF fields; the original
genotype-containing VCF should be retained as the archival source when
genotype round-trip preservation is required.

## Browser-like combined tracks

[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
assembles gene models, genomic signal, interval features, and variants
into one genome-browser-like figure with a shared genomic x-axis. It is
intended for integrated inspection of a locus rather than for replacing
the dedicated plotting functions for each data type.

This module follows one continuous workflow:

1.  prepare annotation, signal, feature, and variant tracks;
2.  build a gene-centered browser view;
3.  build a transcript-centered browser view;
4.  inspect an explicit genomic region;
5.  add BED features and VCF variants;
6.  control layout, panel heights, highlighting, and styling;
7.  choose an appropriate browser workflow for large real-world
    datasets.

The examples use the same `GeneA` / `TxA1` region used throughout the
other GeneTrackR modules.

### Step 1. Prepare browser-track inputs

Locate the bundled demo files:

\`\`\`{r} gp_file \<- system.file( “extdata”, “gtr_demo.genePredExt”,
package = “GeneTrackR” )

signal_files \<- system.file( “extdata”, c(
“gtr_demo_rnaseq_plus.bedgraph”, “gtr_demo_rnaseq_minus.bedgraph”,
“gtr_demo_riboseq_plus.bedgraph”, “gtr_demo_riboseq_minus.bedgraph” ),
package = “GeneTrackR” )

feature_file \<- system.file( “extdata”, “gtr_demo_features.bed”,
package = “GeneTrackR” )

vcf_file \<- system.file( “extdata”, “gtr_demo_variants.vcf”, package =
“GeneTrackR” )


    Read the annotation, four strand-specific signal tracks, BED features, and VCF variants:

    ```{r}
    gp <- read_genepred(
      gp_file,
      format = "genePredExt",
      verbose = FALSE,
      progress = FALSE
    )

    signal_all <- read_bwg(
      signal_files,
      format = "bedgraph",
      sample_names = c(
        "RNA_seq_plus",
        "RNA_seq_minus",
        "Ribo_seq_plus",
        "Ribo_seq_minus"
      ),
      strand = c("+", "-", "+", "-"),
      mode = "memory",
      verbose = FALSE
    )

    features <- read_bed(
      feature_file,
      verbose = FALSE,
      progress = FALSE
    )

    variants <- read_vcf(
      vcf_file,
      mode = "memory",
      keep_genotype = TRUE,
      verbose = FALSE,
      progress = FALSE
    )

The four inputs play different roles in
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md):

- `annotation`: gene/transcript/exon structure used for the gene-model
  panel;
- `signal`: optional `BwgTrack` containing one or more quantitative
  tracks;
- `features`: optional `FeatureTrack` containing BED/GFF/GTF-style
  intervals;
- `variants`: optional `VariantTrack` containing VCF-derived variants.

General GeneTrackR plotting palette arguments default to `"Paired"`. The
deliberate exception is
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md),
whose `color_palette` defaults to the sequential `"Reds"` palette for LD
heatmaps. Supply another RColorBrewer palette explicitly when a
different color scheme is required.

[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
accepts **exactly one genomic locator** per call:

- `gene_id`;
- `transcript_id`;
- or the complete `chrom` + `start` + `end` combination.

Do not supply more than one locator type in the same call.

### Step 2. Build a gene-centered browser view

Start with the annotation alone. When `gene_id` is supplied, GeneTrackR
derives the chromosome and plotting interval from the complete gene
locus:

\`\`\`{r} p_browser_gene_model \<- plot_tracks( annotation = gp, gene_id
= “GeneA”, gene_palette = “Paired”, direction_mode = “transcript” )

p_browser_gene_model


    Add RNA-seq and Ribo-seq signal over the same locus. `GeneA` is on the positive strand, so select the two positive-strand samples explicitly:

    ```{r}
    p_browser_gene <- plot_tracks(
      annotation = gp,
      signal = signal_all,
      gene_id = "GeneA",
      samples = c("RNA_seq_plus", "Ribo_seq_plus"),
      strand = "+",
      signal_type = "bar",
      signal_palette = "Paired",
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

    p_browser_gene

For discrete signal tracks, colors follow the selected sample/group
order and the standard palette order. Unlike
[`plot_signal_gene()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_gene.md),
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
does not use `strand = "auto"`. Use `+` or `-` when a strand-specific
locus is known, `both` when both strand classes should be shown, or
`ignore` when signal strand metadata should not be used for filtering.

### Step 3. Build a transcript-centered browser view

Use `transcript_id` when the plotting interval should correspond to one
transcript rather than the complete gene locus:

\`\`\`{r} p_browser_transcript \<- plot_tracks( annotation = gp, signal
= signal_all, transcript_id = “TxA1”, samples = c(“RNA_seq_plus”,
“Ribo_seq_plus”), strand = “+”, signal_type = “bar”, ribo_signal_type =
“auto”, frame_palette = “Paired”, signal_y_scale = “free”, collapse =
“none”, direction_mode = “transcript”, heights = c( signal = 4, gene =
1, feature = 0.8, variant = 0.7 ) )

p_browser_transcript


    When `transcript_id` is supplied, `ribo_signal_type = "auto"` detects Ribo-seq/RPF-like sample IDs and renders those panels in `frame` mode by default, while RNA-seq panels remain standard browser signal tracks. This makes the difference between RNA-seq coverage and Ribo-seq frame periodicity more obvious without requiring a separate manual step.

    The signal panel remains in **genomic coordinates** so that it stays aligned with the other browser tracks. For region-level browser plots, Ribo-seq defaults back to standard genomic signal tracks because frame rendering requires one specific transcript model.

    ### Step 4. Inspect an explicit genomic region

    For a multi-gene interval or a locus not defined by one annotation ID, provide explicit genomic coordinates. The demo interval below contains `GeneA`, `GeneB`, and `GeneC` and therefore includes both positive- and negative-strand signal:

    ```{r}
    browser_chrom <- "chr1"
    browser_start <- 12339001
    browser_end <- 12374500

    p_browser_region <- plot_tracks(
      annotation = gp,
      signal = signal_all,
      chrom = browser_chrom,
      start = browser_start,
      end = browser_end,
      samples = c(
        "RNA_seq_plus",
        "RNA_seq_minus",
        "Ribo_seq_plus",
        "Ribo_seq_minus"
      ),
      strand = "both",
      signal_type = "bar",
      signal_palette = "Paired",
      signal_y_scale = "free",
      signal_y_ticks = "pretty",
      collapse = "none",
      direction_mode = "gene",
      heights = c(
        signal = 5,
        gene = 1.2,
        feature = 0.8,
        variant = 0.7
      )
    )

    p_browser_region

Use `signal_y_scale = "free"` when assays have substantially different
count ranges, as in the RNA-seq and Ribo-seq demo. Use `fixed` when the
signal panels should share the same y-axis range for direct quantitative
comparison.

### Step 5. Add BED features and VCF variants

Pass `FeatureTrack` and `VariantTrack` objects directly to `features`
and `variants`.
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
retrieves only records overlapping the selected browser interval:

\`\`\`{r} p_browser_complete \<- plot_tracks( annotation = gp, signal =
signal_all, features = features, variants = variants, chrom =
browser_chrom, start = browser_start, end = browser_end, samples =
c(“RNA_seq_plus”, “Ribo_seq_plus”), strand = “+”, signal_type = “bar”,
signal_palette = “Paired”, signal_y_scale = “free”, signal_y_ticks =
“pretty”, gene_palette = “Paired”, feature_color_by = “auto”,
feature_max_legend_levels = 5, variant_palette = “Paired”,
variant_colors = c( SNP = “#1F78B4”, INS = “#33A02C”, DEL = “#E31A1C” ),
direction_mode = “gene”, heights = c( signal = 4, gene = 1.2, feature =
0.9, variant = 0.8 ) )

p_browser_complete


    For BED-like interval tracks, `feature_color_by = "auto"` tries to choose a compact and informative attribute for the legend. In the bundled demo BED file, names such as `GeneA_promoter|promoter` and `GeneB_enhancer|enhancer` are simplified to a small set of feature groups (for example promoter, enhancer, candidate region, and QTL-like intervals), and the legend is automatically capped to at most five groups. This is usually more informative than a single `BED` legend entry, while still keeping the browser legend manageable.

    Variant colors are controlled independently with `variant_palette` and `variant_colors`. `variant_palette = "Paired"` is the default. For stable biological categories, a named vector such as `c(SNP = ..., INS = ..., DEL = ...)` is preferable because the same variant type then retains the same color across browser regions.

    The combined figure uses the following default panel order when `layout = "signal_top"`:

    1. signal;
    2. feature track(s);
    3. variant track(s);
    4. gene model.

    With `layout = "gene_top"`, the gene model is moved to the top, followed by signal, feature, and variant panels.

    `features` and `variants` can also be named lists when several independent interval or variant collections should be drawn as separate panels:

    ```{r}
    # feature_tracks <- list(
    #   Regulatory = regulatory_features,
    #   QTL = qtl_features
    # )
    #
    # variant_tracks <- list(
    #   Natural = natural_variants,
    #   Candidate = candidate_variants
    # )
    #
    # plot_tracks(
    #   annotation = gp,
    #   features = feature_tracks,
    #   variants = variant_tracks,
    #   chrom = browser_chrom,
    #   start = browser_start,
    #   end = browser_end
    # )

Each list entry becomes a separate panel. The `feature` and `variant`
entries of `heights` are reused for every panel of the corresponding
class.

### Step 6. Control layout, heights, highlighting, and style

Panel proportions are controlled through the named `heights` vector.
Increase `signal` when quantitative tracks need more vertical space, or
increase `gene`, `feature`, and `variant` when structural tracks are
crowded.

The following example moves the gene model to the top and highlights the
small LD demonstration interval shared with the LD workflow:

\`\`\`{r} ld_highlight \<- data.frame( start = 12342620, end = 12343180
)

p_browser_styled \<- plot_tracks( annotation = gp, signal = signal_all,
features = features, variants = variants, chrom = browser_chrom, start =
browser_start, end = browser_end, samples = c(“RNA_seq_plus”,
“Ribo_seq_plus”), strand = “+”, signal_type = “bar”, signal_palette =
“Paired”, signal_transform = “sqrt”, signal_y_scale = “free”,
signal_y_ticks = “pretty”, signal_alpha = 0.85, signal_bar_width = 0.9,
gene_palette = “Paired”, gene_border_color = “black”, feature_color_by =
“auto”, feature_max_legend_levels = 5, highlight = ld_highlight, layout
= “gene_top”, heights = c( signal = 4, gene = 1.5, feature = 0.9,
variant = 0.8 ), plot_theme = “classic”, show_panel_border = FALSE,
text_size = 14 )

p_browser_styled


    `highlight` uses genomic coordinates and is currently applied to the signal and gene-model panels. It is useful for marking a candidate interval while preserving the common browser x-axis.

    Gene-model styling is controlled by `gene_palette`, `gene_colors`, `gene_border_color`, `cds_height`, `utr_height`, `direction_mode`, `label_position`, and `label_by`. Signal styling is controlled by the corresponding `signal_*` arguments. Feature styling is controlled by `feature_color_by`, `feature_palette`, `feature_colors`, `feature_border_color`, and `feature_max_legend_levels`. Variant styling in the integrated browser is controlled by `variant_palette` and `variant_colors`.

    The integrated browser intentionally uses simplified feature and variant rendering with labels hidden. When detailed interval labels, variant labels, or more detailed track-specific customization are required, first inspect them with `plot_feature_track()`, `plot_variant()`, or `plot_signal_transcript()`, then use `plot_tracks()` for the compact integrated overview.

    ### Step 7. Use browser tracks efficiently on real datasets

    For large datasets, avoid loading genome-wide signal or VCF records when regional access is available:

    - use lazy `BwgTrack` access for indexed bigWig or supported large signal inputs;
    - use `read_vcf(mode = "lazy")` with a bgzip-compressed, indexed VCF for regional variant retrieval;
    - keep the browser interval biologically focused rather than plotting an entire chromosome with base-level Ribo-seq or dense variants;
    - preselect relevant signal samples with `samples` before plotting;
    - use dedicated functions such as `plot_signal_transcript()`, `plot_variant()`, or `plot_feature_track()` when one track requires detailed analysis beyond the compact browser representation.

    `plot_tracks()` returns the assembled ggplot/patchwork figure directly. Figure export is covered in the export workflow rather than duplicated here.

    ## Haplotype construction and visualization

    GeneTrackR builds haplotypes directly from VCF genotypes and keeps the complete analysis state in a `HapVariant` object. The recommended workflow is to define the biological region first, construct haplotypes once, inspect the resulting groups, and then pass the same object to phenotype association or haplotype-refinement functions.

    This module follows one continuous workflow:

    1. prepare VCF genotype data and gene annotation;
    2. build the primary gene-body haplotypes;
    3. inspect the `HapVariant` object and its core tables;
    4. choose an appropriate genotype representation;
    5. control missing-genotype filtering;
    6. construct transcript, flanking-region, or custom-region haplotypes;
    7. draw the haplotype-variant figure;
    8. pass the haplotype object to downstream phenotype and refinement workflows.

    The main example uses the designed `GeneA` locus. Its gene body contains 11 variants and 36 complete samples that form four balanced haplotypes with nine samples per haplotype.

    ### Step 1. Prepare VCF genotypes and annotation

    Locate and read the bundled demo files:

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

    vcf <- read_vcf(
      vcf_file,
      mode = "memory",
      keep_genotype = TRUE,
      verbose = FALSE,
      progress = FALSE
    )

    gp <- read_genepred(
      gp_file,
      format = "genePredExt",
      verbose = FALSE,
      progress = FALSE
    )

Haplotype construction requires VCF sample genotype columns. If a VCF
was previously read with `keep_genotype = FALSE`, read it again with
`keep_genotype = TRUE` before calling the haplotype functions.

For indexed large VCF files, the input can also be a lazy `VariantTrack`
or a VCF file path. GeneTrackR retrieves only the requested haplotype
region before constructing the genotype matrix.

### Step 2. Build the primary GeneA haplotypes

For gene- or transcript-defined analysis, use
[`hap_gene_variant()`](https://renscq.github.io/GeneTrackR/reference/hap_gene_variant.md).
The primary GeneA workflow intentionally uses the **gene body only** so
that the haplotype definition is based on variants inside the annotated
locus:

\`\`\`{r} hap_gene \<- hap_gene_variant( vcf = vcf, annotation = gp,
gene_id = “GeneA”, genotype_mode = “string” )

hap_gene


    The demo is designed so that this call contains:

    - 11 GeneA gene-body variants;
    - 36 retained samples;
    - four haplotype groups;
    - nine samples in each haplotype group.

    The four genotype patterns are intentionally hierarchical: `Hap1/Hap2` are the closest pair, `Hap3/Hap4` are the second closest pair, and cross-pair comparisons differ at many more GeneA variants. This structure is used later by the phenotype-refinement demo.

    Check these values directly from the returned object rather than assuming a fixed number for real datasets:

    ```{r}
    hap_gene$meta
    hap_gene$haplotypes[, c("hap_id", "sample_n", "samples"), with = FALSE]

`hap_id` is assigned from the observed genotype patterns after grouping.
Treat the genotype pattern and sample membership as the biological
result; do not assume that `Hap1`, `Hap2`, and so on have the same
biological meaning across independently constructed regions.

### Step 3. Inspect the HapVariant object

A `HapVariant` keeps the complete state needed by downstream GeneTrackR
functions:

`{r} hap_gene$region hap_gene$variants hap_gene$genotype_long hap_gene$genotype_wide hap_gene$haplotypes hap_gene$sample_haplotypes`

The components have distinct purposes:

| Component | Content | Typical use |
|----|----|----|
| `region` | locator, chromosome, boundaries, strand, and flank settings | confirm the biological interval |
| `variants` | retained VCF variant records | inspect variant order and REF/ALT alleles |
| `genotype_long` | one sample × variant genotype per row | detailed genotype inspection |
| `genotype_wide` | one sample per row with variant columns | inspect the full genotype matrix |
| `haplotypes` | one row per unique haplotype with sample counts | summarize haplotype groups |
| `sample_haplotypes` | sample-to-haplotype assignment | join with phenotype or other metadata |

The `meta` list records the genotype representation, missing-value
policy, retained sample number, variant number, and haplotype number:

`{r} hap_gene$meta`

This object is already the canonical input for
[`plot_hap_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_variant.md),
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md),
[`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md),
and
[`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md).
There is no need to rebuild a genotype matrix manually.

### Step 4. Choose the genotype representation

GeneTrackR provides two genotype representations through
`genotype_mode`.

#### Allele-string mode

`genotype_mode = "string"` is recommended when the haplotype table
itself is important for interpretation:

\`\`\`{r} hap_gene_string \<- hap_gene_variant( vcf = vcf, annotation =
gp, gene_id = “GeneA”, genotype_mode = “string” )

hap_gene_string\$haplotypes


    Single-base alleles are displayed directly as `A`, `T`, `C`, or `G`. Longer REF/ALT alleles are compressed to `iN`, where `N` is the displayed allele length. The haplotype plot maps allele classes in the stable order `A`, `T`, `C`, `G`, and `indel`.

    String mode is also a compact display representation rather than a phased dosage representation: if a genotype contains an alternate allele, GeneTrackR displays the first observed ALT allele label for that variant. Inspect the original genotype columns in `vcf$data` when zygosity, phase, or multi-allelic dosage must be retained explicitly.

    #### Binary code mode

    `genotype_mode = "code"` collapses each genotype to reference versus alternate presence:

    ```{r}
    hap_gene_code <- hap_gene_variant(
      vcf = vcf,
      annotation = gp,
      gene_id = "GeneA",
      genotype_mode = "code"
    )

    hap_gene_code$haplotypes

In code mode:

- `0` means the genotype contains only the reference allele;
- `1` means at least one alternate allele is present.

Therefore code mode is a **binary REF/ALT-presence representation**, not
allele dosage. A heterozygous genotype and a homozygous alternate
genotype both become `1`. Use allele-string mode when the displayed
allele identity is more important than compact binary grouping.

### Step 5. Control missing-genotype filtering

By default, `min_variant_number = NULL` requires a sample to have
non-missing genotypes at **all retained variants**. This is the safest
default for defining complete haplotypes.

The GeneA gene body has complete demo genotypes, so all 36 samples are
retained:

`{r} hap_gene$meta$sample_n hap_gene$meta$min_variant_number`

Flanking regions can introduce additional variants with missing data.
The demo upstream variant `varAup01` is intentionally designed with
missing/heterozygous genotypes. Adding it changes the haplotype
definition:

\`\`\`{r} hap_gene_upstream \<- hap_gene_variant( vcf = vcf, annotation
= gp, gene_id = “GeneA”, upstream = 500, downstream = 0, strand_aware =
TRUE, genotype_mode = “string” )

hap_gene_upstream$`region
hap_gene_upstream`$meta


    To retain samples that are missing a limited number of variants, set an explicit threshold. The extended demo region contains 12 variants, so `min_variant_number = 11` allows one missing genotype:

    ```{r}
    hap_gene_upstream_relaxed <- hap_gene_variant(
      vcf = vcf,
      annotation = gp,
      gene_id = "GeneA",
      upstream = 500,
      genotype_mode = "string",
      min_variant_number = 11
    )

    hap_gene_upstream_relaxed$sample_haplotypes

Use relaxed thresholds deliberately. Missing values become part of the
observed haplotype pattern and may create additional haplotype groups.
For the main GeneA tutorial, the unextended gene-body object `hap_gene`
remains the recommended analysis object.

A custom missing label can be supplied for display:

`{r} hap_missing_label <- hap_gene_variant( vcf = vcf, annotation = gp, gene_id = "GeneA", upstream = 500, genotype_mode = "code", missing_genotype = "-", min_variant_number = 11 )`

Internally, GeneTrackR still tracks which genotypes are biologically
missing, so a custom display label does not convert missing values into
valid genotype calls.

### Step 6. Construct transcript, flanking-region, and custom-region haplotypes

#### Transcript-defined haplotypes

Use `transcript_id` instead of `gene_id` when variants should be
restricted to one transcript locus:

\`\`\`{r} hap_tx \<- hap_gene_variant( vcf = vcf, annotation = gp,
transcript_id = “TxA1”, genotype_mode = “string” )

hap_tx\$region


    Supply exactly one of `gene_id` or `transcript_id`.

    #### Strand-aware flanking regions

    `upstream` and `downstream` are interpreted relative to transcription direction when `strand_aware = TRUE`:

    ```{r}
    hap_gene_flank <- hap_gene_variant(
      vcf = vcf,
      annotation = gp,
      gene_id = "GeneA",
      upstream = 1000,
      downstream = 500,
      strand_aware = TRUE,
      genotype_mode = "string",
      min_variant_number = 1
    )

    hap_gene_flank$region

For a negative-strand gene, biological upstream is toward increasing
genomic coordinates. GeneTrackR handles this reversal automatically when
`strand_aware = TRUE`.

#### Explicit genomic regions

Use
[`hap_region_variant()`](https://renscq.github.io/GeneTrackR/reference/hap_region_variant.md)
when the interval is defined directly rather than by annotation:

\`\`\`{r} hap_region \<- hap_region_variant( vcf = vcf, chrom = “chr1”,
start = 12340001, end = 12352000, genotype_mode = “code” )

hap_region


    This is useful for GWAS/LD intervals, promoter windows, QTL regions, or any custom locus that does not correspond exactly to one annotated gene or transcript.

    #### Sample and variant-type filtering

    The haplotype can be restricted before grouping. For example, analyze a subset of samples:

    ```{r}
    hap_sample_subset <- hap_gene_variant(
      vcf = vcf,
      annotation = gp,
      gene_id = "GeneA",
      samples = sprintf("S%02d", 1:18),
      genotype_mode = "string"
    )

    hap_sample_subset$meta

Or construct haplotypes using only SNPs:

\`\`\`{r} hap_snp \<- hap_gene_variant( vcf = vcf, annotation = gp,
gene_id = “GeneA”, variant_type = “SNP”, genotype_mode = “string” )

hap_snp\$variants


    Filtering variants changes the haplotype definition and may merge previously distinct haplotypes. Always inspect `haplotypes` and `sample_haplotypes` after changing the retained variant set.

    `hap_variant()` remains available as a compatibility wrapper, but new code should prefer `hap_gene_variant()` for gene/transcript queries and `hap_region_variant()` for direct genomic intervals.

    ### Step 7. Draw the haplotype-variant figure

    `plot_hap_variant()` combines the gene model, natural-variant markers, connector lines, REF/ALT reference rows, and haplotype genotype table:

    ```{r}
    hap_variant_figure <- plot_hap_variant(
      hap_gene,
      annotation = gp,
      min_hap_samples = 3,
      show_reference_row = TRUE,
      variant_label = "pos",
      gene_pos_x_angle = 90,
      gene_track_legend_position = "top",
      direction_mode = "gene",
      table_x_angle = 90
    )

    hap_variant_figure

General plotting palettes default to `Paired`;
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
is the deliberate package-level exception and defaults to `Reds`. The
main haplotype visual controls can still be changed independently:

\`\`\`{r} hap_variant_custom \<- plot_hap_variant( hap_gene, annotation
= gp, min_hap_samples = 3, show_reference_row = TRUE, gene_palette =
“Paired”, table_palette = “Paired”, variant_palette = “Paired”,
table_alpha = 0.7, variant_alpha = 0.8, genotype_text_size = 3,
variant_marker_size = 3, gene_track_height = 1.4, connector_height = 0.4
)

hap_variant_custom


    For allele-string haplotypes, table colors always follow the stable biological class order `A`, `T`, `C`, `G`, and `indel`. Custom colors can be supplied explicitly with those names:

    ```{r}
    plot_hap_variant(
      hap_gene,
      annotation = gp,
      min_hap_samples = 3,
      table_colors = c(
        A = "#1B9E77",
        T = "#D95F02",
        C = "#7570B3",
        G = "#E7298A",
        indel = "#66A61E"
      )
    )

Use `show_variant_marker = FALSE` to hide natural-variant triangles, or
adjust `variant_marker_size` when the locus contains many variants.

### Step 8. Pass the same haplotype object downstream

Do not reconstruct haplotypes separately for each downstream analysis.
The `hap_gene` object created in Step 2 is the input for the following
modules:

\`\`\`{r} \# Haplotype-phenotype association (09-phenotype.qmd) \#
hap_res \<- plot_hap_pheno( \# hap_gene, \# phenotype = pheno, \# traits
= “seed_weight”, \# min_hap_samples = 3 \# )

# Phenotype-guided refinement (11-haplotype-refinement.qmd)

# refined \<- refine_haplotype(

# hap_gene,

# phenotype = pheno,

# traits = “seed_weight”,

# min_hap_samples = 3

# )

# Variant-effect prioritization (12-variant-effect.qmd)

# effect \<- plot_variant_effect(

# hap_gene,

# phenotype = pheno,

# traits = “protein_content”,

# min_group_samples = 3

# )


    Keeping one well-defined `HapVariant` as the common input ensures that the variant set, sample set, missing-data rule, and haplotype assignment remain consistent across visualization, phenotype association, refinement, and variant-effect analysis.

    ## Phenotype association workflow

    GeneTrackR keeps phenotype data as a simple sample-level table and joins it to haplotype or variant genotypes by `sample_id`. The recommended workflow is to inspect phenotype quality first, verify sample matching explicitly, then perform haplotype- or variant-based association tests and keep the returned statistical tables together with the figure.

    This module follows one continuous workflow:

    1. prepare phenotype, VCF, annotation, and the standard GeneA haplotypes;
    2. summarize phenotype types, missingness, and distributions;
    3. verify sample matching before association analysis;
    4. test haplotype-phenotype associations;
    5. inspect multiple traits and statistical-test settings;
    6. control phenotype-association plot appearance;
    7. test a single variant against phenotype;
    8. pass the same phenotype and haplotype objects to refinement and variant-effect workflows.

    The bundled phenotype data contain four numeric traits and one categorical trait. Their designed roles are:

    | Trait | Type | Demo role |
    | --- | --- | --- |
    | `seed_weight` | numeric | strong four-level GeneA haplotype effect with nearest genotype pairs remaining phenotypically closer |
    | `plant_height` | numeric | moderate four-level GeneA haplotype effect with the same hierarchical ordering |
    | `protein_content` | numeric | strong `varA03` effect aligned with the `Hap1/Hap2` versus `Hap3/Hap4` genotype split |
    | `flowering_time` | numeric | negative control for GeneA haplotypes |
    | `flower_color` | categorical | categorical phenotype example |

    ### Step 1. Prepare phenotype and genotype inputs

    Locate the bundled phenotype, VCF, and annotation files:

    ```{r}
    pheno_file <- system.file(
      "extdata",
      "gtr_demo_pheno.tsv",
      package = "GeneTrackR"
    )

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

Read the three inputs and reconstruct the same GeneA gene-body
haplotypes used in the haplotype workflow:

\`\`\`{r} pheno \<- read_pheno( pheno_file, verbose = FALSE, progress =
FALSE )

vcf \<- read_vcf( vcf_file, mode = “memory”, keep_genotype = TRUE,
verbose = FALSE, progress = FALSE )

gp \<- read_genepred( gp_file, format = “genePredExt”, verbose = FALSE,
progress = FALSE )

hap_gene \<- hap_gene_variant( vcf = vcf, annotation = gp, gene_id =
“GeneA”, genotype_mode = “string” )


    `read_pheno()` returns a `data.table`. By default, the first column is treated as the sample column and renamed to `sample_id`. If another column stores sample IDs, specify it explicitly:

    ```{r}
    # pheno <- read_pheno(
    #   "phenotype.tsv",
    #   sample_col = "Taxa"
    # )

Sample IDs must be unique. Duplicate IDs are rejected because they would
make genotype-phenotype matching ambiguous.

### Step 2. Summarize phenotype quality and distributions

Start with
[`summary_pheno()`](https://renscq.github.io/GeneTrackR/reference/summary_pheno.md)
before any association test:

`{r} pheno_summary <- summary_pheno(pheno) pheno_summary`

The summary reports:

- numeric versus categorical trait type;
- total sample number;
- missing number and missing rate;
- number of unique non-missing values;
- minimum, mean, median, and maximum for numeric traits.

Plot numeric and categorical phenotype distributions directly from the
phenotype table:

\`\`\`{r} pheno_figure \<- plot_pheno( pheno, traits = c( “seed_weight”,
“protein_content”, “plant_height”, “flowering_time”, “flower_color” ) )

pheno_figure


    `plot_pheno()` returns the figure directly. Numeric traits are shown as histograms, while categorical traits are shown as count bars. When both trait classes are requested, the result is a combined patchwork figure.

    `plot_hap_pheno()` and `plot_variant_pheno()` are currently intended for **numeric phenotype traits**. Categorical traits such as `flower_color` should therefore be summarized or plotted with `summary_pheno()` / `plot_pheno()` unless a categorical association method is implemented separately.

    ### Step 3. Verify sample matching explicitly

    Phenotype rows do not need to follow the VCF or haplotype sample order. GeneTrackR joins genotype and phenotype records by sample ID.

    The demo phenotype table intentionally uses a different row order from the genotype data. Check the overlap before testing:

    ```{r}
    hap_sample_ids <- as.character(hap_gene$sample_haplotypes$sample_id)
    pheno_sample_ids <- as.character(pheno$sample_id)

    sample_match <- data.frame(
      phenotype_samples = length(unique(pheno_sample_ids)),
      haplotype_samples = length(unique(hap_sample_ids)),
      matched_samples = length(intersect(pheno_sample_ids, hap_sample_ids)),
      phenotype_only = length(setdiff(pheno_sample_ids, hap_sample_ids)),
      haplotype_only = length(setdiff(hap_sample_ids, pheno_sample_ids))
    )

    sample_match

For the bundled demo, all 36 phenotype samples match the 36 GeneA
haplotype samples even though their row orders differ.

In real datasets, samples present on only one side are not available for
the association test. Checking
[`setdiff()`](https://rdrr.io/r/base/sets.html) before plotting is
therefore recommended, especially after genotype missingness filters or
phenotype QC have removed samples.

### Step 4. Test the primary haplotype-phenotype association

`seed_weight` is the main deterministic GeneA haplotype phenotype in the
demo. Test the four GeneA haplotypes with
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md):

\`\`\`{r} hap_seed \<- plot_hap_pheno( hap = hap_gene, phenotype =
pheno, traits = “seed_weight”, min_hap_samples = 3, plot_type =
“violin_boxplot”, test_method = “t.test”, p_adjust = “BH”, p_value_type
= “adjusted”, p_label = “stars”, show_signif_only = TRUE, fill_palette =
“Paired” )

hap_seed\$figure


    `plot_hap_pheno()` returns a `GeneTrackRPhenoPlot`, not only a ggplot. Keep the complete result because it contains the statistical output used to draw the figure:

    ```{r}
    hap_seed$pvalue
    hap_seed$summary
    hap_seed$bracket

The main components are:

| Component | Content |
|----|----|
| `$figure` | phenotype distribution figure |
| `$pvalue` | all pairwise tests, including raw and adjusted p-values |
| `$summary` | sample number, median, mean, minimum, and maximum by trait and haplotype |
| `$bracket` | comparisons actually prepared for display on the figure |
| `$plot_data` | matched long-format phenotype data used for plotting |

Printing the complete `GeneTrackRPhenoPlot` also prints its `$figure`,
but use the named components when saving figures or exporting
statistical tables.

Haplotype groups are ordered by retained sample number. The x-axis label
includes the group sample count, for example `Hap1 (9)`.

### Step 5. Compare traits and choose the statistical settings

The demo provides three useful haplotype-level behaviors in one dataset:

- `seed_weight`: strong four-level differences, while the
  genotype-nearest pairs remain phenotypically closer than cross-cluster
  pairs;
- `plant_height`: the same hierarchical pattern with a more moderate
  effect scale;
- `flowering_time`: negative control without a designed GeneA haplotype
  effect.

Plot them together while preserving the requested trait order:

\`\`\`{r} hap_multi \<- plot_hap_pheno( hap = hap_gene, phenotype =
pheno, traits = c( “seed_weight”, “plant_height”, “flowering_time” ),
min_hap_samples = 3, test_method = “t.test”, p_adjust = “BH”,
p_value_type = “adjusted”, p_label = “stars”, facet_ncol = 3,
fill_palette = “Paired” )

hap_multi$`figure
hap_multi`$pvalue


    Three pairwise test methods are available:

    - `t.test`: compares group means and is appropriate for approximately continuous, reasonably behaved traits;
    - `wilcox.test`: rank-based alternative that is less dependent on normality assumptions;
    - `ks.test`: tests whether the overall distributions differ, not specifically their means or medians.

    For example:

    ```{r}
    hap_seed_wilcox <- plot_hap_pheno(
      hap = hap_gene,
      phenotype = pheno,
      traits = "seed_weight",
      min_hap_samples = 3,
      test_method = "wilcox.test",
      p_adjust = "BH",
      p_value_type = "adjusted"
    )

    hap_seed_wilcox$pvalue

`$pvalue` always retains both `p_value` and `p_adj`. The argument
`p_value_type` controls which of these is used to filter and label
significance brackets; it does not remove the other statistic from the
returned table.

With four haplotypes there are six pairwise comparisons per trait, so
adjusted p-values are generally preferable for the main interpretation.
`p_adjust = "BH"` is the default multiple-testing adjustment used in
these examples.

`show_signif_only = TRUE` only controls bracket display. Non-significant
comparisons remain available in `$pvalue`.

### Step 6. Control phenotype-association plot appearance

General GeneTrackR plotting palettes default to `Paired`;
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
is the deliberate package-level exception and defaults to `Reds`.
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md)
maps fill colors according to the **median phenotype value within each
trait**, so the colors are used as an ordered visual encoding rather
than as permanent haplotype identities.

A more detailed plot can be requested as follows:

\`\`\`{r} hap_seed_styled \<- plot_hap_pheno( hap = hap_gene, phenotype
= pheno, traits = “seed_weight”, min_hap_samples = 3, plot_type =
“violin_boxplot”, test_method = “t.test”, p_adjust = “BH”, p_value_type
= “adjusted”, p_label = “both”, show_points = TRUE, show_outliers =
FALSE, fill_palette = “Paired”, fill_alpha = 0.75, x_text_angle = 90,
text_size = 14 )

hap_seed_styled\$figure


    Use `plot_type = "violin"`, `"boxplot"`, or `"violin_boxplot"` depending on whether the distribution shape, compact summary, or both are most useful.

    If fixed colors for specific haplotype IDs are required, supply a named vector instead of relying on the median-based palette assignment:

    ```{r}
    hap_seed_fixed_colors <- plot_hap_pheno(
      hap = hap_gene,
      phenotype = pheno,
      traits = "seed_weight",
      min_hap_samples = 3,
      fill_colors = c(
        Hap1 = "#A6CEE3",
        Hap2 = "#1F78B4",
        Hap3 = "#B2DF8A",
        Hap4 = "#33A02C"
      )
    )

    hap_seed_fixed_colors$figure

Use fixed haplotype colors only after confirming the actual `hap_id`
values in `hap_gene$haplotypes`; haplotype numbering is specific to the
variant set used to construct each `HapVariant`.

### Step 7. Test a single variant against phenotype

[`plot_variant_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_pheno.md)
is the single-variant counterpart of
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md).
The demo variant `varA03` follows the same genotype-cluster split used
by refinement: `Hap1/Hap2` carry the REF class, `Hap3/Hap4` carry the
ALT class, and ALT carriers have higher `protein_content` values.

Select the variant by ID and use the compact binary REF/ALT
representation:

\`\`\`{r} variant_protein \<- plot_variant_pheno( variant = vcf,
phenotype = pheno, variant_id = “varA03”, traits = “protein_content”,
genotype_mode = “code”, min_group_samples = 3, test_method = “t.test”,
p_adjust = “BH”, p_value_type = “adjusted”, p_label = “stars”,
fill_palette = “Paired” )

variant_protein$`figure
variant_protein`$pvalue variant_protein$`summary
variant_protein`$variant_data


    For `genotype_mode = "code"`:

    - `0` means reference-only genotype;
    - `1` means at least one alternate allele is present.

    This is a binary REF/ALT-presence representation, not allele dosage. Use `genotype_mode = "string"` when allele labels are more informative:

    ```{r}
    variant_protein_string <- plot_variant_pheno(
      variant = vcf,
      phenotype = pheno,
      variant_id = "varA03",
      traits = "protein_content",
      genotype_mode = "string",
      min_group_samples = 3,
      p_value_type = "adjusted"
    )

    variant_protein_string$figure

A variant can also be selected by exact genomic position:

\`\`\`{r} variant_protein_position \<- plot_variant_pheno( variant =
vcf, phenotype = pheno, chrom = “chr1”, pos = 12342550, traits =
“protein_content”, genotype_mode = “code”, min_group_samples = 3,
p_value_type = “adjusted” )

variant_protein_position\$variant_data


    The selected locator must resolve to exactly one variant. If a region contains several variants, provide `variant_id` or an exact chromosome-position locator.

    `plot_variant_pheno()` returns the same `GeneTrackRPhenoPlot` structure as `plot_hap_pheno()`, with one additional component: `$variant_data`, containing the selected variant record.

    ### Step 8. Pass the same objects downstream

    Keep the same `hap_gene` and `pheno` objects for phenotype-guided refinement and variant-effect prioritization. Rebuilding the haplotypes with a different region or missing-data rule would change the biological groups being tested.

    ```{r}
    # Haplotype refinement (11-haplotype-refinement.qmd)
    # refined <- refine_haplotype(
    #   hap_gene,
    #   phenotype = pheno,
    #   traits = "protein_content",
    #   min_hap_samples = 3,
    #   effect_threshold = 0.5
    # )

    # Variant-effect prioritization (12-variant-effect.qmd)
    # effect <- plot_variant_effect(
    #   hap_gene,
    #   phenotype = pheno,
    #   traits = "protein_content",
    #   min_group_samples = 3
    # )

The resulting workflow keeps one consistent chain:

`VCF + annotation -> HapVariant -> phenotype association -> refinement / variant-effect analysis`.

For reproducible analysis, export `$pvalue`, `$summary`, and the matched
sample-level `$plot_data` together with `$figure` rather than saving
only the visual result.

## Linkage disequilibrium analysis and visualization

GeneTrackR separates linkage disequilibrium analysis into two stages.
[`compute_ld_block()`](https://renscq.github.io/GeneTrackR/reference/compute_ld_block.md)
calculates pairwise LD from VCF genotypes and stores the complete result
in an `LDTrack`;
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
then adds a triangular heatmap and optionally a gene/variant region
track to the same object. Keeping calculation and plotting separate
makes the numerical LD table, matrix, region definition, sample
settings, and final figure available for downstream analysis.

The bundled demo now contains three deterministic LD scenarios:

- `varLD01`-`varLD12`: twelve SNPs forming the primary **LD-gradient**
  example. `varLD01`-`varLD06` are a perfect high-LD core, while
  `varLD07`-`varLD12` progressively decorrelate from that core so the
  heatmap contains high, intermediate, and low `r2` values;

- `varLow01`-`varLow04`: four variants designed to have little linkage,
  with most pairwise `r2 = 0` and a maximum of approximately `0.111`;

- `varPair01`-`varPair02`: exactly two variants, retained to test the
  single-diamond LD plot.

The added gradient variants are outside the `GeneA` gene body, so the
canonical GeneA haplotype example still contains 11 variants and four
balanced haplotypes.

This module follows nine steps:

1.  prepare genotype-rich VCF data and gene annotation;

2.  compute the primary twelve-variant LD-gradient example;

3.  inspect the `LDTrack` object and pairwise statistics;

4.  inspect the LD matrix and quantify high/intermediate/low LD pairs;

5.  draw the triangular LD heatmap using the default `Reds` palette;

6.  add gene structure, natural-variant markers, and connector lines;

7.  compare the LD-gradient, low-LD, and two-variant regions;

8.  control samples, variant types, LD methods, and retained genotype
    matrices;

9.  interpret LD results and pass a selected region to haplotype
    analysis.

### Step 1. Prepare VCF genotypes and annotation

Locate and read the same demo VCF and GenePred annotation used in the
variant and haplotype workflows:

\`\`\`{r}

vcf_file \<- system.file(

“extdata”,

“gtr_demo_variants.vcf”,

package = “GeneTrackR”

)

gp_file \<- system.file(

“extdata”,

“gtr_demo.genePredExt”,

package = “GeneTrackR”

)

vcf \<- read_vcf(

vcf_file,

mode = “memory”,

keep_genotype = TRUE,

verbose = FALSE,

progress = FALSE

)

gp \<- read_genepred(

gp_file,

format = “genePredExt”,

verbose = FALSE,

progress = FALSE

)


    LD calculation requires VCF genotype columns. If the VCF was previously read with `keep_genotype = FALSE`, read it again with `keep_genotype = TRUE` before calling `compute_ld_block()`.

    For large bgzip-compressed and indexed VCF files, a lazy `VariantTrack` can be used, but lazy LD queries require a complete `chrom`, `start`, and `end` interval so only the requested region is retrieved.

    ### Step 2. Compute the primary twelve-variant LD-gradient example

    The main demo interval starts inside `GeneA` and extends into the short downstream interval before `GeneB`. It contains 12 designed LD SNPs plus two GeneA indels. Use `variant_type = "snp"` so the LD example contains exactly `varLD01`-`varLD12`:

    ```{r}

    ld_gradient <- compute_ld_block(

      vcf,

      chrom = "chr1",

      start = 12342620,

      end = 12355500,

      variant_type = "snp",

      method = "r2",

      min_pair_samples = 3,

      verbose = FALSE

    )

    ld_gradient

The expected result is:

- 12 retained SNPs;

- 66 unique variant pairs (`12 * 11 / 2`);

- 36 available samples per pair;

- a perfect six-variant core (`varLD01`-`varLD06`, `r2 = 1` within the
  core);

- a progressive decline in LD toward `varLD12`, including weak pairs
  with `r2 < 0.2` and pairs reaching `r2 = 0`.

Check the retained variants before interpreting the heatmap:

\`\`\`{r}

ld_gradient\$variants


    The first six variants remain inside GeneA and use the same genotype pattern required by the GeneA haplotype demo. The six added variants lie downstream of the GeneA gene body, so they increase LD-plot density without changing the 11-variant GeneA haplotype definition.

    ### Step 3. Inspect the LDTrack object and pairwise statistics

    `LDTrack` keeps the complete state of an LD calculation:

    ```{r}

    names(ld_gradient)

The main components are:

Component \| Content \| Typical use \|

— \| — \| — \|

`data` \| one row per unique variant pair \| inspect `r2`, `Dprime`,
distance, sample number, and allele frequencies \|

`matrix` \| symmetric matrix of the selected LD metric \| matrix-based
inspection or export \|

`variants` \| ordered metadata for retained variants \| confirm
positions and IDs \|

`region` \| chromosome and queried boundaries \| record the analyzed
interval \|

`genotype` \| optional variant-by-sample dosage matrix \| detailed QC
when explicitly retained \|

`figure` \| figure created by
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
\| plotting/export \|

`meta` \| method, sample count, ploidy, and calculation notes \|
reproduce the analysis settings \|

Inspect the pairwise table directly:

\`\`\`{r}

ld_gradient\$data\[, .(

variant_i,

variant_j,

distance_bp,

n_samples,

r,

r2,

Dprime,

ld

)\]


    Important columns are:

    - `distance_bp`: physical distance between the two variants;

    - `n_samples`: paired non-missing sample number used for that comparison;

    - `r`: dosage correlation;

    - `r2`: squared dosage correlation;

    - `Dprime`: dosage-based D-prime approximation;

    - `ld`: the metric selected by `method`; for `method = "r2"`, `ld` is identical to `r2`.

    The region and settings used for calculation are preserved in the object:

    ```{r}

    ld_gradient$region

    ld_gradient$meta

### Step 4. Inspect the LD matrix and quantify the LD gradient

The symmetric LD matrix uses variant IDs as row and column names and
contains 1 on the diagonal:

\`\`\`{r}

round(ld_gradient\$matrix, 3)


    Unlike the previous all-identical six-variant example, the new matrix is intentionally heterogeneous. The upper-left core is strongly linked, whereas LD decreases as progressively decorrelated downstream variants are compared with the core.

    Summarize the distribution directly:

    ```{r}

    range(ld_gradient$data$r2, na.rm = TRUE)

    summary(ld_gradient$data$r2)

Define descriptive groups explicitly when useful:

\`\`\`{r}

ld_pair_summary \<- ld_gradient\$data\[, .(

high_ld = sum(r2 \>= 0.8, na.rm = TRUE),

intermediate_ld = sum(r2 \>= 0.2 & r2 \< 0.8, na.rm = TRUE),

low_ld = sum(r2 \< 0.2, na.rm = TRUE)

)\]

ld_pair_summary


    For the deterministic demo, 15 pairs have `r2 >= 0.8` and 16 pairs have `r2 < 0.2`. The remaining pairs form the intermediate gradient.

    GeneTrackR calculates pairwise LD but does **not** automatically call biological LD-block boundaries from a fixed threshold. A rule such as `r2 >= 0.8` is therefore a user-defined interpretation threshold, not an implicit block-calling rule inside `compute_ld_block()`.

    ### Step 5. Draw the triangular LD heatmap

    `plot_ld_block()` now defaults to the sequential `Reds` palette, which is well suited to an `r2` scale from 0 to 1:

    ```{r}

    ld_gradient <- plot_ld_block(

      ld_gradient,

      label_by = "variant_id",

      show_variant_labels = TRUE,

      show_region = FALSE,

      return_object = TRUE

    )

    ld_gradient$figure

The default color interpretation is straightforward: pale cells indicate
low LD and progressively darker red cells indicate stronger LD.

The palette can still be changed explicitly when needed:

\`\`\`{r}

ld_gradient_paired \<- plot_ld_block(

ld_gradient,

color_palette = “Reds”,

show_variant_labels = FALSE,

return_object = TRUE

)

ld_gradient_paired\$figure


    `plot_ld_block()` returns the updated `LDTrack` by default. The calculation results remain unchanged and the generated ggplot/patchwork figure is stored in `$figure`.

    ### Step 6. Add gene structure and natural-variant markers

    For locus interpretation, add the compact genomic region track above the triangular heatmap:

    ```{r}

    ld_gradient_region <- plot_ld_block(

      ld_gradient,

      show_region = TRUE,

      annotation = gp,

      connect_region = TRUE,

      show_variant_marker = TRUE,

      variant_marker_size = 3,

      label_by = "variant_id",

      show_variant_labels = TRUE,

      region_height = 1.25,

      connector_height = 0.35,

      return_object = TRUE

    )

    ld_gradient_region$figure

The combined plot contains three aligned components when
`connect_region = TRUE`:

1.  compact gene/variant region track;

2.  connector lines from genomic variant positions to LD columns;

3.  triangular LD heatmap.

With 12 variants, the region panel is now dense enough to show the
relationship between physical position and the LD gradient. Variant
triangle markers can be hidden with `show_variant_marker = FALSE`, or
resized with `variant_marker_size`.

### Step 7. Compare the LD-gradient, low-LD, and two-variant regions

The dedicated low-LD interval remains useful as an extreme comparison:

\`\`\`{r}

ld_low \<- compute_ld_block(

vcf,

chrom = “chr2”,

start = 1999000,

end = 2008500,

method = “r2”,

verbose = FALSE

)

ld_low\$data\[, .(

variant_i,

variant_j,

distance_bp,

r2

)\]


    The six pairwise values are designed so that most are zero and the maximum is approximately `1 / 9`:

    ```{r}

    range(ld_low$data$r2, na.rm = TRUE)

    max(ld_low$data$r2, na.rm = TRUE)

Plot it using the same default `Reds` scale:

\`\`\`{r}

ld_low \<- plot_ld_block(

ld_low,

label_by = “variant_id”,

return_object = TRUE

)

ld_low\$figure


    The two-variant demo region contains exactly one pairwise comparison:

    ```{r}

    ld_pair <- compute_ld_block(

      vcf,

      chrom = "chr2",

      start = 16995001,

      end = 17006000,

      method = "r2",

      verbose = FALSE

    )

    ld_pair$data

[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
draws this case as one centered diamond:

\`\`\`{r}

ld_pair \<- plot_ld_block(

ld_pair,

label_by = “variant_id”,

return_object = TRUE

)

ld_pair\$figure


    When only the figure is needed for a temporary plotting call, compatibility mode can return it directly:

    ```{r}

    pair_figure <- plot_ld_block(

      ld_pair,

      return_object = FALSE

    )

    pair_figure

### Step 8. Control samples, variant types, LD methods, and genotype retention

#### Select samples

LD is population dependent. Use `samples` when LD should be estimated in
a defined population subset:

\`\`\`{r}

subset_samples \<- sprintf(“S%02d”, 1:24)

ld_subset \<- compute_ld_block(

vcf,

chrom = “chr1”,

start = 12342620,

end = 12355500,

variant_type = “snp”,

samples = subset_samples,

min_pair_samples = 10,

method = “r2”,

verbose = FALSE

)

ld_subset$`meta`$sample_n


    #### Select variant types

    `variant_type` controls the records retained before pairwise calculation. It is important in the primary demo because the genomic interval also contains GeneA indels that are not part of the 12-SNP LD truth set:

    ```{r}

    ld_snp <- compute_ld_block(

      vcf,

      chrom = "chr1",

      start = 12342620,

      end = 12355500,

      variant_type = "snp",

      method = "r2",

      verbose = FALSE

    )

Available choices are `both`, `snp`, and `ind`. At least two retained
variants are required.

#### Retain the dosage matrix when needed

The genotype dosage matrix is omitted by default to reduce object size.
Retain it only when detailed genotype QC or downstream matrix analysis
is needed:

\`\`\`{r}

ld_with_gt \<- compute_ld_block(

vcf,

chrom = “chr1”,

start = 12342620,

end = 12355500,

variant_type = “snp”,

method = “r2”,

keep_genotype_matrix = TRUE,

verbose = FALSE

)

dim(ld_with_gt\$genotype)

ld_with_gt\$genotype\[, 1:6, drop = FALSE\]


    For the demo, the retained genotype matrix is `12 x 36`.

    #### Use D-prime cautiously

    `method = "Dprime"` is also available:

    ```{r}

    ld_dprime <- compute_ld_block(

      vcf,

      chrom = "chr2",

      start = 1999000,

      end = 2008500,

      method = "Dprime",

      verbose = FALSE

    )

    ld_dprime$data[, .(

      variant_i,

      variant_j,

      r2,

      Dprime,

      ld

    )]

The current `Dprime` implementation is a **dosage-based approximation
for unphased VCF genotypes**. It should not be described as an exact
haplotype-based D-prime estimate. For ordinary unphased genotype VCF
workflows, `r2` is the recommended primary GeneTrackR LD metric.

### Step 9. Interpret LD and pass the selected region downstream

LD measures statistical association between variants; it does not by
itself identify a causal variant or prove that all variants in a
visually strong region should be treated as one biological unit.
Interpret the pattern together with the focal variant, physical
distance, gene structure, population, variant quality, and
phenotype/association evidence.

Once an interval has been selected biologically, the same region can be
passed to the haplotype workflow. To reproduce the 12-SNP LD definition,
retain SNPs only:

\`\`\`{r}

hap_ld_region \<- hap_region_variant(

vcf,

chrom = “chr1”,

start = 12342620,

end = 12355500,

variant_type = “SNP”,

genotype_mode = “string”

)

hap_ld_region


    This keeps the analysis chain explicit:

    ```text

    VCF genotypes

        -> pairwise LD calculation

        -> inspect the LD gradient and regional gene context

        -> define the biological interval of interest

        -> construct regional haplotypes

        -> phenotype/refinement analysis if appropriate

For real GWAS loci, it is usually preferable to define the LD query
around a lead SNP or candidate interval first, inspect the pairwise
pattern, and only then decide which linked variants or genes should move
into the haplotype analysis.

## Haplotype refinement

[`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md)
performs **phenotype-guided grouping of existing haplotypes**. It does
not discover new variants and it does not rebuild the original genotype
patterns. Instead, it starts from an existing `HapVariant`, compares the
phenotype distributions of the original haplotypes, and groups
haplotypes that are sufficiently similar for the selected trait or
traits.

This module uses the same `GeneA` haplotypes created in the haplotype
workflow and the same phenotype table used in the phenotype workflow.
The demo is deliberately structured so that the refinement behavior can
be checked against known expectations:

- `protein_content`: the genotype-nearest pair `Hap1/Hap2` shares the
  low phenotype class and the genotype-nearest pair `Hap3/Hap4` shares
  the high phenotype class, providing the main 4 -\> 2 refinement
  example;
- `seed_weight`: all four haplotypes retain distinct means, but
  genotype-nearest haplotypes remain phenotypically closer than
  cross-cluster haplotypes;
- `flowering_time`: no GeneA haplotype effect was designed, so the four
  haplotypes can collapse into one phenotype group under the default
  significance criterion.

The workflow follows nine steps:

1.  prepare the original haplotypes and phenotype data;
2.  inspect genotype similarity and phenotype correspondence before
    refinement;
3.  refine GeneA haplotypes using `protein_content`;
4.  inspect the `HapRefined` result object;
5.  inspect pairwise tests and the original-to-refined mapping;
6.  plot the grouped and collapsed refined haplotype structures;
7.  plot phenotype distributions for refined groups;
8.  evaluate trait choice, multiple traits, and refinement thresholds;
9.  interpret refined groups and pass them downstream.

### Step 1. Prepare the original haplotypes and phenotype data

Locate and read the same demo files used in the previous workflows:

\`\`\`{r} gp_file \<- system.file( “extdata”, “gtr_demo.genePredExt”,
package = “GeneTrackR” )

vcf_file \<- system.file( “extdata”, “gtr_demo_variants.vcf”, package =
“GeneTrackR” )

pheno_file \<- system.file( “extdata”, “gtr_demo_pheno.tsv”, package =
“GeneTrackR” )

gp \<- read_genepred( gp_file, format = “genePredExt”, verbose = FALSE,
progress = FALSE )

vcf \<- read_vcf( vcf_file, mode = “memory”, keep_genotype = TRUE,
verbose = FALSE, progress = FALSE )

pheno \<- read_pheno( pheno_file, verbose = FALSE, progress = FALSE )


    Rebuild the standard GeneA gene-body haplotypes exactly as in `08-haplotype.qmd`:

    ```{r}
    hap_gene <- hap_gene_variant(
      vcf,
      annotation = gp,
      gene_id = "GeneA",
      genotype_mode = "string"
    )

    hap_gene
    hap_gene$haplotypes

For the bundled demo, the expected starting point is:

- 11 GeneA variants;
- 36 samples;
- 4 original haplotypes;
- 9 samples per haplotype.

Refinement is therefore performed on a fixed biological definition of
the original haplotypes. Do not silently change the gene region, variant
filter, missing-data rule, or sample set immediately before refinement,
because those changes would redefine the groups being compared.

### Step 2. Inspect genotype similarity and phenotype correspondence before refinement

Refinement is easiest to interpret when the phenotype grouping is
compatible with the underlying genotype geometry. First calculate the
number of allele-state differences between the four original GeneA
haplotypes:

\`\`\`{r} hap_variant_ids \<-
as.character(hap_gene$`variants`$variant_id)

hap_alleles \<- as.matrix( hap_gene\$haplotypes\[, ..hap_variant_ids\] )

hap_ids \<- as.character(hap_gene$`haplotypes`$hap_id)

hap_distance \<- outer( seq_along(hap_ids), seq_along(hap_ids),
Vectorize(function(i, j) { sum(hap_alleles\[i, \] != hap_alleles\[j, \])
}) )

dimnames(hap_distance) \<- list(hap_ids, hap_ids) hap_distance


    The deterministic GeneA design is hierarchical:

    ```text
              Hap1  Hap2  Hap3  Hap4
    Hap1         0     1     9    11
    Hap2         1     0    10    10
    Hap3         9    10     0     2
    Hap4        11    10     2     0

Therefore:

- `Hap1` and `Hap2` are the closest genotype pair and differ at only one
  of the 11 GeneA variants;
- `Hap3` and `Hap4` are the second closest pair and differ at only two
  variants;
- comparisons between these two clusters differ at 9-11 variants.

The phenotype design follows the same structure. Inspect
`protein_content` before refinement:

\`\`\`{r} protein_before \<- plot_hap_pheno( hap_gene, phenotype =
pheno, traits = “protein_content”, min_hap_samples = 3, test_method =
“t.test”, p_adjust = “BH”, p_value_type = “adjusted”, show_points =
TRUE, fill_palette = “Paired” )

protein_before$`figure
protein_before`$summary protein_before\$pvalue


    The expected phenotype means are:

    ```text
    Hap1 ≈ 38
    Hap2 ≈ 38
    Hap3 ≈ 44
    Hap4 ≈ 44

This correspondence is intentional:

``` text
genotype-nearest pair        phenotype class
Hap1 + Hap2                  low protein_content
Hap3 + Hap4                  high protein_content
```

`varA03` follows the same cluster split: `Hap1/Hap2` carry the REF class
and `Hap3/Hap4` carry the ALT class. The main refinement example
therefore groups haplotypes that are similar in both multi-variant
genotype structure and phenotype, rather than joining genetically
distant haplotypes only because their phenotype means happen to match.

### Step 3. Refine haplotypes using `protein_content`

Run phenotype-guided refinement:

\`\`\`{r} refined_protein \<- refine_haplotype( hap_gene, phenotype =
pheno, traits = “protein_content”, min_hap_samples = 3, test_method =
“t.test”, p_adjust = “BH”, alpha = 0.05, effect_threshold = 0.5 )

refined_protein refined_protein\$refined_haplotypes


    The expected demo result is **4 original haplotypes -> 2 refined haplotypes**, with 18 samples in each refined group:

    ```text
    Hap1 + Hap2  -> one refined group
    Hap3 + Hap4  -> one refined group

This is the key truth condition for the example: refinement follows the
two nearest genotype clusters rather than crossing between them.

The two refinement criteria have different roles:

- `alpha = 0.05`: an original haplotype pair can be merged only when its
  adjusted pairwise P value is greater than the cutoff;
- `effect_threshold = 0.5`: the absolute difference between the two
  haplotype phenotype means must also be no larger than 0.5
  protein-content units.

Using an effect threshold is recommended when the biological
interpretation depends on similarity. A non-significant P value alone
means that a difference was not detected; it is not proof that the two
groups are biologically equivalent.

### Step 4. Inspect the `HapRefined` object

[`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md)
returns a `HapRefined` object rather than only a plotting result:

`{r} names(refined_protein)`

The main components are:

- `original_hap`: the original `HapVariant` supplied to
  [`refine_haplotype()`](https://renscq.github.io/GeneTrackR/reference/refine_haplotype.md);
- `refined_hap`: a refined `HapVariant` compatible with the standard
  haplotype plotting/phenotype functions;
- `refined_haplotypes`: one row per refined haplotype group;
- `sample_refined_haplotypes`: sample-level refined group assignments
  while retaining the original haplotype ID;
- `haplotype_map`: mapping from each original haplotype to its refined
  haplotype;
- `trait_group_map`: trait-specific phenotype grouping before multiple
  traits are combined;
- `phenotype_summary`: phenotype summary statistics for the original
  haplotypes;
- `pairwise_test`: pairwise phenotype tests and merge decisions;
- `plot_data`: matched sample-level phenotype data used during
  refinement;
- `parameters`: the refinement settings used to generate the object.

Inspect the original-to-refined mapping directly:

`{r} refined_protein$haplotype_map refined_protein$sample_refined_haplotypes`

`refined_protein$refined_hap` is intentionally a normal
`HapVariant`-compatible object:

`{r} refined_protein$refined_hap refined_protein$refined_hap$haplotypes`

This makes it possible to reuse the standard haplotype plotting and
phenotype-analysis infrastructure without rebuilding variants or
genotypes.

### Step 5. Inspect pairwise tests and merge decisions

Do not treat refinement as a black box. Inspect the pairwise evidence
used to generate the refined groups:

`{r} refined_protein$pairwise_test[, c( "trait", "group1", "group2", "mean_group1", "mean_group2", "abs_mean_diff", "p_value", "p_adj", "can_merge" )]`

For the `protein_content` demo, the two mergeable comparisons are
expected to be `Hap1` versus `Hap2` and `Hap3` versus `Hap4`. Both pairs
have:

- nearly identical means;
- adjusted P values greater than 0.05;
- `abs_mean_diff <= 0.5`;
- `can_merge = TRUE`.

All comparisons between the `Hap1/Hap2` and `Hap3/Hap4` genotype
clusters have a large phenotype difference and should not be merged.

The trait-specific grouping can be inspected separately:

`{r} refined_protein$trait_group_map`

The final original-to-refined relationship is stored in:

`{r} refined_protein$haplotype_map`

This separation is useful when multiple traits are used, because each
trait is grouped independently before GeneTrackR combines the
trait-specific group labels into the final refined haplotypes.

### Step 6. Plot the refined haplotype structure

[`plot_refined_hap_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_refined_hap_variant.md)
provides two complementary views.

#### Keep the original haplotypes inside refined blocks

The default `collapse_refined = FALSE` keeps the original haplotype rows
and arranges them into refined groups:

\`\`\`{r} refined_variant_grouped \<- plot_refined_hap_variant(
refined_protein, annotation = gp, min_hap_samples = 3, collapse_refined
= FALSE, show_reference_row = TRUE, variant_label = “pos”,
gene_pos_x_angle = 90, gene_track_legend_position = “top”, gene_palette
= “Paired”, table_palette = “Paired”, variant_palette = “Paired” )

refined_variant_grouped


    This is usually the most informative refinement plot because it shows **which original genotype patterns were grouped together** while retaining the original variant states.

    #### Collapse each refined group to one row

    Set `collapse_refined = TRUE` when a compact summary is required:

    ```{r}
    refined_variant_collapsed <- plot_refined_hap_variant(
      refined_protein,
      annotation = gp,
      min_hap_samples = 3,
      collapse_refined = TRUE,
      show_reference_row = TRUE,
      variant_label = "pos",
      gene_palette = "Paired",
      table_palette = "Paired",
      variant_palette = "Paired"
    )

    refined_variant_collapsed

When original haplotypes within one refined group carry different
alleles at a variant, the collapsed refined row uses the configured
`mixed_label` (default `"mixed"`). A `mixed` state does **not** mean
that an individual sample is heterozygous. It means that the refined
phenotype group contains more than one original haplotype state at that
variant.

For mechanistic interpretation, keep the grouped view available even
when the collapsed view is used in a summary figure.

### Step 7. Plot phenotype distributions after refinement

Use
[`plot_refined_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_refined_hap_pheno.md)
to visualize the phenotype distributions of the refined groups:

\`\`\`{r} refined_protein_pheno \<- plot_refined_hap_pheno(
refined_protein, phenotype = pheno, traits = “protein_content”,
min_hap_samples = 3, test_method = “t.test”, p_adjust = “BH”,
p_value_type = “adjusted”, show_points = TRUE, fill_palette = “Paired” )

refined_protein_pheno$`figure
refined_protein_pheno`$summary refined_protein_pheno\$pvalue


    The return value is the same `GeneTrackRPhenoPlot` structure used by `plot_hap_pheno()`:

    - `$figure`;
    - `$pvalue`;
    - `$summary`;
    - `$bracket`;
    - `$plot_data`.

    This is useful for exporting the numerical evidence together with the refined-group figure rather than saving only the plot.

    ### Step 8. Evaluate trait choice and refinement thresholds

    Refinement is **trait dependent**. The same original haplotypes can produce very different refined structures when a different phenotype is selected.

    #### A strongly haplotype-associated trait may preserve all groups

    `seed_weight` was designed with four distinct means while preserving the same genotype hierarchy: `Hap1/Hap2` remain closer to one another than to `Hap3/Hap4`, and the same is true within the second pair. The within-pair mean differences are nevertheless larger than `effect_threshold = 0.5`, so refinement should retain all four groups:

    ```{r}
    refined_seed_weight <- refine_haplotype(
      hap_gene,
      phenotype = pheno,
      traits = "seed_weight",
      min_hap_samples = 3,
      test_method = "t.test",
      p_adjust = "BH",
      alpha = 0.05,
      effect_threshold = 0.5
    )

    refined_seed_weight$refined_haplotypes

The expected result is four refined groups, showing that refinement does
not automatically force haplotypes to merge.

#### A negative-control trait can collapse all groups

The demo `flowering_time` trait contains no designed GeneA haplotype
effect:

\`\`\`{r} refined_flowering \<- refine_haplotype( hap_gene, phenotype =
pheno, traits = “flowering_time”, min_hap_samples = 3, test_method =
“t.test”, p_adjust = “BH”, alpha = 0.05, effect_threshold = 0.5 )

refined_flowering\$refined_haplotypes


    For this deterministic demo, all four original haplotypes have the same mean flowering time and can collapse into one refined group.

    This example is also a reminder that **a refined group is conditional on the selected phenotype**. It should not be interpreted as a universal genetic equivalence class.

    #### Multiple traits make the refinement criterion stricter

    Several numeric traits can be supplied together:

    ```{r}
    refined_multi_trait <- refine_haplotype(
      hap_gene,
      phenotype = pheno,
      traits = c("protein_content", "seed_weight"),
      min_hap_samples = 3,
      test_method = "t.test",
      p_adjust = "BH",
      alpha = 0.05,
      effect_threshold = 0.5
    )

    refined_multi_trait$haplotype_map

Each trait is grouped independently. GeneTrackR then combines the
trait-specific group labels, so two original haplotypes remain together
in the final result only when their complete multi-trait grouping
pattern is the same. In this demo, adding the strongly discriminating
`seed_weight` trait prevents the 4 -\> 2 collapse produced by
`protein_content` alone.

For real data, evaluate the sensitivity of the result to:

- `min_hap_samples`;
- `test_method`;
- `p_adjust`;
- `alpha`;
- `effect_threshold`;
- the selected trait set.

Do not choose these parameters only to obtain a desired number of
refined groups.

### Step 9. Interpret and reuse refined haplotypes

Phenotype-guided refinement changes the **group label**, not the
underlying VCF or the original haplotype definition. Keep both levels of
information:

`{r} refined_protein$original_hap$haplotypes refined_protein$haplotype_map refined_protein$refined_haplotypes`

For the `protein_content` demo, the final grouping is expected to be
`Hap1/Hap2` versus `Hap3/Hap4`. The two groups retain different allele
classes at `varA03` and at the broader `p13` cluster, including
`varA01`, `varLD01`-`varLD06`, and `varA04`. Within the low group,
`varA02` is heterogeneous; within the high group, `varA02` and `varA05`
are heterogeneous. These positions appear as `mixed` in the collapsed
view. This is the intended refinement logic: genetically close
haplotypes with the same phenotype class are combined, while the
variants that consistently separate the two genotype/phenotype clusters
remain visible.

The refined object can be passed directly to the refined plotting
wrappers, while the **original** `hap_gene` remains the preferred input
for variant-level effect prioritization:

\`\`\`{r} \# Refined-group visualization \#
plot_refined_hap_variant(refined_protein, annotation = gp) \#
plot_refined_hap_pheno(refined_protein, phenotype = pheno, traits =
“protein_content”)

# Variant-level prioritization keeps the original sample genotypes

# effect \<- plot_variant_effect(

# hap_gene,

# phenotype = pheno,

# traits = “protein_content”,

# min_group_samples = 3

# )


    For reproducible analysis, retain and export at least:

    - `haplotype_map`;
    - `pairwise_test`;
    - `phenotype_summary`;
    - `refined_haplotypes`;
    - `sample_refined_haplotypes`;
    - the refinement `parameters`;
    - the grouped and/or collapsed figure.

    This preserves the evidence needed to explain why each original haplotype was or was not merged.

    ## Variant-effect prioritization

    `plot_variant_effect()` evaluates every natural variant retained in a `HapVariant` object against one or more numeric phenotype traits. It is intended for **regional prioritization**: estimate genotype-associated phenotype effects across all variants, rank the strongest candidates, and then interpret those candidates together with LD, gene annotation, single-variant phenotype plots, and haplotype refinement.

    The GeneA demo is deliberately useful for this workflow. `varA03` is the designed `protein_content` effect variant, but it shares the same `p13` genotype partition with several neighboring variants. Consequently, multiple linked variants have the same large phenotype effect. This demonstrates an important principle: a variant-effect scan can identify an associated variant set, but effect size alone cannot distinguish a causal variant from tightly linked passenger variants.

    This module follows nine steps:

    1. prepare a binary-coded GeneA `HapVariant` and phenotype table;
    2. understand the genotype coding used for effect direction;
    3. calculate the primary `protein_content` variant-effect scan;
    4. inspect the `GeneTrackRVariantEffectPlot` result object;
    5. rank and summarize candidate variants;
    6. interpret tied effects in the context of shared genotype patterns and LD;
    7. inspect signed effects and their interpretation limits;
    8. validate the lead example with a direct single-variant phenotype analysis;
    9. compare multiple traits and connect variant effects to refinement and downstream validation.

    ### Step 1. Prepare GeneA variants and phenotype data

    Locate and read the same demo VCF, annotation, and phenotype files used in the haplotype and refinement workflows:

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

    pheno_file <- system.file(
      "extdata",
      "gtr_demo_pheno.tsv",
      package = "GeneTrackR"
    )

    vcf <- read_vcf(
      vcf_file,
      mode = "memory",
      keep_genotype = TRUE,
      verbose = FALSE,
      progress = FALSE
    )

    gp <- read_genepred(
      gp_file,
      format = "genePredExt",
      verbose = FALSE,
      progress = FALSE
    )

    pheno <- read_pheno(
      pheno_file,
      verbose = FALSE,
      progress = FALSE
    )

For effect screening, construct a second GeneA `HapVariant` with
`genotype_mode = "code"`:

\`\`\`{r} hap_effect \<- hap_gene_variant( vcf, annotation = gp, gene_id
= “GeneA”, genotype_mode = “code” )

hap_effect


    The region is the same canonical GeneA gene body used in the haplotype workflow: 11 variants and 36 samples. Only the genotype representation changes.

    ### Step 2. Understand genotype coding before interpreting effect direction

    With `genotype_mode = "code"`, GeneTrackR uses a binary presence representation:

    - `0`: reference genotype with no ALT allele present;
    - `1`: at least one ALT allele is present.

    This is useful for variant-effect screening because a two-group signed effect becomes interpretable as:

    ```text
    mean(ALT-present group) - mean(reference group)

However, this coding is **not an additive 0/1/2 dosage model**.
Heterozygous and homozygous-ALT genotypes are both represented as `1`.
If heterozygous and homozygous-ALT effects must be separated, inspect
the variant directly with
`plot_variant_pheno(..., genotype_mode = "string")` rather than
interpreting the binary screen as an additive genetic model.

The original string-coded `hap_gene` from the haplotype workflow remains
preferable for visualizing exact allele patterns. The binary-coded
`hap_effect` is created specifically to make the direction of two-group
variant effects easier to interpret.

### Step 3. Calculate the primary protein-content effect scan

The main demo trait is `protein_content`. Run one effect calculation for
every GeneA variant:

\`\`\`{r} effect_protein \<- plot_variant_effect( hap_effect, phenotype
= pheno, traits = “protein_content”, min_group_samples = 3, test_method
= “t.test”, p_adjust = “BH”, effect_type = “absolute”, top_n = 10,
variant_label = “variant_id”, x_axis = “position” )

effect_protein\$figure


    The default screening view uses absolute effect size on the y-axis. Point size also reflects absolute effect size, point transparency reflects `-log10(adjusted P)`, and point color retains the signed effect direction.

    Using `x_axis = "position"` places variants at their genomic coordinates. Use `x_axis = "index"` when evenly spaced variant labels are easier to read than physical distance.

    ### Step 4. Inspect the result object

    `plot_variant_effect()` returns a `GeneTrackRVariantEffectPlot` object rather than a bare ggplot:

    ```{r}
    class(effect_protein)
    names(effect_protein)

The three main components are:

| Component | Content | Typical use |
|----|----|----|
| `figure` | effect-size ggplot | visualization and export |
| `effect` | one row per trait × variant | candidate ranking and statistical interpretation |
| `plot_data` | long sample-level genotype/phenotype table | detailed QC or custom plotting |

Inspect the effect table directly:

`{r} effect_protein$effect[, .( variant_id, chrom, pos, group_n, sample_n, effect, abs_effect, low_group, high_group, low_group_mean, high_group_mean, p_value, p_adj )]`

For a binary `code`-mode variant, `effect` is the mean difference
between genotype group `1` and genotype group `0`. `abs_effect` is its
absolute magnitude. `p_adj` is adjusted within each trait using the
method supplied through `p_adjust`.

### Step 5. Rank candidate variants

Sort variants by absolute phenotype effect:

\`\`\`{r} protein_rank \<- effect_protein\$effect\[ order(-abs_effect,
p_adj), .( variant_id, pos, effect, abs_effect, p_value, p_adj )\]

protein_rank


    The deterministic GeneA demo has the following designed effect structure:

    ```text
    varA01 / varA03 / varLD01-varLD06 / varA04 : absolute effect ~= 6
    varA05                                       : absolute effect ~= 4
    varA02                                       : absolute effect ~= 0

Retrieve the strongest class explicitly:

\`\`\`{r} protein_top \<- effect_protein\$effect\[ abs_effect \>= 5.9,
.(variant_id, pos, effect, abs_effect, p_adj)\]

protein_top\[order(pos)\]


    `varA03` is therefore correctly recovered as a high-effect variant, but it is **not uniquely ranked first**. That is expected and biologically informative rather than an error.

    ### Step 6. Interpret tied effects together with genotype structure and LD

    The variants with an effect of approximately `+6` share the same `p13` genotype split in the deterministic demo: `Hap1/Hap2` are in the reference class and `Hap3/Hap4` are in the ALT-present class. They consequently partition the same 36 samples into the same two phenotype groups and receive the same effect estimate.

    Inspect their variant metadata:

    ```{r}
    retrieve_vcf(
      vcf,
      variant_id = protein_top$variant_id,
      verbose = FALSE
    )[, .(
      variant_id,
      chrom,
      pos,
      ref,
      alt,
      variant_type,
      info
    )]

This is the key interpretation boundary of the effect scan:

``` text
large phenotype effect
        !=
proof of causal variation
```

Variants that are identical or nearly identical in genotype pattern
cannot be separated by association effect alone. Combine the effect
ranking with:

- the LD workflow in `10-ld.qmd`;
- gene position and functional annotation;
- coding/regulatory consequence information;
- independent populations or recombinants that break the linkage;
- experimental validation when causal interpretation is required.

In the demo, `varA03` is marked as the designed `protein_effect`
variant, while `varLD01`-`varLD06` are explicitly retained as linked LD
markers. The equal phenotype effects demonstrate why biological
annotation and recombination information remain necessary after
statistical prioritization.

### Step 7. Inspect signed effects carefully

A signed-effect plot can be generated directly:

\`\`\`{r} effect_protein_signed \<- plot_variant_effect( hap_effect,
phenotype = pheno, traits = “protein_content”, min_group_samples = 3,
test_method = “t.test”, p_adjust = “BH”, effect_type = “signed”, top_n =
10, variant_label = “variant_id”, x_axis = “position” )

effect_protein_signed\$figure


    For the present binary `code`-mode workflow, positive values mean that genotype group `1` has a higher phenotype mean than genotype group `0`. For `varA03`, this corresponds to ALT presence being associated with higher `protein_content`.

    Check the designed site directly:

    ```{r}
    effect_protein_signed$effect[
      variant_id == "varA03",
      .(
        variant_id,
        effect,
        abs_effect,
        low_group,
        high_group,
        low_group_mean,
        high_group_mean,
        p_adj
      )
    ]

The expected effect is approximately `+6`.

Do not generalize this interpretation blindly to every possible
`HapVariant`. With string-coded genotypes, the two-group order can
follow genotype-label ordering rather than a universal REF-to-ALT
direction. With more than two genotype groups,
[`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md)
reports the maximum-minus-minimum mean difference, so the result is
primarily an effect magnitude rather than a simple ALT-minus-REF
contrast. For region-wide prioritization, `effect_type = "absolute"` is
therefore the safer default.

### Step 8. Validate varA03 with a direct single-variant phenotype plot

The regional effect scan tells us that `varA03` belongs to the
high-effect set. Confirm the genotype/phenotype relationship directly
with
[`plot_variant_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_pheno.md):

\`\`\`{r} varA03_pheno \<- plot_variant_pheno( variant = vcf, phenotype
= pheno, variant_id = “varA03”, traits = “protein_content”,
genotype_mode = “code”, min_group_samples = 3, p_adjust = “BH”,
p_value_type = “adjusted” )

varA03_pheno$`figure
varA03_pheno`$summary varA03_pheno\$pvalue


    Retrieve the variant itself to connect the binary genotype classes to its alleles:

    ```{r}
    retrieve_vcf(
      vcf,
      variant_id = "varA03",
      verbose = FALSE
    )[, .(
      variant_id,
      chrom,
      pos,
      ref,
      alt,
      variant_type,
      info
    )]

The demo truth is:

``` text
varA03: C > G
code 0: REF class, Hap1/Hap2, protein_content ~= 38
code 1: ALT-present class, Hap3/Hap4, protein_content ~= 44
ALT-present minus REF effect ~= +6
```

This direct check is preferable to inferring allele direction from the
regional ranking plot alone.

### Step 9. Compare multiple traits and connect downstream evidence

[`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md)
can evaluate several numeric traits in one call:

\`\`\`{r} effect_multi \<- plot_variant_effect( hap_effect, phenotype =
pheno, traits = c( “protein_content”, “seed_weight”, “plant_height”,
“flowering_time” ), min_group_samples = 3, test_method = “t.test”,
p_adjust = “BH”, effect_type = “absolute”, top_n = 5, variant_label =
“variant_id”, x_axis = “position” )

effect_multi\$figure


    Numeric traits may be stored as either integer or double columns. `plot_variant_effect()` normalizes selected numeric traits internally before reshaping them, so mixed numeric storage modes can be analyzed together without `data.table::melt()` coercion warnings.

    Inspect the same variant across traits:

    ```{r}
    effect_multi$effect[
      variant_id == "varA03",
      .(
        trait,
        effect,
        abs_effect,
        p_value,
        p_adj
      )
    ]

For the deterministic demo, `varA03` has approximately:

``` text
protein_content : +6
seed_weight     : +8
plant_height    : +6
flowering_time  :  0
```

The first three effects arise because these traits follow the broader
GeneA genotype hierarchy, whereas `flowering_time` is the designed
negative control. P-value adjustment is performed separately within each
trait.

The complete evidence chain is therefore:

``` text
GeneA HapVariant
    -> region-wide variant-effect ranking
    -> identify high-effect genotype partitions
    -> compare with LD and functional annotation
    -> validate specific variants with plot_variant_pheno()
    -> compare with haplotype refinement
    -> prioritize variants for biological validation
```

For the demo, `protein_content` refinement groups `Hap1/Hap2` and
`Hap3/Hap4`, and `varA03` separates exactly those two refined
genotype/phenotype clusters. The surrounding `p13` variants receive the
same statistical effect because they carry the same genotype partition.
That agreement is useful evidence for the regional signal, but it should
not be mistaken for proof that every tied variant is causal.

## Export tracks, figures, tables, and analysis objects

GeneTrackR produces several different kinds of output, and they should
not all be saved in the same way. Genomic tracks should be written with
the corresponding GeneTrackR writer, figures should be exported from the
returned ggplot/patchwork object, result tables should be written as
tabular files, and complete S3 analysis objects should be retained with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) when the full
calculation state is needed later.

A reproducible export therefore separates four output layers:

1.  genomic track files for use by GeneTrackR or genome browsers;
2.  PDF/PNG figures for reports and manuscripts;
3.  tabular results for inspection and downstream statistics;
4.  complete R objects for exact reloading without recalculation.

This module follows eight steps:

1.  prepare a structured output directory and representative GeneTrackR
    objects;
2.  export annotation, signal, and variant tracks;
3.  save direct and class-contained figures as PDF/PNG;
4.  export haplotype, phenotype, LD, refinement, and variant-effect
    tables;
5.  export matrix-shaped results such as the LD matrix;
6.  optionally collect related tables into an Excel workbook;
7.  save complete GeneTrackR objects with
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html);
8.  record a manifest and session information and verify that saved
    objects reload correctly.

### Step 1. Prepare output directories and representative objects

Create separate directories for track files, figures, tables, and
serialized R objects:

\`\`\`{r} export_root \<- file.path(tempdir(), “GeneTrackR_export”)
export_dirs \<- file.path( export_root, c(“tracks”, “figures”, “tables”,
“objects”) )

invisible(lapply(export_dirs, dir.create, recursive = TRUE, showWarnings
= FALSE))

track_dir \<- export_dirs\[1\] figure_dir \<- export_dirs\[2\] table_dir
\<- export_dirs\[3\] object_dir \<- export_dirs\[4\]


    Read the demo inputs and build representative downstream objects. These are the same object classes used in the preceding workflow modules:

    ```{r}
    gp <- read_genepred(
      system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
      format = "genePredExt",
      verbose = FALSE,
      progress = FALSE
    )

    vcf <- read_vcf(
      system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR"),
      mode = "memory",
      keep_genotype = TRUE,
      verbose = FALSE,
      progress = FALSE
    )

    pheno <- read_pheno(
      system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR"),
      verbose = FALSE,
      progress = FALSE
    )

    hap_gene <- hap_gene_variant(
      vcf,
      annotation = gp,
      gene_id = "GeneA",
      genotype_mode = "string"
    )

    hap_effect <- hap_gene_variant(
      vcf,
      annotation = gp,
      gene_id = "GeneA",
      genotype_mode = "code"
    )

    hap_pheno <- plot_hap_pheno(
      hap_gene,
      phenotype = pheno,
      traits = "seed_weight",
      min_hap_samples = 3,
      p_value_type = "adjusted"
    )

    ld_result <- compute_ld_block(
      vcf,
      chrom = "chr1",
      start = 12342620,
      end = 12355500,
      variant_type = "snp",
      method = "r2",
      verbose = FALSE
    )

    ld_result <- plot_ld_block(
      ld_result,
      annotation = gp,
      show_region = TRUE,
      show_variant_labels = FALSE
    )

    refined <- refine_haplotype(
      hap_gene,
      phenotype = pheno,
      traits = "protein_content",
      min_hap_samples = 3,
      effect_threshold = 0.5
    )

    refined_variant_figure <- plot_refined_hap_variant(
      refined,
      annotation = gp,
      min_hap_samples = 3
    )

    refined_pheno <- plot_refined_hap_pheno(
      refined,
      phenotype = pheno,
      traits = "protein_content",
      min_hap_samples = 3,
      p_value_type = "adjusted"
    )

    effect_result <- plot_variant_effect(
      hap_effect,
      phenotype = pheno,
      traits = "protein_content",
      min_group_samples = 3,
      effect_type = "absolute",
      x_axis = "position"
    )

Keeping these objects separate is useful because their public contracts
differ.
[`plot_hap_variant()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_variant.md)
returns a figure directly, whereas
[`plot_hap_pheno()`](https://renscq.github.io/GeneTrackR/reference/plot_hap_pheno.md)
and
[`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md)
return result objects containing `$figure`;
[`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
returns an updated `LDTrack` with its figure stored in `$figure`.

### Step 2. Export genomic track files

#### Annotation tracks

Use
[`write_feature()`](https://renscq.github.io/GeneTrackR/reference/write_feature.md)
for annotation conversion. Annotation coordinate conventions and
format-specific details are covered in the annotation module; the export
module only centralizes representative calls:

\`\`\`{r} write_feature( gp, file = file.path(track_dir,
“GeneA_demo.genePredExt”), format = “genepredext”, overwrite = TRUE )

write_feature( gp, file = file.path(track_dir, “GeneA_demo.bed12”),
format = “bed12”, overwrite = TRUE )


    #### Signal tracks

    Read an in-memory signal track and export it with `write_bwg()`:

    ```{r}
    rnaseq_plus <- read_bwg(
      system.file(
        "extdata",
        "gtr_demo_rnaseq_plus.bedgraph",
        package = "GeneTrackR"
      ),
      format = "bedgraph",
      sample_names = "RNA_seq_plus",
      strand = "+",
      mode = "memory",
      verbose = FALSE
    )

    bedgraph_files <- write_bwg(
      rnaseq_plus,
      outdir = track_dir,
      format = "bedgraph",
      overwrite = TRUE
    )

    bedgraph_files

bigWig output uses the bundled libBigWig backend and requires chromosome
sizes:

\`\`\`{r} chrom_sizes_file \<- system.file( “extdata”,
“gtr_demo.chrom.sizes”, package = “GeneTrackR” )

bigwig_files \<- write_bwg( rnaseq_plus, outdir = track_dir, format =
“bigwig”, chrom_sizes = chrom_sizes_file, overwrite = TRUE )

bigwig_files


    #### Variant tracks

    Use `write_vcf()` when a standardized site-level VCF is required:

    ```{r}
    write_vcf(
      vcf,
      file = file.path(track_dir, "gtr_demo.sites.vcf"),
      overwrite = TRUE
    )

The current
[`write_vcf()`](https://renscq.github.io/GeneTrackR/reference/write_vcf.md)
writes the standard site columns `CHROM`, `POS`, `ID`, `REF`, `ALT`,
`QUAL`, `FILTER`, and `INFO`. It does **not** guarantee round-trip
preservation of `FORMAT` and sample genotype columns. Keep the original
genotype-containing VCF as the archival source whenever genotypes must
be preserved exactly.

### Step 3. Save figures as PDF and PNG

Use the actual figure object rather than passing an entire GeneTrackR
result class to `ggsave()`.

Define a small local helper when both PDF and PNG are required:

\`\`\`{r} save_gtr_figure \<- function(plot, stem, outdir, width = 8,
height = 6, dpi = 300) { ggplot2::ggsave( filename = file.path(outdir,
paste0(stem, “.pdf”)), plot = plot, width = width, height = height )

ggplot2::ggsave( filename = file.path(outdir, paste0(stem, “.png”)),
plot = plot, width = width, height = height, dpi = dpi ) }


    A direct-return figure can be saved immediately:

    ```{r}
    hap_variant_figure <- plot_hap_variant(
      hap_gene,
      annotation = gp,
      min_hap_samples = 3
    )

    save_gtr_figure(
      hap_variant_figure,
      stem = "GeneA_haplotype_variant",
      outdir = figure_dir,
      width = 10,
      height = 7
    )

For structured result classes, extract `$figure` explicitly:

\`\`\`{r} save_gtr_figure( hap_pheno\$figure, stem =
“GeneA_seed_weight_haplotype”, outdir = figure_dir, width = 7, height =
5 )

save_gtr_figure( ld_result\$figure, stem = “GeneA_LD_gradient”, outdir =
figure_dir, width = 10, height = 8 )

save_gtr_figure( refined_variant_figure, stem =
“GeneA_refined_haplotype_variant”, outdir = figure_dir, width = 10,
height = 7 )

save_gtr_figure( refined_pheno\$figure, stem =
“GeneA_refined_protein_content”, outdir = figure_dir, width = 7, height
= 5 )

save_gtr_figure( effect_result\$figure, stem =
“GeneA_protein_variant_effect”, outdir = figure_dir, width = 9, height =
5 )


    The same rule applies to `plot_variant_pheno()` and `plot_refined_hap_pheno()`: save their `$figure` component. `plot_refined_hap_variant()` returns its figure directly.

    ### Step 4. Export analysis tables

    A `HapVariant` contains several biologically different tables. Export them separately instead of flattening the complete class:

    ```{r}
    data.table::fwrite(
      hap_gene$variants,
      file.path(table_dir, "GeneA.variants.tsv"),
      sep = "\t"
    )

    data.table::fwrite(
      hap_gene$haplotypes,
      file.path(table_dir, "GeneA.haplotypes.tsv"),
      sep = "\t"
    )

    data.table::fwrite(
      hap_gene$sample_haplotypes,
      file.path(table_dir, "GeneA.sample_haplotypes.tsv"),
      sep = "\t"
    )

    data.table::fwrite(
      hap_gene$genotype_long,
      file.path(table_dir, "GeneA.genotype_long.tsv"),
      sep = "\t"
    )

For phenotype association results, retain the statistical tables used to
create the figure:

\`\`\`{r} data.table::fwrite( hap_pheno\$pvalue, file.path(table_dir,
“GeneA.seed_weight.pvalue.tsv”), sep = “ )

data.table::fwrite( hap_pheno\$summary, file.path(table_dir,
“GeneA.seed_weight.summary.tsv”), sep = “ )

data.table::fwrite( hap_pheno\$plot_data, file.path(table_dir,
“GeneA.seed_weight.plot_data.tsv”), sep = “ )


    Export pairwise LD rather than only the heatmap image:

    ```{r}
    data.table::fwrite(
      ld_result$data,
      file.path(table_dir, "GeneA.LD_pairs.tsv"),
      sep = "\t"
    )

    data.table::fwrite(
      ld_result$variants,
      file.path(table_dir, "GeneA.LD_variants.tsv"),
      sep = "\t"
    )

For phenotype-guided refinement, the mapping and decision tables are
essential for explaining why original haplotypes were merged:

\`\`\`{r} data.table::fwrite( refined\$haplotype_map,
file.path(table_dir, “GeneA.refined_haplotype_map.tsv”), sep = “ )

data.table::fwrite( refined\$refined_haplotypes, file.path(table_dir,
“GeneA.refined_haplotypes.tsv”), sep = “ )

data.table::fwrite( refined\$sample_refined_haplotypes,
file.path(table_dir, “GeneA.sample_refined_haplotypes.tsv”), sep = “ )

data.table::fwrite( refined\$pairwise_test, file.path(table_dir,
“GeneA.refinement_pairwise_test.tsv”), sep = “ )

data.table::fwrite( refined\$phenotype_summary, file.path(table_dir,
“GeneA.refinement_pheno_summary.tsv”), sep = “ )

data.table::fwrite( refined_pheno\$pvalue, file.path(table_dir,
“GeneA.refined_protein.pvalue.tsv”), sep = “ )

data.table::fwrite( refined_pheno\$summary, file.path(table_dir,
“GeneA.refined_protein.summary.tsv”), sep = “ )


    For regional variant-effect prioritization, retain both the summarized effects and the sample-level plotting table:

    ```{r}
    data.table::fwrite(
      effect_result$effect,
      file.path(table_dir, "GeneA.protein_variant_effect.tsv"),
      sep = "\t"
    )

    data.table::fwrite(
      effect_result$plot_data,
      file.path(table_dir, "GeneA.protein_variant_effect_plot_data.tsv"),
      sep = "\t"
    )

### Step 5. Export the LD matrix cleanly

A matrix should be converted to a rectangular table with an explicit row
identifier before writing it:

\`\`\`{r} ld_matrix_table \<- data.table::as.data.table(
ld_result\$matrix, keep.rownames = “variant_id” )

data.table::fwrite( ld_matrix_table, file.path(table_dir,
“GeneA.LD_matrix.tsv”), sep = “ )


    This preserves both row and column variant IDs and is easier to reopen in R, Python, spreadsheet software, or downstream matrix tools than a matrix written without row names.

    ### Step 6. Optionally collect tables in an Excel workbook

    For collaborators who prefer one multi-sheet workbook, the same result tables can be collected with the optional `openxlsx` package. This is an interoperability step rather than a GeneTrackR dependency:

    ```{r}
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(
        list(
          variants = as.data.frame(hap_gene$variants),
          haplotypes = as.data.frame(hap_gene$haplotypes),
          sample_haps = as.data.frame(hap_gene$sample_haplotypes),
          hap_pvalue = as.data.frame(hap_pheno$pvalue),
          hap_summary = as.data.frame(hap_pheno$summary),
          LD_pairs = as.data.frame(ld_result$data),
          LD_matrix = as.data.frame(ld_matrix_table),
          refined_map = as.data.frame(refined$haplotype_map),
          refined_test = as.data.frame(refined$pairwise_test),
          variant_effect = as.data.frame(effect_result$effect)
        ),
        file = file.path(table_dir, "GeneA.analysis.xlsx"),
        overwrite = TRUE
      )
    }

TSV remains the dependency-free default and is preferable for very large
tables. Excel is most useful for compact result summaries intended for
manual inspection.

### Step 7. Save complete GeneTrackR objects with RDS

Tabular exports are easy to inspect, but they do not preserve the
complete object state. Save the corresponding S3 objects when an
analysis may need to be reopened without recomputation:

\`\`\`{r} saveRDS( hap_gene, file.path(object_dir,
“GeneA.HapVariant.rds”) )

saveRDS( hap_pheno, file.path(object_dir,
“GeneA.seed_weight.GeneTrackRPhenoPlot.rds”) )

saveRDS( ld_result, file.path(object_dir, “GeneA.LDTrack.rds”) )

saveRDS( refined, file.path(object_dir, “GeneA.HapRefined.rds”) )

saveRDS( effect_result, file.path(object_dir, “GeneA.VariantEffect.rds”)
)


    The RDS file preserves class information, metadata, parameters, intermediate tables, and stored figures. It complements rather than replaces human-readable tables and standard genomic track files.

    ### Step 8. Record provenance and verify the export

    Record the package/session environment together with the result directory:

    ```{r}
    writeLines(
      capture.output(sessionInfo()),
      con = file.path(export_root, "sessionInfo.txt")
    )

    export_manifest <- data.table::data.table(
      file = list.files(export_root, recursive = TRUE)
    )

    data.table::fwrite(
      export_manifest,
      file.path(export_root, "manifest.tsv"),
      sep = "\t"
    )

    export_manifest

Finally, verify that a serialized analysis object can be restored with
its expected class and results:

\`\`\`{r} ld_reloaded \<- readRDS( file.path(object_dir,
“GeneA.LDTrack.rds”) )

class(ld_reloaded) head(ld_reloaded\$data)


    For a complete project, keep the **original inputs**, **standard-format track exports**, **PDF/PNG figures**, **tabular results**, **RDS analysis objects**, and **session information** together. This makes the analysis both easy to inspect and reproducible at the object level.

    ## End-to-end GeneTrackR workflow

    The preceding modules describe individual GeneTrackR tasks in detail. This module condenses them into one recommended end-to-end analysis for a candidate gene/locus. The goal is to keep one consistent set of input objects and pass them forward rather than repeatedly rebuilding slightly different regions or sample sets.

    The demo workflow uses `GeneA` and follows ten steps:

    1. define and read the input files;
    2. validate and summarize the core objects;
    3. inspect the locus with annotation, signal, features, and variants;
    4. retrieve the gene-body variants used for haplotype analysis;
    5. construct and visualize GeneA haplotypes;
    6. test haplotype and single-variant phenotype associations;
    7. inspect the regional LD gradient;
    8. refine genetically similar haplotypes using phenotype evidence;
    9. rank regional variant effects and integrate the evidence;
    10. export the results and adapt the workflow to large real datasets.

    ### Step 1. Define and read the inputs

    Use one explicit set of source files at the start of the analysis:

    ```{r}
    library(GeneTrackR)

    gp_file <- system.file(
      "extdata",
      "gtr_demo.genePredExt",
      package = "GeneTrackR"
    )

    feature_file <- system.file(
      "extdata",
      "gtr_demo_features.bed",
      package = "GeneTrackR"
    )

    vcf_file <- system.file(
      "extdata",
      "gtr_demo_variants.vcf",
      package = "GeneTrackR"
    )

    pheno_file <- system.file(
      "extdata",
      "gtr_demo_pheno.tsv",
      package = "GeneTrackR"
    )

    signal_files <- system.file(
      "extdata",
      c(
        "gtr_demo_rnaseq_plus.bedgraph",
        "gtr_demo_rnaseq_minus.bedgraph",
        "gtr_demo_riboseq_plus.bedgraph",
        "gtr_demo_riboseq_minus.bedgraph"
      ),
      package = "GeneTrackR"
    )

Read the objects once:

\`\`\`{r} gp \<- read_genepred( gp_file, format = “genePredExt”, verbose
= FALSE, progress = FALSE )

features \<- read_bed( feature_file, verbose = FALSE, progress = FALSE )

vcf \<- read_vcf( vcf_file, mode = “memory”, keep_genotype = TRUE,
verbose = FALSE, progress = FALSE )

pheno \<- read_pheno( pheno_file, verbose = FALSE, progress = FALSE )

signal_all \<- read_bwg( signal_files, format = “bedgraph”, sample_names
= c( “RNA_seq_plus”, “RNA_seq_minus”, “Ribo_seq_plus”, “Ribo_seq_minus”
), strand = c(“+”, “-”, “+”, “-”), mode = “memory”, verbose = FALSE )


    For real projects, keep these source objects unchanged and create filtered/retrieved objects for each downstream question. This makes it clear which transformations were applied at each step.

    ### Step 2. Validate and summarize the core objects

    Inspect the annotation and VCF before defining haplotypes or LD regions:

    ```{r}
    validate_feature(gp)
    validate_vcf(vcf)

    summary_feature(gp)
    summary_vcf(vcf)
    summary_pheno(pheno)
    summary_bwg(signal_all)

The most important consistency checks are:

- chromosome names agree among annotation, variants, and signal tracks;
- VCF genotype columns are retained for LD and haplotype analysis;
- phenotype sample IDs overlap the VCF/haplotype sample IDs;
- coordinate conventions have already been normalized by the GeneTrackR
  readers.

The demo phenotype order is intentionally different from the VCF sample
order. Downstream phenotype functions align records by `sample_id`, not
by row position.

### Step 3. Inspect the locus in a browser-like view

Start by viewing the candidate gene in its genomic context.
[`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
can combine the gene model, RNA/Ribo signal, BED-like features, and VCF
variants:

\`\`\`{r} browser_geneA \<- plot_tracks( annotation = gp, signal =
signal_all, features = features, variants = vcf, gene_id = “GeneA”,
samples = c(“RNA_seq_plus”, “Ribo_seq_plus”), strand = “+”, signal_type
= “bar”, feature_color_by = “auto”, feature_max_legend_levels = 5,
direction_mode = “transcript” )

browser_geneA


    All general plotting palettes default to `Paired`; `plot_ld_block()` is the deliberate exception and defaults to the sequential `Reds` palette. For transcript-centered Ribo-seq inspection, supply `transcript_id` to `plot_tracks()` or use `plot_signal_transcript()` so reading-frame information can be displayed explicitly.

    ### Step 4. Retrieve the gene-body variants

    Inspect the exact variants that will define the GeneA haplotypes:

    ```{r}
    geneA_variants <- retrieve_vcf(
      vcf,
      annotation = gp,
      gene_id = "GeneA",
      as = "data.table",
      verbose = FALSE,
      progress = FALSE
    )

    geneA_variants[, .(
      variant_id,
      chrom,
      pos,
      ref,
      alt,
      variant_type
    )]

The canonical GeneA gene body contains 11 demo variants. Do not silently
add flanking variants to the primary haplotype definition: the upstream
demo variant contains missing/heterogeneous genotypes and intentionally
changes the resulting haplotype structure.

### Step 5. Construct and visualize the GeneA haplotypes

Use the gene body and allele-string representation for the main
biological haplotype display:

\`\`\`{r} hap_gene \<- hap_gene_variant( vcf, annotation = gp, gene_id =
“GeneA”, genotype_mode = “string” )

hap_gene\$haplotypes\[, c( “hap_id”, “sample_n”, “samples” ), with =
FALSE\]


    The deterministic demo contains four balanced haplotypes with nine samples each.

    Plot their variant patterns together with the gene structure:

    ```{r}
    hap_gene_figure <- plot_hap_variant(
      hap_gene,
      annotation = gp,
      min_hap_samples = 3,
      variant_label = "pos",
      show_reference_row = TRUE
    )

    hap_gene_figure

Keep the full `HapVariant` object. Its sample assignments and genotype
tables are the input for phenotype association, refinement, and
variant-effect analysis.

### Step 6. Test phenotype associations

Test the original GeneA haplotypes against a trait with a designed
haplotype effect:

\`\`\`{r} hap_seed \<- plot_hap_pheno( hap_gene, phenotype = pheno,
traits = “seed_weight”, min_hap_samples = 3, test_method = “t.test”,
p_adjust = “BH”, p_value_type = “adjusted” )

hap_seed$`figure
hap_seed`$pvalue


    For a specific natural variant, use the VCF genotype directly. The demo `varA03` separates the two major GeneA genotype/phenotype clusters for `protein_content`:

    ```{r}
    varA03_protein <- plot_variant_pheno(
      variant = vcf,
      phenotype = pheno,
      variant_id = "varA03",
      traits = "protein_content",
      genotype_mode = "code",
      min_group_samples = 3,
      p_adjust = "BH",
      p_value_type = "adjusted"
    )

    varA03_protein$figure

These functions return structured result objects. Retain `$pvalue`,
`$summary`, and `$plot_data` together with `$figure` rather than saving
only the visual result.

### Step 7. Inspect the regional LD gradient

The primary demo LD interval contains 12 SNPs. Six form a perfect
high-LD core and the remaining variants progressively decorrelate,
producing high, intermediate, and low `r2` values:

\`\`\`{r} ld_gradient \<- compute_ld_block( vcf, chrom = “chr1”, start =
12342620, end = 12355500, variant_type = “snp”, method = “r2”, verbose =
FALSE )

ld_gradient \<- plot_ld_block( ld_gradient, annotation = gp, show_region
= TRUE, connect_region = TRUE, show_variant_marker = TRUE,
show_variant_labels = FALSE )

ld_gradient\$figure


    Inspect the numerical result as well as the heatmap:

    ```{r}
    ld_gradient$data[, .(
      variant_i,
      variant_j,
      distance_bp,
      r2
    )]

GeneTrackR calculates pairwise LD; it does not automatically declare a
biological block boundary from a hard-coded threshold. Define the
downstream interval using the lead variant, LD pattern, gene context,
and the biological question.

### Step 8. Refine haplotypes with phenotype evidence

The GeneA genotype structure contains two close haplotype pairs:
`Hap1/Hap2` and `Hap3/Hap4`. Their `protein_content` values follow the
same two-cluster structure, making this trait appropriate for
phenotype-guided refinement:

\`\`\`{r} refined_protein \<- refine_haplotype( hap_gene, phenotype =
pheno, traits = “protein_content”, min_hap_samples = 3, test_method =
“t.test”, p_adjust = “BH”, alpha = 0.05, effect_threshold = 0.5 )

refined_protein$`haplotype_map
refined_protein`$refined_haplotypes


    The expected deterministic grouping is:

    ```text
    Hap1 + Hap2 -> refined group 1
    Hap3 + Hap4 -> refined group 2

Inspect the decision evidence rather than reporting only the final group
number:

`{r} refined_protein$pairwise_test refined_protein$phenotype_summary`

A non-significant pairwise test alone is not proof of biological
equivalence. `effect_threshold` adds an explicit phenotype-distance
requirement and should be chosen in the context of the trait scale and
experimental precision.

### Step 9. Rank variant effects and integrate the evidence

For interpretable two-group effect direction, rebuild the same GeneA
region in binary REF/ALT-presence mode:

\`\`\`{r} hap_effect \<- hap_gene_variant( vcf, annotation = gp, gene_id
= “GeneA”, genotype_mode = “code” )

effect_protein \<- plot_variant_effect( hap_effect, phenotype = pheno,
traits = “protein_content”, min_group_samples = 3, test_method =
“t.test”, p_adjust = “BH”, effect_type = “absolute”, top_n = 10,
variant_label = “variant_id”, x_axis = “position” )

effect_protein\$figure


    Rank the numerical results directly:

    ```{r}
    protein_candidates <- effect_protein$effect[
      trait == "protein_content"
    ][order(-abs_effect, p_adj)]

    protein_candidates[, .(
      variant_id,
      effect,
      abs_effect,
      p_adj
    )]

Several GeneA variants share the same large effect because they encode
the same `Hap1/Hap2` versus `Hap3/Hap4` genotype partition. A large
phenotype effect therefore identifies an associated **variant set**, not
necessarily a unique causal variant.

The preferred interpretation combines all preceding evidence:

``` text
variant effect
    + pairwise LD
    + gene/feature annotation
    + single-variant phenotype association
    + haplotype structure
    + phenotype-guided refinement
    -> candidate variants for biological validation
```

### Step 10. Export results and scale the workflow to real datasets

At the end of the analysis, retain the main objects together so the
result can be reopened without reconstructing the workflow:

\`\`\`{r} analysis_bundle \<- list( annotation = gp, features =
features, haplotype = hap_gene, haplotype_pheno = hap_seed,
variant_pheno = varA03_protein, ld = ld_gradient, refined =
refined_protein, variant_effect = effect_protein )

saveRDS( analysis_bundle, file = file.path(tempdir(),
“GeneA.GeneTrackR_workflow.rds”) )


    Use the export module for a complete project directory containing standard genomic tracks, PDF/PNG figures, tabular results, optional Excel workbooks, complete RDS objects, and session information.

    For large real-world datasets, adapt the same workflow rather than changing its logic:

    - read bgzip-compressed indexed VCF files with `read_vcf(mode = "lazy")` and always query a defined genomic interval before LD or haplotype analysis;
    - use lazy bigWig access when genome-wide signal should not be loaded into memory;
    - select only relevant signal samples before plotting;
    - use `bin_size`, `heatmap_bin_size`, or `heatmap_max_bins` for long/high-density signal regions;
    - prefer `hap_gene_variant()` and `hap_region_variant()` for new code; `hap_variant()` remains a compatibility wrapper;
    - remember that direct plotting functions return ggplot/patchwork figures, whereas phenotype/variant-effect results and `LDTrack` store figures inside structured result objects.

    The complete recommended chain is therefore:

    ```text
    annotation + signal + VCF + phenotype
        -> locus inspection
        -> gene/region variant definition
        -> haplotype construction
        -> phenotype association
        -> LD interpretation
        -> phenotype-guided refinement
        -> variant-effect prioritization
        -> biological candidate selection
        -> reproducible export
