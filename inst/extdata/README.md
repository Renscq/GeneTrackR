# GeneTrackR demo data model (v0.5.3)

This directory contains one deterministic demo genome shared across GeneTrackR examples.
All `gtr_demo_*` files are generated from the canonical model tables in `inst/scripts/demo_model/`.

## Coordinate conventions

- Canonical model tables: 1-based closed.
- GTF/GFF3/VCF: 1-based coordinates.
- GenePredExt/BED/bedGraph: written using their standard 0-based conventions and converted by GeneTrackR readers.

## Core dimensions

- Chromosomes: 2 (`chr1`, `chr2`).
- Genes: 20.
- Transcripts: 24.
- Samples: 36 (`S01`-`S36`).
- Variants: 50.
- Primary gene: `GeneA` (`+` strand, two transcripts).
- Negative-strand gene: `GeneB` (`-` strand, two transcripts).

## Designed truth

`GeneA` contains four balanced genotype-defined design groups (9 samples each). `seed_weight` and `plant_height` have designed haplotype effects. `protein_content` is designed around `varA03`, whose ALT state is present in DesignHap2/DesignHap3. `flowering_time` is a negative-control phenotype with the same within-group value pattern in all four haplotypes.

`varLD01`-`varLD06` share the same genotype pattern and therefore form a perfect high-LD block (`r2 = 1` in complete samples). The `GeneT` region (`chr2:16995001-17006000`) contains exactly `varPair01` and `varPair02`, providing a stable two-variant LD plotting case.

The phenotype rows are deliberately stored in an order different from the VCF sample columns. Analyses must align samples by `sample_id`, not by row position.

## Signal tracks

Four strand-specific signal files are provided: RNA-seq plus/minus and Ribo-seq plus/minus. Their purpose is to represent two biologically different sequencing assays while also preserving transcript strand information in the demo tracks.

### `gtr_demo_rnaseq_plus.bedgraph`

Synthetic RNA-seq coverage generated from exons of `+`-strand transcripts.

- Exonic regions, including UTRs, contain positive coverage.
- Introns and intergenic regions have no bedGraph records and therefore represent zero coverage.
- Protein-coding and lncRNA transcripts can both have RNA-seq coverage.
- Coverage values from overlapping `+`-strand isoforms are summed.

### `gtr_demo_rnaseq_minus.bedgraph`

Synthetic RNA-seq coverage generated from exons of `-`-strand transcripts.

- Uses the same exon-enriched logic as the plus track.
- Makes strand-aware plotting examples more realistic when plus and minus genes are near each other.

### `gtr_demo_riboseq_plus.bedgraph`

Synthetic Ribo-seq P-site-like density generated from the primary protein-coding `+`-strand transcript of each gene.

- Every bedGraph record is exactly 1 bp wide.
- Signal is restricted to CDS positions; UTR, intron, intergenic, and lncRNA regions have no Ribo-seq records.
- Internal CDS density has a strong three-nucleotide periodic pattern: phase 0 > phase 1 > phase 2.
- Initiation and termination regions contain pronounced peaks.

### `gtr_demo_riboseq_minus.bedgraph`

Synthetic Ribo-seq P-site-like density generated from the primary protein-coding `-`-strand transcript of each gene.

- Uses the same 1-bp CDS-only logic as the plus track.
- Three-nucleotide periodicity is calculated in transcript direction and then written back to genomic coordinates.
- Initiation and termination peaks are preserved for negative-strand coding genes.

The transcript-specific RNA-seq and Ribo-seq weights are defined in `inst/scripts/demo_model/signal_design.tsv`.

## Generated input files

- `gtr_demo.genePredExt`, `gtr_demo.gtf`, `gtr_demo.gff3`: the same canonical annotation in three formats.
- `gtr_demo_features.bed`: promoters, enhancers, QTL/candidate regions, repeats, and conserved intervals.
- `gtr_demo_variants.vcf`: one shared VCF for variant, haplotype, LD, refinement, and variant-effect workflows.
- `gtr_demo_pheno.tsv`: numeric and categorical phenotypes for all 36 VCF samples.
- `gtr_demo_rnaseq_plus.bedgraph`, `gtr_demo_rnaseq_minus.bedgraph`: strand-specific exon-enriched RNA-seq coverage.
- `gtr_demo_riboseq_plus.bedgraph`, `gtr_demo_riboseq_minus.bedgraph`: strand-specific single-base CDS Ribo-seq density with 3-nt periodicity and start/stop peaks.
- `gtr_demo.chrom.sizes`: chromosome lengths.

## Legacy files

Files beginning with `example_` are retained temporarily during the 0.5.x migration so existing examples and tests do not break. They will be removed after all documentation and tests use the new `gtr_demo_*` data.
