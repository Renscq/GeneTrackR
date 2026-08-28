# GeneTrackR demo data model (v0.5.33)

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
- Variants: 56.
- Primary gene: `GeneA` (`+` strand, two transcripts).
- Negative-strand gene: `GeneB` (`-` strand, two transcripts).

## Designed truth

`GeneA` contains four balanced genotype-defined design groups (9 samples each). The genotype geometry is intentionally hierarchical: GeneTrackR Hap1/Hap2 are the nearest pair and Hap3/Hap4 are the second nearest pair, while cross-cluster haplotypes differ at many loci. `protein_content` follows this same two-cluster structure through `varA03`: Hap1/Hap2 carry the low class and Hap3/Hap4 carry the high class. `seed_weight` and `plant_height` retain four distinct means, but the nearest genotype pairs also remain phenotypically closer than cross-cluster pairs. `flowering_time` is a negative-control phenotype with the same within-group value pattern in all four haplotypes.

The primary LD example now contains 12 SNPs (`varLD01`-`varLD12`). `varLD01`-`varLD06` retain the same `p13` genotype pattern and form a perfect high-LD core. `varLD07`-`varLD12` are placed immediately downstream of `GeneA` and progressively decorrelate from that core, producing a deterministic LD gradient with strong, intermediate, and weak pairwise `r2` values. The added variants remain outside the `GeneA` gene body, so the canonical 11-variant GeneA haplotype truth is unchanged. The `GeneT` region (`chr2:16995001-17006000`) still contains exactly `varPair01` and `varPair02`, providing a stable two-variant LD plotting case.

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

Synthetic Ribo-seq P-site-like counts generated from the primary protein-coding `+`-strand transcript of each gene.

- Every bedGraph record is exactly 1 bp wide and stores a positive integer count.
- Zero-count bases are omitted, so the track is sparse rather than a continuous CDS waveform.
- Signal is restricted to CDS positions; UTR, intron, intergenic, and lncRNA regions have no Ribo-seq records.
- Total counts are distributed at approximately 80%/10%/10% across frame 0/frame 1/frame 2. Frame 0 is broadly occupied with variable peak heights, while frame 1 and frame 2 occur at different subsets of codons with lower, irregular counts. The signal remains moderately dense while zero-count codons are retained.
- Initiation and termination frame-0 positions contain pronounced peaks.

### `gtr_demo_riboseq_minus.bedgraph`

Synthetic Ribo-seq P-site-like counts generated from the primary protein-coding `-`-strand transcript of each gene.

- Uses the same sparse integer 1-bp CDS-only logic as the plus track.
- Frame enrichment is calculated in transcript direction and then written back to genomic coordinates.
- The same approximately two-fold initiation/termination behavior is preserved for negative-strand coding genes.

The transcript-specific RNA-seq and Ribo-seq weights are defined in `inst/scripts/demo_model/signal_design.tsv`.

## Generated input files

- `gtr_demo.genePredExt`, `gtr_demo.gtf`, `gtr_demo.gff3`: the same canonical annotation in three formats.
- `gtr_demo_features.bed`: promoters, enhancers, QTL/candidate regions, repeats, and conserved intervals.
- `gtr_demo_variants.vcf`: one shared VCF for variant, haplotype, LD, refinement, and variant-effect workflows.
- `gtr_demo_pheno.tsv`: numeric and categorical phenotypes for all 36 VCF samples.
- `gtr_demo_rnaseq_plus.bedgraph`, `gtr_demo_rnaseq_minus.bedgraph`: strand-specific exon-enriched RNA-seq coverage.
- `gtr_demo_riboseq_plus.bedgraph`, `gtr_demo_riboseq_minus.bedgraph`: strand-specific moderately dense heterogeneous integer P-site-like Ribo-seq counts with approximately 80%/10%/10% frame0/frame1/frame2 totals and approximately two-fold initiation/termination frame-0 counts.
- `gtr_demo.chrom.sizes`: chromosome lengths.

