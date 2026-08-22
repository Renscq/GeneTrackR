# GeneTrackR 0.5.3

- Split the demo RNA-seq bedGraph into strand-specific `gtr_demo_rnaseq_plus.bedgraph` and `gtr_demo_rnaseq_minus.bedgraph` files.
- Split the demo Ribo-seq bedGraph into strand-specific `gtr_demo_riboseq_plus.bedgraph` and `gtr_demo_riboseq_minus.bedgraph` files.
- Updated the demo-data generator, validator, and tests to use strand-specific signal tracks.

# GeneTrackR 0.5.2

- Replaced the four generic demo control/treatment strand bedGraph tracks with two assay-specific signal examples: `gtr_demo_rnaseq.bedgraph` and `gtr_demo_riboseq.bedgraph`.
- Added `inst/scripts/demo_model/signal_design.tsv` as the canonical transcript-level source for demo RNA-seq and Ribo-seq signal amplitudes.
- RNA-seq demo coverage is restricted to transcript exons, including UTRs and lncRNA exons; designed intronic and intergenic regions have zero coverage.
- Ribo-seq demo density is restricted to protein-coding CDS positions, uses one-base bedGraph intervals, follows transcript-oriented three-nucleotide periodicity, and contains pronounced initiation and termination peaks.
- Updated demo-data generation, validation, regression tests, and data-model documentation for the new signal design.

# GeneTrackR 0.5.1

- Started the 0.5.x example-data redesign with one deterministic demo genome shared across annotation, signal, variant, haplotype, phenotype, LD, refinement, variant-effect, and export workflows.
- Added canonical model tables under `inst/scripts/demo_model/` for 2 chromosomes, 20 genes, 24 transcripts, 36 samples, 50 variants, interval features, and explicit designed truth.
- Added `generate_demo_data.R` and `validate_demo_data.R` so user-facing `gtr_demo_*` files are generated from one source model instead of being maintained independently.
- Added equivalent GenePredExt, GTF, and GFF3 annotations, a shared VCF, phenotype table, four strand-labelled bedGraph tracks, BED features, and chromosome sizes.
- Designed GeneA with four balanced haplotypes, a perfect six-variant LD block, phenotype-associated and negative-control traits, and GeneT with exactly two variants for stable LD plotting regression.
- Retained legacy `example_*` files temporarily so existing documentation and tests remain functional during the staged 0.5.x migration.

## GeneTrackR 0.4.10

- Fixed LD block region track x-axis labels to rotate and align with genomic ticks.

# GeneTrackR 0.4.8

- Updated the exactly-two-variant case in `plot_ld_block()` so the single LD cell is drawn as a 45-degree diamond rather than a horizontal rectangle.
- The two-variant cell now uses the same rotated-square geometry as cells in the standard triangular LD matrix, with a matching diamond-shaped outer frame.

# GeneTrackR 0.4.7

- Fixed `plot_ld_block()` for exactly two variants. The single pairwise LD value is now drawn as one complete rectangular heatmap cell with a rectangular frame instead of a clipped/empty triangular geometry.
- LD heatmaps with three or more variants retain the existing inverted triangular layout.

# GeneTrackR 0.4.6

- Added strand-direction arrows to the compact gene model drawn by `plot_hap_variant()`. Positive-strand transcripts point right and negative-strand transcripts point left.
- Added `direction_mode` to `plot_hap_variant()` with the same `transcript`, `gene`, `end`, and `none` semantics used by `plot_region()`.

# GeneTrackR 0.4.5

- Updated `plot_hap_pheno()` so multi-trait figures use facet strips without a redundant overall title.
- Added `facet_ncol` (default `3`) so multiple phenotype panels wrap across rows instead of relying on automatic layout.
- Preserved the exact user-specified `traits` order in phenotype facets instead of allowing facet labels to be reordered.

# GeneTrackR 0.4.4

- Rotated `plot_variant_effect()` x-axis tick labels by 90 degrees and aligned them with their ticks to keep large genomic coordinates readable.

# GeneTrackR 0.4.3

- Updated `plot_variant_effect()` so scatter-point color encodes the signed phenotype effect: negative effects are blue, values near zero are light, and positive effects are red. Point size and transparency continue to encode absolute effect size and adjusted-P significance, respectively.

# GeneTrackR 0.4.2

- Standardized all package R script development headers to `dev001`; package versioning is now controlled only by `DESCRIPTION`.
- Reset `GeneTrackR.Rproj` to the standard RStudio project format version (`1.0`) instead of mirroring the package version.
- Removed zero-reference internal helper functions that were no longer used by package code.
- Removed redundant roxygen/NAMESPACE imports for functions already called with explicit package namespaces.
- Removed duplicated package-level global variable declarations that were no longer referenced.
- Removed the empty `src/Makevars.cpp` translation unit and stale local compilation artifacts.

# GeneTrackR 0.4.1

## Bug fixes

- Preserve interval widths when exporting in-memory signal tracks to WIG by writing the correct `span` for each variableStep block.
- Preserve `.gz`/`.bgz` compression suffixes when lazily copying bedGraph or WIG source files, and fail explicitly when the source file is missing or copying fails.
- Keep custom missing-genotype display labels separate from genotype missingness so `non_missing_variant_n` and `min_variant_number` filtering remain correct.
- Exclude missing genotypes from phenotype genotype groups even when a non-NA display label is requested.
- Declare `Rsamtools` as an optional dependency used by indexed VCF and bedGraph queries.
