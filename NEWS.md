# GeneTrackR 0.5.23

- Fix `summary_vcf()` so an in-memory `VariantTrack` can be summarized without supplying a genomic region.
- Fix `retrieve_vcf()` so region, gene, and transcript filters are optional for in-memory/full-file queries; ID, type, and pattern filters can now be used independently as documented.
- Add validation for incomplete direct regions and flank arguments used without a gene/transcript locator.
- Add regression tests for whole-object summaries and non-regional VCF retrieval.

# GeneTrackR 0.5.22

- Reorganized `docs/06-variant.qmd` into a continuous VCF workflow covering reading, validation, summaries, region/gene/transcript retrieval, ID/type/pattern filtering, and plotting.
- Clarified `VariantTrack` versus `data.table` return semantics for `retrieve_vcf()`.
- Added explicit guidance for indexed lazy VCF access and strand-aware gene/transcript retrieval.
- Documented the current site-level scope of `write_vcf()` to avoid implying genotype-preserving round-trip export.

# GeneTrackR 0.5.21

- Changed discrete signal palette assignment so sample, group, and frame levels use RColorBrewer class colors in strict level order instead of interpolating between palette endpoints; this also fixes `plot_signal_region()` color ordering.
- Kept heatmap palettes continuous while separating continuous-gradient generation from discrete signal color assignment.
- Reworked demo Ribo-seq counts so frame 0, frame 1, and frame 2 have visibly heterogeneous heights and different occupied codons while preserving approximately 80%/10%/10% total counts.
- Increased frame-0 codon occupancy while reducing off-frame occupancy, producing a dense but visually varied RPF pattern; initiation and termination frame-0 counts remain approximately two-fold over the internal mean.
- Added regression checks for discrete palette order and within-frame Ribo-seq count variability.

# GeneTrackR 0.5.20

- Fixed signal qualitative palettes so discrete levels use RColorBrewer colors strictly in palette order; `frame0`, `frame1`, and `frame2` now map to the first three `Set1` colors by default.
- Reworked demo Ribo-seq counts to moderate density with approximately 80%/10%/10% frame0/frame1/frame2 total-count proportions.
- Reduced demo Ribo-seq initiation and termination peaks to approximately two times the internal frame-0 mean.

# GeneTrackR 0.5.19

- Rebuilt the demo Ribo-seq tracks as sparse integer P-site-like counts instead of assigning a positive synthetic density to every CDS base.
- Ribo-seq demo counts are now strongly enriched in frame 0, retain only rare frame 1/2 background counts, omit zero-count bases from bedGraph, and keep pronounced frame-0 initiation/termination peaks.
- Corrected all protein-coding demo transcript CDS lengths to multiples of three so frame and stop-boundary examples are biologically coherent.
- Updated the demo generator, validator, regression tests, annotation files, signal documentation, and data-model truth table for the revised Ribo-seq model.
- Clarified that whole-gene plots cannot visually resolve a 3-nt pattern across kilobase-scale genomic spans; transcript frame plots or short windows should be used to inspect periodicity.

# GeneTrackR 0.5.18

- Removed the external `bedGraphToBigWig` compatibility backend from `write_bwg()`.
- Removed the public `bigwig_backend` and `bedGraphToBigWig` arguments; `format = "bigwig"` now always uses the bundled libBigWig writer.
- Removed the UCSC-only temporary chromosome-size-file conversion helper and related documentation examples.
- BigWig export now has one deterministic implementation path: GeneTrackR validation followed by bundled libBigWig writing.

# GeneTrackR 0.5.17

- Fixed a Windows access-violation crash in bundled libBigWig writing by opening output bigWig files in binary update mode (`w+b`) instead of text update mode (`w+`).
- Fixed optional libBigWig callback handling in local read/type-detection paths so `NULL` callbacks are passed through safely instead of dereferenced.
- Removed the immediate C++ reopen check from the writer bridge; output validation is performed after the writer returns through the normal GeneTrackR reader path.
- Added Windows-focused regression coverage for multi-chromosome and strand-specific bigWig round trips.

# GeneTrackR 0.5.16

- Fixed Windows compilation of the bundled libBigWig write bridge with current Rcpp and GCC 14 by removing ambiguous comparisons between `CharacterVector` string proxies and `NA_STRING`.
- Character-vector elements in the C++ writer are now extracted through the R C API (`STRING_ELT`) before missing-value validation and string conversion, avoiding Rcpp proxy overload ambiguity across compiler/Rcpp versions.
- Re-audited the BigWig write path from `write_bwg()` through coordinate validation, chromosome ordering, the Rcpp bridge, and bundled libBigWig calls; no API change was required.
- Added an explicit math-library link (`-lm`) for the vendored libBigWig sources, which use `sqrt`, `pow`, and `ceil`; this prevents later linker failures on Unix-like build systems.
- Added a two-sample plus/minus BigWig round-trip regression test matching the RNA-seq documentation workflow.

# GeneTrackR 0.5.15

- Fixed `write_bwg(format = "bigwig")` so the default path writes bigWig files through the bundled third-party libBigWig library in `src/`, without requiring the UCSC `bedGraphToBigWig` executable.
- Added `bigwig_backend = c("libbigwig", "ucsc")`; the UCSC converter remains available as an explicit compatibility backend.
- Added chromosome-order, coordinate, bounds, finite-value, and non-overlap validation before bigWig export.
- Added a bigWig round-trip regression test through the bundled libBigWig reader/writer backend.

# GeneTrackR 0.5.14

- Standardized all R code fences in QMD documentation from ` ```r ` to executable Quarto/knitr chunks using ` ```{r} `.
- No analysis code or public API changes.

# GeneTrackR 0.5.13

- Removed the experimental native bigWig writer added in 0.5.11, including the Rcpp write bridge and R-level native writer helpers.
- Kept the bundled libBigWig third-party source and existing C++ reader backend unchanged for bigWig reading.
- Restored `write_bwg(format = "bigwig")` to the external UCSC `bedGraphToBigWig` conversion workflow until a separate bwtools package is introduced.
- Removed the temporary roxygen2 `load = "source"` workaround from DESCRIPTION so documentation/build behavior returns to the normal package workflow.
- Removed native-writer-specific regression tests and restored the signal documentation to state the UCSC converter requirement.

# GeneTrackR 0.5.12

- Changed the roxygen2 code-loading strategy from the default `pkgload` mode to `source` mode. This allows `devtools::document()` to regenerate documentation on Windows without compiling the package DLL first.
- Native C/C++ compilation is still validated separately by `devtools::load_all()`, `devtools::install()`, or `devtools::check()`.

# GeneTrackR 0.5.11

- Added native bigWig writing through the bundled libBigWig backend; `write_bwg(format = "bigwig")` no longer requires UCSC `bedGraphToBigWig` by default.
- Kept `bedGraphToBigWig` as an optional legacy backend when an executable path is supplied explicitly.
- Added native bigWig round-trip regression coverage and updated the signal documentation to remove the external-tool requirement.

# GeneTrackR 0.5.10

- Reorganized `docs/05-signal.qmd` into a six-step RNA-seq/Ribo-seq workflow: read bedGraph, write bigWig, gene-level plotting, transcript-level plotting, integrated `plot_tracks()` plotting, and region-level plotting.
- Made the signal QMD module self-contained so it can later be used directly as a MkDocs/Quarto documentation page.
- Added explicit output objects for bigWig export and each plotting step, and documented the external UCSC `bedGraphToBigWig` requirement.

# GeneTrackR 0.5.9

- Fixed `signal_palette_direction` for discrete signal plots when only one sample/group is visible. Palette direction is now applied to the palette gradient before colors are sampled, so `signal_palette_direction = -1` also changes a one-track `Blues` plot.
- Audited the shared signal color path used by `plot_signal_transcript()`, `plot_signal_gene()`, `plot_signal_region()`, and `plot_tracks()` for bar, line, heatmap, sample, and group coloring.
- Explicit `signal_colors` are now treated as exact user mappings and are no longer reversed by `signal_palette_direction`.
- Added `signal_track_height` and `gene_track_height` to `plot_signal_transcript()`, `plot_signal_gene()`, and `plot_signal_region()`; defaults remain 3 and 1.
- Updated signal documentation and regression tests for palette direction and panel-height control.

# GeneTrackR 0.5.8

- Converted the README source to Quarto (`README.qmd`) while retaining `README.md` as the synchronized GitHub-facing output.
- Split the long README into portable QMD modules under `docs/` for future MkDocs/Quarto navigation.
- Removed fixed numeric prefixes from documentation headings so navigation order can be managed externally.
- Kept R package vignettes in R Markdown for stable package-vignette compatibility.

# GeneTrackR 0.5.7

- Reorganized the first README workflow around annotation formats and operations.
- Added format-specific read, subset, merge, and write examples for GenePred, GenePredExt, GTF, GFF3, and BED.
- Documented annotation coordinate conventions, BED interval semantics, retrieval modes, merge conflict strategies, cross-format conversion, and annotation object contracts.
- Removed duplicated annotation-export examples from the later export section.

# GeneTrackR 0.5.6

- Reworked README and vignette workflows against the current public API and return-object contracts.
- Fixed obsolete `plot_gene()` example arguments and standardized the documentation object flow across annotation, signal, variant, haplotype, phenotype, LD, refinement, and variant-effect modules.
- Clarified which plotting functions return figures directly, which return result objects with `$figure`, and that `plot_ld_block()` returns an updated `LDTrack` by default.
- Stabilized the recommended GeneA haplotype workflow by keeping the designed upstream missing/heterozygous variant out of the primary four-haplotype example.
- Added documentation workflow regression tests for core object classes and return fields.

# GeneTrackR 0.5.5

- Migrated annotation, variant, haplotype, phenotype, LD, haplotype-refinement, and variant-effect examples to the unified `gtr_demo_*` dataset.
- Added runnable demo examples for GTF, GFF3, BED, VCF, `compute_ld_block()`, and `plot_ld_block()`.
- Added a core-workflows vignette covering annotation, variants, haplotypes, phenotypes, LD, refinement, and variant effects.
- Updated README examples to use valid designed demo variants, coordinates, phenotype traits, and current plotting parameters.
- Migrated regression tests from legacy `example_*` inputs to the deterministic demo model and removed legacy example files.

# GeneTrackR 0.5.4

- Migrated README signal examples to the four strand-specific demo tracks: RNA-seq plus/minus and Ribo-seq plus/minus.
- Added the first signal-track vignette using the unified `gtr_demo_*` dataset.
- Updated signal-related roxygen and `.Rd` examples to use the new strand-specific demo tracks and matching demo genome coordinates.
- Migrated signal plotting tests from legacy `example_signal_*` files to the unified demo data model.

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
