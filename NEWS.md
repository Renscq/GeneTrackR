# GeneTrackR 0.8.1

- Expanded Bioconductor interoperability so `as_granges()` now accepts annotation, in-memory signal (`BwgTrack`), and in-memory variant (`VariantTrack`) objects while preserving the existing annotation-level interface.
- Added `retrieve_vcf(..., as = "GRanges")` so indexed or in-memory variant retrieval can feed directly into standard Bioconductor genomic interval workflows.
- Added `S4Vectors`, `grDevices`, and `tools` to explicit runtime imports where GeneTrackR calls their namespaces directly; retained `Rsamtools` in Suggests because tabix-backed querying remains an optional acceleration path.
- Added dependency-placement and cross-track `GRanges` regression coverage, including lazy-track conversion diagnostics.
- Documented the boundary between GeneTrackR's native pure-R bedGraph/WIG/BigWig I/O and Bioconductor's `rtracklayer` infrastructure while keeping the native backend and public signal APIs unchanged.
- Added `tools/bioc_preflight.002.R` with dependency duplication, namespace-declaration, and interoperability checks while retaining the 0.8.0 preflight script for development history.

# GeneTrackR 0.8.0

- Started the Bioconductor-readiness phase while keeping the 0.7.8 pure-R signal implementation and public analysis APIs unchanged.
- Added Bioconductor-oriented DESCRIPTION metadata, including software `biocViews`, a broader package title, and a three-sentence functional description that reflects annotation, signal, variant, haplotype, phenotype, LD, and variant-effect workflows.
- Kept `NeedsCompilation: no`, removed stale generated DESCRIPTION fields and the obsolete roxygen configuration pin, and added `openxlsx` to Suggests for the optional workbook workflow used in an executable vignette.
- Added RiboParser-style README navigation links for pkgdown documentation, issues, citation, and license, and added explicit Citation and License sections.
- Converted all pkgdown article Rmd files into formal R package vignettes with `VignetteIndexEntry`, `VignetteEngine`, and UTF-8 metadata so they can participate in `R CMD build` and `R CMD check`.
- Added a Bioconductor-readiness regression test for metadata and GRanges interoperability and a versioned developer preflight script for repository, vignette, compiled-code, and submission-boundary checks.
- Kept the package at version 0.8.0 during preflight development; the formal Bioconductor submission candidate will use the required 0.99.0 version series.

# GeneTrackR 0.7.8

- Finalized the 0.7.x pure-R signal migration documentation without changing public APIs or signal calculations.
- Removed obsolete Rcpp, bundled libBigWig, `src/`, and compiled-backend wording from installation guides, signal/export tutorials, pkgdown vignettes, generated README content, and public signal help pages.
- Restored `README.qmd` as the canonical include-based source used by `tools/render_readme.R`; the clean 0.7.7 deployment artifact had retained the rendered README but stripped its include directives.
- Updated package-level documentation to state that bedGraph, WIG, and BigWig reading, regional querying, and writing use the built-in pure-R signal I/O architecture.
- Synchronized the public signal Rd pages with their roxygen sources and renamed legacy writer-test labels to match the pure-R implementation.
- Retained schema-v2 `BwgTrack`, 1-based closed internal coordinates, lazy/memory semantics, tabix bedGraph behavior, BigWig chromosome-size requirements, and all 0.7.7 signal behavior unchanged.

# GeneTrackR 0.7.7

- Added full signal-subsystem regression coverage for lazy-versus-memory BigWig retrieval, chromosome-edge intervals, strand semantics, normalization parity, and cross-format round trips.
- Fixed the pure-R package boundary regression test so data.table key attributes are not mistaken for signal-value differences, while keeping compiled-helper absence as a hard requirement.
- Improved the compiled-helper failure diagnostic to identify stale namespace bindings and remind developers to remove `R/bw_cpp_backend.R` and restart R before testing.
- Prevented `summary_bwg()` and `bin_bwg()` from mutating caller-owned data.tables through data.table reference semantics.
- Fixed `merge_bwg(sample_conflict = "keep_first")` so duplicate sample IDs retain records and sequence metadata only from the first source, including lazy objects without in-memory signal data.
- Retained the pure-R bedGraph, WIG, and BigWig I/O architecture, schema-v2 `BwgTrack` contract, and public signal APIs unchanged.

# GeneTrackR 0.7.6

- Removed the compiled libBigWig/Rcpp compatibility backend after native-R BigWig reading, querying, and writing had been migrated and unified in 0.7.2-0.7.5.
- Removed `R/bw_cpp_backend.R` and the complete `src/` tree, including libBigWig C sources, the Rcpp bridge, platform Makevars files, and compiled build artifacts.
- Removed `Rcpp` from `Imports` and `LinkingTo`, removed the generated Rcpp import and GeneTrackR `useDynLib()` registration, and changed package metadata to `NeedsCompilation: no`.
- Replaced the obsolete compiled-reader parity test with native low-level-versus-adapter parity and added package-level regression checks that compiled BigWig helpers are absent from the namespace.
- Updated the BigWig writer architecture test to follow the 0.7.5 dispatcher design (`write_bwg()` -> `write_signal_file_memory()` -> `write_bigwig_native()`) instead of requiring a direct writer call from `write_bwg()`.
- Retained the existing public signal APIs, schema-v2 `BwgTrack` contract, 1-based closed coordinates, tabix bedGraph behavior, and native-R BigWig round-trip behavior.

# GeneTrackR 0.7.5

- Added a unified native-R signal I/O layer for bedGraph, WIG, and BigWig formats while preserving the existing public `read_bwg()`, `retrieve_bwg()`, and `write_bwg()` APIs.
- Centralized in-memory format dispatch through `read_signal_file_memory()` and `write_signal_file_memory()`, and centralized lazy per-sample retrieval through `query_signal_file_region()`.
- Replaced duplicate bedGraph parsing in full-file and tabix query paths with one canonical 0-based half-open to 1-based closed parser.
- Reworked WIG text reading around the stabilized bwTools parser pattern, including fixedStep, variableStep, bedGraph-style records, comments, and gzip/bgzip text connections.
- Centralized signal-region clipping and sequence-metadata construction so all supported formats follow the same canonical six-column signal schema and schema-v2 `seqinfo` contract.
- Retained compatibility wrappers for existing internal text I/O helper names while removing their duplicate implementations from `read_bwg.R` and `write_bwg.R`.
- Added regression coverage for unified dispatch, compressed text inputs, interval clipping, and bedGraph/WIG round trips without coordinate drift.
- Kept Rcpp and `src/` temporarily for compiled-backend parity code only; their final removal remains isolated to the next migration stage.

# GeneTrackR 0.7.4

- Fixed schema-v2 `seqinfo` chromosome subsetting in `retrieve_bwg(..., as = "BwgTrack")` and `slice_bwg(..., as = "BwgTrack")`; the previous helper allowed the `data.table` `chrom` column to mask the function argument and therefore retained unrelated chromosomes.
- Added the native-R BigWig writer adapted from the stabilized bwTools 0.8.10 implementation, including binary packing, chromosome B+ tree construction, R-tree indexing, compressed data blocks, total summaries, and zoom-level generation.
- Switched `write_bwg(..., format = "bigwig")` from the bundled libBigWig/Rcpp writer to the native R writer without changing the public `write_bwg()` argument contract.
- Kept `chrom_sizes` mandatory for GeneTrackR BigWig export in this release to avoid mixing backend migration with API changes.
- Added native writer regression coverage for round-trip signal values, multiple samples, overlap and chromosome-boundary validation, and multi-level chromosome indexes.
- Retained the compiled backend only for temporary parity tests and later cleanup; public BigWig reading and writing now both use native R implementations.

# GeneTrackR 0.7.3

- Switched the public BigWig read path from the bundled libBigWig/Rcpp reader to the native R reader introduced in 0.7.2.
- Added internal `query_bigwig_native()` and `read_bigwig_whole_native()` adapters that return the standard GeneTrackR signal schema.
- Updated `read_bwg()` memory loading and sequence metadata collection, lazy `retrieve_bwg()` queries, and legacy `seqinfo_bwg()` fallback to use the native R backend.
- Updated `BwgTrack$meta$backend` to `GeneTrackR-native-R` for newly read signal objects.
- Retained the compiled reader temporarily only for regression comparison; BigWig writing remains on the compiled backend until the writer migration.
- Added regression coverage that verifies public BigWig read functions no longer reference the compiled reader while preserving interval values, clipping, and schema behavior.

# GeneTrackR 0.7.2

- Added the first native-R BigWig reader layer adapted from the stabilized bwTools 0.8.10 implementation.
- Added internal binary decoding for BigWig headers, chromosome B+ trees, R-tree indexes, compressed data blocks, and 32-bit floating-point signal values.
- Added internal `bigwig_metadata_native()`, `bigwig_seqinfo_native()`, and `bigwig_query_native()` functions using 1-based closed coordinates without Rcpp or `.Call()`.
- Added metadata caching keyed by normalized path, file size, and modification time while keeping the native reader local-file only.
- Added a deterministic BigWig fixture and regression tests for chromosome metadata, interval decoding, clipping, empty ranges, invalid files, and parity with the existing compiled backend.
- Kept `read_bwg()` and `retrieve_bwg()` on the existing libBigWig/Rcpp backend in this release so the new reader can be verified before the public backend switch in 0.7.3.

# GeneTrackR 0.7.1

- Established the schema-v2 `BwgTrack` object foundation for the 0.7.x pure-R signal subsystem migration.
- Added a stable `seqinfo` slot alongside `samples`, `data`, `meta`, and `validation`; chromosome lengths can be known or `NA` when the source format does not encode sequence lengths.
- Standardized `BwgTrack` metadata to `coordinate = "1-based closed"` and `schema_version = "2"` while preserving the existing public constructor argument order.
- Added sequence metadata collection for bigWig inputs and in-memory bedGraph/WIG inputs, and preserved relevant `seqinfo` through BwgTrack retrieval, slicing, and merging.
- Updated `seqinfo_bwg()` to use stored schema-v2 sequence metadata first while retaining the existing libBigWig fallback for legacy objects.
- Added regression tests for schema construction, sequence metadata validation, text-track sequence metadata, subset propagation, and merge propagation.
- Kept the existing libBigWig/Rcpp backend unchanged in this release; native BigWig replacement begins in the following 0.7.x development stages.

# GeneTrackR 0.6.17

- Optimized pkgdown article figure layouts by increasing the graphics canvas according to plot complexity while keeping rendered figures at full article width.
- Standardized vignette web graphics to 150 dpi with `fig.retina = 1`; publication-oriented export examples continue to use 500 dpi.
- Reduced documentation-only plotting text sizes to 9-10 pt for functions that expose `text_size` or `font`, without changing GeneTrackR plotting defaults or package APIs.
- Enlarged annotation-specific figure heights and retained all R implementations, staged `docs/*.qmd` sources, tests, and statistical behavior unchanged.

# GeneTrackR 0.6.16

- Removed the `Documentation QA` GitHub Actions workflow.
- Retained the pkgdown workflow as the single documentation build and deployment workflow.
- Kept README rendering, pkgdown articles, package code, tests, and generated documentation unchanged.

# GeneTrackR 0.6.15

* Fixed pkgdown homepage code rendering by generating `README.md` from `README.qmd` and staged `docs/*.qmd` sources as display-only GitHub-Flavored Markdown; Quarto/knitr `\`\`\`{r}` fences are converted to standard `\`\`\`r` fences instead of being exposed as plain text.
* Added `tools/render_readme.R` as the single deterministic README assembly entry point. The script expands Quarto include directives without executing the 202 documentation code chunks and supports `--check` for CI synchronization checks.
* Added README synchronization to Documentation QA and made the pkgdown workflow regenerate `README.md` before site construction, preventing stale homepage content.
* No R implementation, API, vignette, statistical method, or plotting behavior was changed.

# GeneTrackR 0.6.14

* Removed the noisy GitHub Pages 404 annotation from the pkgdown workflow by probing repository Pages availability before invoking `actions/configure-pages`.
* `actions/configure-pages` now runs only when the Pages REST endpoint is available; repositories without Pages enabled receive a clean pkgdown preview artifact and Step Summary without an Error or Warning annotation.
* Kept pkgdown site construction as a hard CI check and retained automatic native Pages deployment once repository-level Pages is enabled.
* No R implementation, API, vignette, test, or pkgdown content was changed.

# GeneTrackR 0.6.13

* Synchronized the package-level Rd file with the roxygen2 8.1.0 output generated from `DESCRIPTION` URL and BugReports metadata.
* Added the automatically generated Useful links `\seealso{}` section so the documentation QA roxygen drift check remains clean.
* No R implementation, API, vignette, test, pkgdown workflow, or analysis behavior was changed.

# GeneTrackR 0.6.12

* Made the pkgdown GitHub Pages deployment tolerant of repositories where Pages has not yet been enabled: pkgdown site construction remains a hard check, while Pages configuration/deployment is skipped with a clear warning instead of failing the entire workflow.
* Added a fallback `pkgdown-site-preview` workflow artifact when GitHub Pages is unavailable, so the generated site can still be inspected from the Actions run.
* Retained the native GitHub Pages artifact/deploy path automatically when `actions/configure-pages` succeeds; no R implementation, API, vignette, or pkgdown content was changed.

# GeneTrackR 0.6.11

* Re-synchronized roxygen2 8.1.0 generated metadata, `NAMESPACE`, and package-level Rd documentation after the 0.6.10 GitHub commit omitted part of the generated-file update.
* Retained the documentation QA roxygen drift guard and the explicit `roxygen2@8.1.0` CI pin.

# GeneTrackR 0.6.10

- Synchronized roxygen-generated package metadata, `NAMESPACE`, and package-level Rd documentation with roxygen2 8.1.0, matching the generator used by the documentation QA workflow.
- Pinned roxygen2 8.1.0 in the documentation-QA dependency set so future CI runs do not fail solely because a newer roxygen2 release reformats generated documentation.
- Preserved all R implementations, tests, pkgdown articles, staged Quarto sources, and package APIs unchanged; this release is limited to generated-documentation synchronization and CI reproducibility.

# GeneTrackR 0.6.9

- Added a dedicated GitHub Actions documentation-QA workflow that rebuilds the pkgdown site in a clean R process, verifies roxygen-generated documentation is synchronized, checks internal site links, inventories generated figure dimensions, and runs `R CMD check --as-cran`.
- Added a Bioconductor release preflight for main-branch and manually triggered runs using the official `bioconductor/bioconductor_docker:RELEASE_3_23` image and `BiocCheck`; package-specific documentation/coding errors fail CI, while submission-only `biocViews` checks are deferred until the formal Bioconductor submission phase.
- Modernized pkgdown deployment to the native GitHub Pages artifact workflow while retaining `pkgdown-site/` as the generated-site directory.
- Kept the 0.6.8 article, Reference, API, and analysis implementations unchanged; this release is limited to documentation QA and CI integration.

# GeneTrackR 0.6.8

- Added a dedicated export and reproducibility article covering standard genomic track writers, PDF/PNG figure export, tabular result preservation, LD matrix export, optional `openxlsx` workbooks, complete RDS objects, manifests, and `sessionInfo()` provenance.
- Added a production-oriented complete workflow article connecting annotation, interval features, RNA-seq/Ribo-seq signals, variants, haplotypes, phenotype association, LD, phenotype-guided refinement, variant-effect prioritization, and reproducible export through one consistent candidate-locus analysis.
- Added a **Workflow & reproducibility** section to the pkgdown Articles index while retaining `core-workflows` as the shorter Getting Started workflow.
- Kept the staged `docs/13-export.qmd` and `docs/14-workflow.qmd` sources unchanged while migrating their current workflows into formal pkgdown vignettes.

# GeneTrackR 0.6.7

- Added a dedicated linkage-disequilibrium article covering genotype-dosage LD calculation, `LDTrack` result inspection, `r2` and dosage-based D-prime interpretation, triangular heatmaps, regional gene/variant context, population subsets, variant filtering, and downstream regional haplotype handoff.
- Added a dedicated phenotype-guided haplotype-refinement article covering original-versus-refined haplotype state, pairwise merge evidence, `effect_threshold`, multi-trait grouping, grouped and collapsed refined-variant plots, and refined phenotype visualization.
- Added a dedicated variant-effect article covering binary genotype coding, absolute versus signed effects, adjusted P values, linked tied-effect interpretation, direct single-variant phenotype validation, and multi-trait prioritization.
- Extended the pkgdown **Haplotype & association** article sequence to `Haplotype -> Phenotype -> LD -> Haplotype refinement -> Variant effect` while keeping the staged `docs/10-ld.qmd`, `docs/11-haplotype-refinement.qmd`, and `docs/12-variant-effect.qmd` sources unchanged.

# GeneTrackR 0.6.6

- Added a dedicated haplotype article covering gene-, transcript-, and region-defined `HapVariant` construction, genotype representations, missing-genotype filtering, sample/variant filtering, object inspection, and `plot_hap_variant()` visualization.
- Added a dedicated phenotype-association article covering phenotype import and QC, explicit sample-ID alignment checks, haplotype-phenotype testing, multi-trait inference, single-variant phenotype testing, and preservation of returned statistical tables.
- Added a task-oriented **Haplotype & association** section to the pkgdown Articles index, positioned after variant/browser visualization and before the later LD/refinement/variant-effect modules.
- Kept the staged `docs/08-haplotype.qmd` and `docs/09-phenotype.qmd` sources unchanged while migrating their current workflows into pkgdown vignettes.

# GeneTrackR 0.6.5

- Added a dedicated pkgdown variant-track article covering `VariantTrack` construction, memory versus indexed lazy VCF access, validation, summaries, genomic/gene/transcript retrieval, ID/type/pattern filtering, variant plotting, site-level VCF export, and downstream handoff.
- Added a dedicated integrated browser-track article covering `plot_tracks()` locator semantics, gene/transcript/region views, strand-aware RNA-seq and Ribo-seq display, transcript-level Ribo-seq frame rendering, feature and variant panels, panel layout, highlighting, styling, and large-dataset usage.
- Extended the pkgdown Visualization section from gene models and signal tracks to the complete gene-model -> signal -> variant -> integrated-browser documentation sequence.
- Kept all staged `docs/*.qmd` sources unchanged for subsequent haplotype, phenotype, LD, refinement, variant-effect, export, and end-to-end workflow migration.

# GeneTrackR 0.6.4

- Added a dedicated pkgdown gene-model visualization article covering gene, transcript, and regional structural views, genomic versus spliced transcript coordinates, transcript collapse, strand-direction arrows, labels, highlights, and explicit feature colors.
- Rebuilt the signal-track article around the complete `BwgTrack` workflow: reading, lazy versus in-memory access, strand semantics, regional retrieval, normalization, merging, binning, writing, gene/transcript/region visualization, Ribo-seq frame views, and integrated `plot_tracks()` figures.
- Separated gene-model structure from continuous-signal documentation so later variant/browser-track articles can build on a stable visualization foundation without duplicating annotation or signal basics.
- Added the gene-model article to the pkgdown Visualization section while retaining all staged `docs/*.qmd` sources unchanged for subsequent migration steps.

# GeneTrackR 0.6.3

- Added a dedicated pkgdown annotation article covering GenePred/GenePredExt, GTF, GFF3, and BED inputs through one unified GeneTrackR workflow.
- Reorganized annotation documentation around native coordinate handling, object contracts, validation and summaries, feature retrieval, cross-format conversion, merging, writing, and round-trip checks.
- Added annotation-focused visualization examples for genes, transcripts, genomic regions, interval feature tracks, and gene-length distributions.
- Added a Data & annotation section to the pkgdown Articles index while retaining the existing Getting Started and Visualization navigation.
- Kept all staged `docs/*.qmd` sources unchanged for later module-specific migration.

# GeneTrackR 0.6.2

- Rebuilt the pkgdown Get Started path around a dedicated first-run article that covers installation, package loading, bundled demo inputs, and a minimal annotation-to-haplotype/phenotype workflow.
- Added a deterministic example-data article documenting the shared file inventory, coordinate conventions, core dimensions, sample alignment rule, GeneA haplotype/phenotype truth, LD gradient, and strand-specific RNA-seq/Ribo-seq design.
- Updated the pkgdown article order to `Getting started -> Example data -> Core workflow`, while retaining signal-track documentation under Visualization.
- Added prerequisite navigation to the existing core workflow without changing its executable analysis workflow.
- Kept all staged `docs/*.qmd` sources unchanged for later module-by-module migration.

# GeneTrackR 0.6.1

- Redesigned the pkgdown information architecture around a dedicated **Get Started** entry, task-oriented **Articles**, and a functionally grouped **Reference** index.
- Organized all current Rd topics into explicit Reference sections for core objects, annotation I/O and conversion, gene-model visualization, signal tracks, variant tracks, haplotype analysis, phenotype association, linkage disequilibrium, variant effect/refinement, and advanced diagnostics.
- Positioned `core-workflows` as the top-level Get Started workflow and retained `signal-tracks` under the Visualization article section.
- Added article descriptions to the existing vignettes so the pkgdown article index explains the purpose of each workflow.
- Kept the staged `docs/*.qmd` tutorial sources unchanged; additional task-oriented articles will be migrated into the established information architecture in subsequent 0.6.x releases.

# GeneTrackR 0.6.0

- Initialized the GeneTrackR package website with `pkgdown`, replacing the abandoned MkDocs migration plan with an R-native documentation architecture.
- Added `_pkgdown.yml` with the canonical GitHub Pages URL, Bootstrap 5, the existing workflow vignettes, and a temporary `pkgdown-site/` output directory.
- Added a GitHub Actions workflow to build the package website on pushes, pull requests, releases, and manual dispatches, and to deploy non-PR builds to the `gh-pages` branch.
- Added package website and issue-tracker metadata to `DESCRIPTION`, together with the website-only `pkgdown` dependency declaration.
- Preserved the existing `docs/*.qmd` documentation sources during staged migration; they remain the source material used by `README.qmd` and are not overwritten by pkgdown in this release.

# GeneTrackR 0.5.37

- Fixed haplotype allele fill color helpers so the stable `A`, `T`, `C`, `G`, and `indel` mapping retains its names after alpha adjustment.
- Fixed variant marker color helpers so the stable `SNP`, `Ind`, and `...` mapping retains its names after alpha adjustment.
- Updated plotting-default regression tests to preserve the intentional `Reds` default for `plot_ld_block()` while other plotting palettes remain `Paired`.
- Changed full-memory bigWig loading to emit one memory warning per `read_bwg()` call instead of one warning per input file, and made round-trip tests explicitly expect this warning.

# GeneTrackR 0.5.36

- Rebuilt `docs/13-export.qmd` as a structured export workflow for genomic tracks, PDF/PNG figures, analysis tables, optional Excel workbooks, and complete RDS objects.
- Rebuilt `docs/14-workflow.qmd` as the end-to-end GeneTrackR candidate-locus workflow using the current GeneA haplotype, phenotype, LD-gradient, refinement, and variant-effect examples.
- Removed `docs/15-notes.qmd` and incorporated its large-data and object-return guidance into the end-to-end workflow.
- Synchronized `README.qmd` and `README.md` with the revised documentation sequence.

# GeneTrackR 0.5.35

- Fixed `plot_variant_effect()` multi-trait reshaping so mixed integer and double phenotype columns no longer trigger `data.table::melt()` coercion warnings.
- Standardized selected numeric traits to double before long-format conversion without changing their values or effect statistics.
- Added regression coverage requiring the documented multi-trait variant-effect workflow to run silently.
- Clarified mixed numeric trait handling in `docs/12-variant-effect.qmd`.
- Bumped the package version to 0.5.35.

# GeneTrackR 0.5.34

- Rebuilt `docs/12-variant-effect.qmd` as a complete regional variant-effect prioritization workflow.
- Added explicit guidance for binary genotype coding, absolute versus signed effects, adjusted P values, and multi-trait interpretation.
- Documented the deterministic GeneA effect truth: the p13-linked variant set has an effect of approximately +6 on `protein_content`, `varA05` has an intermediate effect, and `varA02` is a zero-effect control.
- Clarified that tied phenotype effects among linked variants do not identify a unique causal variant and should be interpreted together with LD, annotation, direct variant-phenotype plots, and refinement.
- Added regression coverage for the redesigned variant-effect documentation workflow.
- Bumped the package version to 0.5.34.

# GeneTrackR 0.5.33

- Replaced `docs/10-ld.qmd` with the user-updated LD workflow while restoring valid executable QMD syntax.
- Redesigned the GeneA genotype/phenotype truth so refinement follows genotype-nearest haplotype pairs.
- Changed `varA03` from the cross-cluster `p23` pattern to the cluster-consistent `p13` pattern.
- Redesigned phenotype means so Hap1/Hap2 are locally similar, Hap3/Hap4 are locally similar, and cross-cluster differences are larger.
- Updated the `protein_content` refinement truth to `Hap1 + Hap2` versus `Hap3 + Hap4`.
- Updated phenotype, refinement, variant-effect, workflow, vignette, generator, validator, help examples, and regression tests.
- Bumped the package version to 0.5.33.

# GeneTrackR 0.5.32

- Expanded the primary LD demo from 6 to 12 natural SNPs while preserving the 11-variant GeneA haplotype truth.
- Redesigned the primary LD demo to contain a perfect high-LD core plus progressively decorrelated downstream variants, producing a clear r2 gradient with both high- and low-LD pairs.
- Changed the default `plot_ld_block()` heatmap palette from `Paired` to `Reds`.
- Updated LD documentation, demo truth metadata, generated VCF data, and regression tests.

# GeneTrackR 0.5.31

- Rebuilt `docs/10-ld.qmd` as a complete LD calculation, inspection, visualization, and downstream handoff workflow.
- Added deterministic high-LD, low-LD, and two-variant demo examples with explicit expected pairwise behavior.
- Documented the `LDTrack` object contract, pairwise statistics, LD matrix, genotype retention, sample/variant filtering, and `plot_ld_block()` return semantics.
- Clarified that GeneTrackR calculates pairwise LD but does not automatically call biological block boundaries from a fixed threshold.
- Added guidance on dosage-based `Dprime` limitations for unphased VCF data and recommended `r2` as the primary generic metric.
- Added documentation regression coverage for the designed LD truth sets.

# GeneTrackR 0.5.30

- Rebuilt `docs/11-haplotype-refinement.qmd` as a phenotype-guided refinement workflow.
- Changed the primary demo refinement trait to `protein_content`, which deterministically refines four GeneA haplotypes into two phenotype groups.
- Added explicit inspection of pairwise tests, merge decisions, trait-specific groups, original-to-refined mappings, grouped/collapsed refined variant plots, and refined phenotype plots.
- Added guidance on `effect_threshold`, trait dependence, negative-control refinement, and multi-trait refinement.
- Added documentation regression coverage for the designed 4 -> 2, 4 -> 4, and 4 -> 1 refinement outcomes.

# GeneTrackR 0.5.29

- Reorganized `docs/09-phenotype.qmd` into a self-contained phenotype-association workflow covering phenotype QC, sample matching, haplotype association, multiple testing, single-variant association, result objects, and downstream reuse.
- Documented the designed demo phenotype roles for `seed_weight`, `plant_height`, `protein_content`, `flowering_time`, and `flower_color`.
- Clarified numeric versus categorical phenotype support and the relationship among raw p-values, adjusted p-values, and displayed significance brackets.
- Added documentation regression checks for sample-ID matching and the phenotype workflow.

# GeneTrackR 0.5.28

- Reorganized `docs/08-haplotype.qmd` into a self-contained, stepwise haplotype workflow from VCF input through downstream analysis handoff.
- Documented the `HapVariant` object structure, genotype representations, missing-genotype filtering, flanking-region semantics, sample/variant filtering, and haplotype plotting controls.
- Established the GeneA gene body as the primary tutorial haplotype definition (11 variants, 36 samples, four balanced haplotypes) and separated the upstream missing-genotype example from the main analysis.

# GeneTrackR 0.5.27

- Added `variant_palette` and `variant_colors` controls to `plot_tracks()` and forwarded them to integrated variant panels.
- Standardized all plotting palette defaults to `Paired`, including signal, frame, feature, variant, LD, haplotype-table, and phenotype fill palettes.
- Updated browser-track documentation and regression tests for variant color control and package-wide palette defaults.

# GeneTrackR 0.5.26

- Improved browser-track feature legends by adding automatic compact feature grouping in `plot_tracks()` and `plot_feature_track()`.
- Added transcript-centered automatic Ribo-seq frame rendering in `plot_tracks()` via `ribo_signal_type = "auto"`.
- Updated `docs/07-browser-tracks.qmd` to explain feature legend grouping and transcript-centered Ribo-seq frame behavior.

# GeneTrackR 0.5.25

- Reorganized `docs/07-browser-tracks.qmd` into a self-contained browser-track workflow covering input preparation, gene/transcript/region locators, signal selection, feature/variant integration, layout, heights, highlighting, and large-data usage.
- Clarified the exact-one-locator rule used by `plot_tracks()` and the distinction between browser-level genomic signal plotting and transcript-coordinate/frame-specific signal visualization.
- Updated browser examples to use the current public `read_vcf()` API and explicit strand/sample selection.
- Documented current panel ordering, named feature/variant list behavior, and the compact styling scope of integrated feature/variant tracks.

# GeneTrackR 0.5.24

- Fix fixed-string pattern filtering so `fixed = TRUE` no longer passes the incompatible `ignore.case` argument to `grepl()`.
- Preserve case-insensitive fixed-string matching by normalizing both the pattern and candidate fields before matching.
- Add regression tests for warning-free fixed matching and case-insensitive fixed-string VCF retrieval.

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
