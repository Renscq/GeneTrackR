# Changelog

## GeneTrackR 0.6.0

- Initialized the GeneTrackR package website with `pkgdown`, replacing
  the abandoned MkDocs migration plan with an R-native documentation
  architecture.
- Added `_pkgdown.yml` with the canonical GitHub Pages URL, Bootstrap 5,
  the existing workflow vignettes, and a temporary `pkgdown-site/`
  output directory.
- Added a GitHub Actions workflow to build the package website on
  pushes, pull requests, releases, and manual dispatches, and to deploy
  non-PR builds to the `gh-pages` branch.
- Added package website and issue-tracker metadata to `DESCRIPTION`,
  together with the website-only `pkgdown` dependency declaration.
- Preserved the existing `docs/*.qmd` documentation sources during
  staged migration; they remain the source material used by `README.qmd`
  and are not overwritten by pkgdown in this release.

## GeneTrackR 0.5.37

- Fixed haplotype allele fill color helpers so the stable `A`, `T`, `C`,
  `G`, and `indel` mapping retains its names after alpha adjustment.
- Fixed variant marker color helpers so the stable `SNP`, `Ind`, and
  `...` mapping retains its names after alpha adjustment.
- Updated plotting-default regression tests to preserve the intentional
  `Reds` default for
  [`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
  while other plotting palettes remain `Paired`.
- Changed full-memory bigWig loading to emit one memory warning per
  [`read_bwg()`](https://renscq.github.io/GeneTrackR/reference/read_bwg.md)
  call instead of one warning per input file, and made round-trip tests
  explicitly expect this warning.

## GeneTrackR 0.5.36

- Rebuilt `docs/13-export.qmd` as a structured export workflow for
  genomic tracks, PDF/PNG figures, analysis tables, optional Excel
  workbooks, and complete RDS objects.
- Rebuilt `docs/14-workflow.qmd` as the end-to-end GeneTrackR
  candidate-locus workflow using the current GeneA haplotype, phenotype,
  LD-gradient, refinement, and variant-effect examples.
- Removed `docs/15-notes.qmd` and incorporated its large-data and
  object-return guidance into the end-to-end workflow.
- Synchronized `README.qmd` and `README.md` with the revised
  documentation sequence.

## GeneTrackR 0.5.35

- Fixed
  [`plot_variant_effect()`](https://renscq.github.io/GeneTrackR/reference/plot_variant_effect.md)
  multi-trait reshaping so mixed integer and double phenotype columns no
  longer trigger
  [`data.table::melt()`](https://rdrr.io/pkg/data.table/man/melt.data.table.html)
  coercion warnings.
- Standardized selected numeric traits to double before long-format
  conversion without changing their values or effect statistics.
- Added regression coverage requiring the documented multi-trait
  variant-effect workflow to run silently.
- Clarified mixed numeric trait handling in
  `docs/12-variant-effect.qmd`.
- Bumped the package version to 0.5.35.

## GeneTrackR 0.5.34

- Rebuilt `docs/12-variant-effect.qmd` as a complete regional
  variant-effect prioritization workflow.
- Added explicit guidance for binary genotype coding, absolute versus
  signed effects, adjusted P values, and multi-trait interpretation.
- Documented the deterministic GeneA effect truth: the p13-linked
  variant set has an effect of approximately +6 on `protein_content`,
  `varA05` has an intermediate effect, and `varA02` is a zero-effect
  control.
- Clarified that tied phenotype effects among linked variants do not
  identify a unique causal variant and should be interpreted together
  with LD, annotation, direct variant-phenotype plots, and refinement.
- Added regression coverage for the redesigned variant-effect
  documentation workflow.
- Bumped the package version to 0.5.34.

## GeneTrackR 0.5.33

- Replaced `docs/10-ld.qmd` with the user-updated LD workflow while
  restoring valid executable QMD syntax.
- Redesigned the GeneA genotype/phenotype truth so refinement follows
  genotype-nearest haplotype pairs.
- Changed `varA03` from the cross-cluster `p23` pattern to the
  cluster-consistent `p13` pattern.
- Redesigned phenotype means so Hap1/Hap2 are locally similar, Hap3/Hap4
  are locally similar, and cross-cluster differences are larger.
- Updated the `protein_content` refinement truth to `Hap1 + Hap2` versus
  `Hap3 + Hap4`.
- Updated phenotype, refinement, variant-effect, workflow, vignette,
  generator, validator, help examples, and regression tests.
- Bumped the package version to 0.5.33.

## GeneTrackR 0.5.32

- Expanded the primary LD demo from 6 to 12 natural SNPs while
  preserving the 11-variant GeneA haplotype truth.
- Redesigned the primary LD demo to contain a perfect high-LD core plus
  progressively decorrelated downstream variants, producing a clear r2
  gradient with both high- and low-LD pairs.
- Changed the default
  [`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
  heatmap palette from `Paired` to `Reds`.
- Updated LD documentation, demo truth metadata, generated VCF data, and
  regression tests.

## GeneTrackR 0.5.31

- Rebuilt `docs/10-ld.qmd` as a complete LD calculation, inspection,
  visualization, and downstream handoff workflow.
- Added deterministic high-LD, low-LD, and two-variant demo examples
  with explicit expected pairwise behavior.
- Documented the `LDTrack` object contract, pairwise statistics, LD
  matrix, genotype retention, sample/variant filtering, and
  [`plot_ld_block()`](https://renscq.github.io/GeneTrackR/reference/plot_ld_block.md)
  return semantics.
- Clarified that GeneTrackR calculates pairwise LD but does not
  automatically call biological block boundaries from a fixed threshold.
- Added guidance on dosage-based `Dprime` limitations for unphased VCF
  data and recommended `r2` as the primary generic metric.
- Added documentation regression coverage for the designed LD truth
  sets.

## GeneTrackR 0.5.30

- Rebuilt `docs/11-haplotype-refinement.qmd` as a phenotype-guided
  refinement workflow.
- Changed the primary demo refinement trait to `protein_content`, which
  deterministically refines four GeneA haplotypes into two phenotype
  groups.
- Added explicit inspection of pairwise tests, merge decisions,
  trait-specific groups, original-to-refined mappings, grouped/collapsed
  refined variant plots, and refined phenotype plots.
- Added guidance on `effect_threshold`, trait dependence,
  negative-control refinement, and multi-trait refinement.
- Added documentation regression coverage for the designed 4 -\> 2, 4
  -\> 4, and 4 -\> 1 refinement outcomes.

## GeneTrackR 0.5.29

- Reorganized `docs/09-phenotype.qmd` into a self-contained
  phenotype-association workflow covering phenotype QC, sample matching,
  haplotype association, multiple testing, single-variant association,
  result objects, and downstream reuse.
- Documented the designed demo phenotype roles for `seed_weight`,
  `plant_height`, `protein_content`, `flowering_time`, and
  `flower_color`.
- Clarified numeric versus categorical phenotype support and the
  relationship among raw p-values, adjusted p-values, and displayed
  significance brackets.
- Added documentation regression checks for sample-ID matching and the
  phenotype workflow.

## GeneTrackR 0.5.28

- Reorganized `docs/08-haplotype.qmd` into a self-contained, stepwise
  haplotype workflow from VCF input through downstream analysis handoff.
- Documented the `HapVariant` object structure, genotype
  representations, missing-genotype filtering, flanking-region
  semantics, sample/variant filtering, and haplotype plotting controls.
- Established the GeneA gene body as the primary tutorial haplotype
  definition (11 variants, 36 samples, four balanced haplotypes) and
  separated the upstream missing-genotype example from the main
  analysis.

## GeneTrackR 0.5.27

- Added `variant_palette` and `variant_colors` controls to
  [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
  and forwarded them to integrated variant panels.
- Standardized all plotting palette defaults to `Paired`, including
  signal, frame, feature, variant, LD, haplotype-table, and phenotype
  fill palettes.
- Updated browser-track documentation and regression tests for variant
  color control and package-wide palette defaults.

## GeneTrackR 0.5.26

- Improved browser-track feature legends by adding automatic compact
  feature grouping in
  [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
  and
  [`plot_feature_track()`](https://renscq.github.io/GeneTrackR/reference/plot_feature_track.md).
- Added transcript-centered automatic Ribo-seq frame rendering in
  [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
  via `ribo_signal_type = "auto"`.
- Updated `docs/07-browser-tracks.qmd` to explain feature legend
  grouping and transcript-centered Ribo-seq frame behavior.

## GeneTrackR 0.5.25

- Reorganized `docs/07-browser-tracks.qmd` into a self-contained
  browser-track workflow covering input preparation,
  gene/transcript/region locators, signal selection, feature/variant
  integration, layout, heights, highlighting, and large-data usage.
- Clarified the exact-one-locator rule used by
  [`plot_tracks()`](https://renscq.github.io/GeneTrackR/reference/plot_tracks.md)
  and the distinction between browser-level genomic signal plotting and
  transcript-coordinate/frame-specific signal visualization.
- Updated browser examples to use the current public
  [`read_vcf()`](https://renscq.github.io/GeneTrackR/reference/read_vcf.md)
  API and explicit strand/sample selection.
- Documented current panel ordering, named feature/variant list
  behavior, and the compact styling scope of integrated feature/variant
  tracks.

## GeneTrackR 0.5.24

- Fix fixed-string pattern filtering so `fixed = TRUE` no longer passes
  the incompatible `ignore.case` argument to
  [`grepl()`](https://rdrr.io/r/base/grep.html).
- Preserve case-insensitive fixed-string matching by normalizing both
  the pattern and candidate fields before matching.
- Add regression tests for warning-free fixed matching and
  case-insensitive fixed-string VCF retrieval.

## GeneTrackR 0.5.23

- Fix
  [`summary_vcf()`](https://renscq.github.io/GeneTrackR/reference/summary_vcf.md)
  so an in-memory `VariantTrack` can be summarized without supplying a
  genomic region.
- Fix
  [`retrieve_vcf()`](https://renscq.github.io/GeneTrackR/reference/retrieve_vcf.md)
  so region, gene, and transcript filters are optional for
  in-memory/full-file queries; ID, type, and pattern filters can now be
  used independently as documented.
- Add validation for incomplete direct regions and flank arguments used
  without a gene/transcript locator.
- Add regression tests for whole-object summaries and non-regional VCF
  retrieval.

## GeneTrackR 0.5.22

- Reorganized `docs/06-variant.qmd` into a continuous VCF workflow
  covering reading, validation, summaries, region/gene/transcript
  retrieval, ID/type/pattern filtering, and plotting.
- Clarified `VariantTrack` versus `data.table` return semantics for
  [`retrieve_vcf()`](https://renscq.github.io/GeneTrackR/reference/retrieve_vcf.md).
- Added explicit guidance for indexed lazy VCF access and strand-aware
  gene/transcript retrieval.
- Documented the current site-level scope of
  [`write_vcf()`](https://renscq.github.io/GeneTrackR/reference/write_vcf.md)
  to avoid implying genotype-preserving round-trip export.

## GeneTrackR 0.5.21

- Changed discrete signal palette assignment so sample, group, and frame
  levels use RColorBrewer class colors in strict level order instead of
  interpolating between palette endpoints; this also fixes
  [`plot_signal_region()`](https://renscq.github.io/GeneTrackR/reference/plot_signal_region.md)
  color ordering.
- Kept heatmap palettes continuous while separating continuous-gradient
  generation from discrete signal color assignment.
- Reworked demo Ribo-seq counts so frame 0, frame 1, and frame 2 have
  visibly heterogeneous heights and different occupied codons while
  preserving approximately 80%/10%/10% total counts.
- Increased frame-0 codon occupancy while reducing off-frame occupancy,
  producing a dense but visually varied RPF pattern; initiation and
  termination frame-0 counts remain approximately two-fold over the
  internal mean.
- Added regression checks for discrete palette order and within-frame
  Ribo-seq count variability.

## GeneTrackR 0.5.20

- Fixed signal qualitative palettes so discrete levels use RColorBrewer
  colors strictly in palette order; `frame0`, `frame1`, and `frame2` now
  map to the first three `Set1` colors by default.
- Reworked demo Ribo-seq counts to moderate density with approximately
  80%/10%/10% frame0/frame1/frame2 total-count proportions.
- Reduced demo Ribo-seq initiation and termination peaks to
  approximately two times the internal frame-0 mean.
