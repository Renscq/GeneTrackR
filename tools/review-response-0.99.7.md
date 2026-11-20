# Preliminary-review response draft: GeneTrackR 0.99.7

Status: local draft, not posted. The candidate has not been synchronized to
BiocStaging. Engineering implementation and test results do not imply reviewer
acceptance. Review source:
https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197

Thank you for the preliminary assessment. The changes below preserve existing
user code while improving documentation, interoperability and test evidence.

## DESCRIPTION and public entry points

- Removed roxygen2 from Suggests in 0.99.3. Documentation generation remains a
  development task.
- Clarified Feature() and read_vcf() as the recommended entry points in the
  README and website reference organization. FeatureTrack() and read_vcf_track()
  remain compatibility aliases; their removal would break existing scripts.
  We request feedback on whether this compatibility approach is acceptable.
- The optional fnd role has not been added because no verified grant metadata
  is available. No funder or author identifier is inferred.

## Vignette organization and ecosystem integration

- Consolidated fifteen vignettes into six numbered executable workflows:
  getting-started, annotation, haplotype, phenotype, complete-workflow and
  export-results. They all use BiocStyle. README/pkgdown navigation follows the
  same order and retired website article names redirect to the relevant section.
- Added the requested short Gviz/rtracklayer/VariantAnnotation comparison in
  vignettes/getting-started.Rmd. It explains complementary roles, not superiority.
- Added an existing-GRanges entry example in vignettes/annotation.Rmd with
  as_feature() and feature-level as_granges(). The round trip retains supported
  metadata, sequence information, names and typed columns; regressions also
  exercise metadata edits, empty inputs and explicit lossy-conversion rejection.

## Classes and existing Bioconductor readers

- GRanges is now a bidirectional interoperability boundary for Feature.
  vignettes/annotation.Rmd documents each class's current workflow state and
  remaining standard-container mapping. Existing S3 classes are retained during
  this compatibility-preserving milestone; a convenience wrapper alone is not
  presented as justification for a separate parser or class.
- Added tools/compare_standard_io.R, deterministic shared fixtures and recorded
  comparisons with rtracklayer, Bioc.gff and VariantAnnotation. All eight
  comparison pairs match every input record. GTF/BED6/BigWig agree on the tested
  fields. Full differences and the limited field scope are reported explicitly.
- The assessment exposes remaining GFF3 multi-parent normalization and flat
  BED12 import limitations, documented in the public reader help and tutorial.
  It also distinguishes VariantTrack point/ID conventions from standard VCF
  ranges and typed-array padding. These are not claimed to be solved by more
  tests or documentation.
- Rsamtools is already reused for indexed queries. The native backends and the
  sites-only VCF writer remain unchanged. A complete VCF-class adapter,
  GenomicFeatures hierarchy mapping, standard-backend qualification and
  genotype-preserving writing remain open. We request reviewer guidance on
  the priority and acceptable scope of those follow-up changes.

## Test coverage and validation evidence

- Added behavior-focused tests for the previously uncovered files and for
  variants, validation, indexed VCF queries and corrected length summaries.
  The 0.99.7 additions cover named signal normalization, inclusive region
  summaries, slicing, conflict aggregation, export formats and plot data/scales.
- Windows covr 4.0.0 tests-only line coverage is 9400/11456 = 82.05%, compared
  with the same-method 0.99.6 baseline of 9184/11456 = 80.17%. No exclusions are
  requested. The reviewer-reported 74.34% used a different measurement, so it
  is not treated as a directly comparable baseline.
- The existing CI runs complete package checks on Linux/Windows/macOS plus
  coverage and website checks. The last published 0.99.6 commit 8952c34 passed
  all four workflows/nine jobs. This does not claim remote validation of the
  unpushed 0.99.7 candidate. Local 0.99.7 checker evidence must accompany the
  final authorized submission.

Local candidate verification on 2026-11-20 completed successfully: full
R CMD check reported Status: OK, including all six vignette rebuilds and the
PDF manual. BiocCheck 1.49.30 reported 0 ERROR / 0 WARNING / 6 NOTE, and
BiocCheckGitClone reported 0/0/0. The package suite passed 1144 assertions with
eight existing dependency warnings and no failures/skips; the standalone
development tools passed 111 assertions with no failures/skips. A separate
read-only review found no blocking issues. The root 0.99.7 tarball is copied
from that checked archive; see the local delivery report for its SHA256.

## Existing BiocCheck advisories

The 0.99.7 local checker reports ORCID, warning/message suppression, long
functions, line length, indentation and unverifiable Bioc-Devel subscription
advisories. ORCID and subscription need maintainer-supplied verification, not
invented metadata. Function decomposition, diagnostic handling and broad style
changes are deferred to focused follow-up work with regression coverage; this
release does not reformat or restructure the analysis/I/O implementations.
These explanations are requests for review, not assertions that the notes have
been waived. Official candidate-build results still need to be linked when
posting this draft after authorized synchronization.

Before posting: confirm the candidate archive/checksums and current checker
results, verify maintainer-only metadata, obtain authorization to synchronize
BiocStaging, and then link the official candidate build. No issue reply,
staging synchronization, commit or push is performed by this document.
