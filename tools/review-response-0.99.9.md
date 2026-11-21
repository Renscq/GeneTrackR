# GeneTrackR 0.99.9 preliminary review response draft

This draft tracks the November 13 preliminary review. The 0.99.9 candidate
is local. The published CI evidence below corresponds to 0.99.8.
[Review comment](https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197)

## Package metadata, API and documentation

We removed `roxygen2` from `Suggests`. The recommended entry points are
`Feature()` and `read_vcf()`. The existing aliases remain available for
downstream callers. The reference layout groups aliases with their recommended
entry points. Adding an `fnd` role depends on verified funding information.

The original 15 vignettes were consolidated into six numbered BiocStyle
workflows. The introductory vignette compares Gviz, rtracklayer and
VariantAnnotation. The annotation workflow demonstrates starting from an
existing `GRanges` object, preserving its row metadata and converting back.
The haplotype and export workflows now demonstrate the corresponding standard
VCF entry and typed export. [Haplotype workflow](../vignettes/haplotype.Rmd) · [Export workflow](../vignettes/export-results.Rmd)

## Standard representation and I/O reuse

The `CollapsedVCF` adapters use VariantAnnotation for the typed header, ranges,
INFO and genotype data. `as_variant_track()` projects the fields needed for
existing POS-based analysis. `as_vcf()` selects source records and samples by
stable row and sample keys. Regression cases include repeated IDs and positions,
allele order, Number=A/R/G fields, Flags, variable ploidy, phased and partially
missing genotypes, filtering, sample reordering, and RDS restoration.

The existing eight-column writer serves sites-only export. Complete standard
objects go through `VariantAnnotation::writeVcf()`. Its file-formatting behavior
includes Float FORMAT rounding and `fileDate` updates. The in-memory boundary
retains exact values from the supplied typed object. The adapters reject missing
declarations, ambiguous empty ALT elements, foreign rows and edits to the
projected VCF fields. [Adapter help](../man/as_variant_track.Rd)

GFF3 retains all declared Parent links and the existing primary assignment.
BED12 retains block geometry, thick intervals and RGB in flat Feature records.
The TxDb evaluation matches all seven transcript/exon/CDS memberships after
mapping display names through the fixture's GFF3 `Name` and `ID` declarations.
The native primary assignment includes six of those memberships. The remaining
membership belongs to a shared exon and is retained in `parent_ids`.

The two indexed VCF readers match 50/50 expected POS records and the eight
samples on a 50,000-site synthetic input. The assessment records namespace load
time, query time and Windows process high-water memory. Results support retaining
the current defaults pending larger and cross-platform qualification. The
assessment identifies separate acceptance criteria for a multi-parent model,
BED12 gene-model conversion and a backend migration. [Backend assessment](io-evaluation/0.99.9-assessment.md)

## Checks and submission status

The 0.99.8 `devel` commit passed package checks on Linux, Windows and macOS,
coverage and site builds. Linux tests-only line coverage reached 82.42%.
The coverage workflow retains the same scope and 80% threshold for 0.99.9.
[Package checks](https://github.com/Renscq/GeneTrackR/actions/runs/44754728222) · [Coverage](https://github.com/Renscq/GeneTrackR/actions/runs/44754727951) · [pkgdown](https://github.com/Renscq/GeneTrackR/actions/runs/44754728005) · [Pages](https://github.com/Renscq/GeneTrackR/actions/runs/44754727623)

The current test and archive evidence is recorded in the
[0.99.9 progress table](review-progress-0.99.9.md). BiocStaging synchronization,
official build results and acceptance of the retained API aliases remain open.
