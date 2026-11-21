# GeneTrackR 0.99.10 preliminary review response draft

The changes below address the November 13 preliminary review.
[Review comment](https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197)

## DESCRIPTION and public API

We removed `roxygen2` from `Suggests`. The optional funding role `fnd` is
inapplicable to this package. `Feature()` and `read_vcf()` are the recommended
entry points in the README, tutorials and reference. The existing aliases
remain available for downstream callers.

## Vignettes and workflow integration

We consolidated 15 vignettes into six numbered BiocStyle workflows. The
introduction compares GeneTrackR with Gviz, rtracklayer and VariantAnnotation.
The annotation workflow starts from an existing `GRanges`, performs a
GeneTrackR analysis step and converts back with row metadata preserved.
[Introduction](../vignettes/getting-started.Rmd) · [Annotation workflow](../vignettes/annotation.Rmd)

## Standard classes and I/O

The optional transcript-grouped `GRangesList` preserves all exon and CDS
memberships while retaining the existing primary-parent statistics and saved
objects. In the TxDb fixture, all three exon and one CDS memberships match
after explicit GFF3 identifier mapping. A shared exon occurs once per source
row within each transcript. [TxDb counts](io-evaluation/0.99.10-txdb-counts.csv)

The bidirectional `CollapsedVCF` adapters use VariantAnnotation for typed VCF
metadata and genotype data. The export workflow demonstrates standard writing
with `VariantAnnotation::writeVcf()`, which can round Float FORMAT values.
The existing GeneTrackR writer continues to export the eight site columns.
[VCF adapter](../man/as_variant_track.Rd) · [Export workflow](../vignettes/export-results.Rmd)

The annotation workflow explains how the existing classes and native I/O
support compatibility with saved objects, analysis conventions and defaults.
It documents the standard GRanges and VCF boundaries and indexed use of
Rsamtools. We compared annotation, signal and VCF fixtures with standard
readers. The indexed VCF comparison matched 50 of 50 expected POS records
across eight samples on a Windows synthetic fixture. Large real files and
Linux and macOS benchmarks are the next evidence needed for a default reader
decision. [Class and I/O rationale](../vignettes/annotation.Rmd#design-decisions) · [Backend assessment](io-evaluation/0.99.9-assessment.md)

## Coverage and validation

On the same Windows tests-only coverage scope, line coverage rose from 82.59%
in 0.99.9 to 82.77% in 0.99.10. The release checks recorded 1,362 passing
package assertions and 121 passing developer tool assertions. `R CMD check`
reported `Status: OK` for the release source archive. The standalone BiocCheck
rerun on the same archive reported 0 errors, 0 warnings and 6 notes and
confirmed Support Site registration.
[0.99.10 validation record](review-progress-0.99.10.md) · [Coverage summary](coverage/0.99.10-summary.csv)
