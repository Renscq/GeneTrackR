# Preliminary review response for GeneTrackR 0.99.12

Thank you for reviewing GeneTrackR. The points below describe the changes
and validation in version 0.99.12. File paths refer to the GeneTrackR source
tree. [Preliminary review](https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197)

## 1. Development dependency

We removed `roxygen2` from `Suggests`. The package builds with the generated
reference documentation and namespace supplied in the source tree.

## 2. Funding role

The optional `fnd` role does not apply to this submission.

## 3. Recommended public API

The README, reference, and tutorials direct users to `Feature()` to
construct annotations and `read_vcf()` to read variants. The compatibility
aliases `FeatureTrack()` and `read_vcf_track()` still support existing
scripts. We added `as_variant_track()` and `as_vcf()` as explicit adapters
to standard Bioconductor VCF objects.

## 4. Vignette organization

We consolidated 15 vignettes into six numbered workflows. They cover
getting started, annotation and signals, haplotypes and linkage
disequilibrium, phenotypes, a complete locus workflow, and export.

## 5. Formatting

All six vignettes use `BiocStyle::html_document`.

## 6. Comparison with existing tools

The introduction in `vignettes/getting-started.Rmd` compares Gviz for track
visualization, rtracklayer for genomic file exchange, and VariantAnnotation
for typed variants. It describes how GeneTrackR connects locus views to its
haplotype and phenotype analysis workflow.

## 7. Entry from a Bioconductor workflow

The annotation vignette starts with an existing `GRanges`, uses
`as_feature()` to enter the GeneTrackR workflow, runs an interval query,
and returns to `GRanges`. The example demonstrates how coordinates, names
and supported metadata are retained. See the existing-GRanges section in
`vignettes/annotation.Rmd`.

## 8. Standard classes

The package uses standard objects at its integration boundaries.
`as_granges()` can return flat feature ranges or an exon/CDS
`GRangesList` grouped by transcript. A shared exon can belong to every
annotated transcript, and duplicate memberships of the same source row
within a transcript are counted once. The TxDb comparison and its ID
mapping are recorded in `tools/io-evaluation/0.99.10-txdb-rows.csv`.

The bidirectional `CollapsedVCF` adapters reuse VariantAnnotation for
typed VCF and support regional and sample selection. The native S3
objects carry workflow-specific data such as legacy point-based variant
coordinates, file-backed signal metadata, and per-sample scaling and
strand. They support compatibility with saved objects and downstream
analysis scripts.

## 9. Standard I/O

Rsamtools handles indexed regional queries when selected by the existing
backend. The export tutorial demonstrates typed reading and writing
with `VariantAnnotation::readVcf()` and `writeVcf()`. The annotation
assessment includes rtracklayer and Bioc.gff comparisons on fixtures.

The 0.99.11 real-data assessment compared the 1000 Genomes chr22 indexed
VCF subset with VariantAnnotation and GENCODE v19 chr22 annotation with
rtracklayer on Windows and Linux. The VCF comparison matched 8,771
overlaps and 8,770 records by POS, with matching header order for all
five samples and matching tested genotype fields. The GFF3 comparison
matched 58,721 records and 57,458 Parent links. The inputs, field scope,
reader versions, query times and peak memory are documented in
`tools/io-evaluation/0.99.11-assessment.md`. The macOS real-data run
remains pending.

The native readers support established file and analysis conventions.
Keeping those readers available with the standard adapters supports
compatibility while the evidence for alternative defaults develops.

## 10. Coverage

Version 0.99.12 adds tests for LD statistics and object construction,
signal retrieval and merging, and the coordinates, groups and colors
used in signal and feature plots. Under the same Windows R 4.6.1 and
covr 4.0.0 tests-only line metric, coverage rose from 83.77%
(9,895/11,812 lines) in 0.99.11 to 88.02% (10,400/11,816 lines) in 0.99.12.
The covered-line count rose by 505. The per-file report in
`tools/coverage/0.99.11-to-0.99.12-files.csv` shows increases in LDTrack
construction and reporting (54.84% to 93.55%), signal operations
(68.69% to 92.38%), signal plots (67.59% to 87.47%), and feature plots
(65.91% to 96.97%). The versioned coverage reports record source and
test checksums.

## Release validation

The 0.99.12 source archive completed Windows R 4.6.1 R CMD check with
`Status: OK`. All six vignettes rebuilt and the PDF manual passed.
The test log records 1,640 passing assertions and four existing
Bioconductor dependency initialization warnings. BiocCheck 1.49.45
reported zero errors, zero warnings, and six notes for the same
archive. In the Bioconductor 3.23 Linux container, the archive
passed 447 assertions in the eight affected test files.
