# Standard-reader assessment for GeneTrackR 0.99.8

The follow-up [0.99.9 assessment](0.99.9-assessment.md) covers TxDb membership
mapping and an indexed VCF workload with process peak memory.
The [0.99.11 assessment](0.99.11-assessment.md) uses a public indexed VCF and
GENCODE GFF3 and includes cross-platform reproduction commands.

This is developer evidence for preliminary Bioconductor review, not a runtime
backend switch or a claim that existing wrappers replace Bioconductor classes.
All fixtures are synthetic. The binary BigWig input is the independent existing
`tests/testthat/fixtures/native_reader_test.bigwig` fixture.

## Reproduction and scope

Install the matching GeneTrackR version plus the optional development tools
`rtracklayer`, `Bioc.gff` and `VariantAnnotation` in a development library, then run:

```sh
Rscript tools/compare_standard_io.R --path . --output-dir temp/standard-io --iterations 3
```

The recorded assessment uses Windows, R 4.6.1, rtracklayer 1.73.0,
Bioc.gff 1.3.0 and VariantAnnotation 1.59.4. The original input files,
metadata/source checksums and exact projected fields are retained in the reports.
The versioned `0.99.7-*` and `0.99.8-*` CSV/DCF files are compact evidence; raw RDS objects and
full session reports remain local validation artifacts, outside the source tarball.

Annotation projections use 1-based closed coordinates and normalize only the
documented missing-value (`.`/empty to NA) and unstranded (`.` to `*`) spellings.
BED thick coordinates and block starts are compared in the source BED convention
(0-based half-open/relative offsets). Numeric comparison has an absolute 1e-6
tolerance; NA never equals zero. Composite keys must be unique and nonmissing.
Every unmatched record is written separately, and equality requires both zero
unmatched records and zero field differences. Original rows and reader objects
are preserved; no silent intersection is used.

VCF projections compare site fields, the declared DP/AF/END/FLAG INFO fields and
GT/AD/PL per sample. INFO/FORMAT values are rendered for comparison only; the
assessment adds no typed VCF adapter to GeneTrackR. Samples remain in input
order (`S2`, `S1`). PL rectangular padding is retained as an observed difference.
The projection is specific to these declared synthetic fields, not a generic
normalizer for arbitrary VCF headers.

## Observed results

| Input | Standard reader | Matched records | Unmatched records | Result on recorded fields |
|---|---|---:|---:|---|
| GFF3 | rtracklayer | 8/8 | 0 | Equivalent complete Parent metadata |
| GFF3 | Bioc.gff | 8/8 | 0 | Equivalent complete Parent metadata |
| GTF | rtracklayer | 8/8 | 0 | Equivalent |
| GTF | Bioc.gff | 8/8 | 0 | Equivalent |
| BED6 | rtracklayer | 3/3 | 0 | Equivalent |
| BED12 | rtracklayer | 2/2 | 0 | Equivalent block/thick metadata |
| VCF | VariantAnnotation | 4/4 | 0 | Three representation differences |
| BigWig | rtracklayer | 4/4 | 0 | Equivalent |

These outcomes mean equality only on the fields listed in
`0.99.8-comparison.csv`; they do not establish whole-object equivalence. The
unchanged fixtures and archived 0.99.7 results provide the before comparison:
two one-value GFF3 discrepancies and eight missing BED12 values are resolved.

- **GFF3:** the shared exon has `Parent=t1,t2`. The normalized `parent_id` still
  keeps `t1`, while the new `parent_ids` list preserves both identifiers. The
  projection now compares that complete list to the standard Parent field.
  Original attributes remain intact and regenerated attributes retain both
  links. The primary hierarchy does not expand records or implement a graph;
  lossy primary-model retrieval is guarded, with a GRanges filtering route.
- **BED12:** flat intervals now retain block counts/widths/relative offsets,
  thick spans and RGB. The projection reads the newly retained values instead
  of filling them with NA. Package tests additionally cover RGB, empty thick
  spans, repeated names, round trips and malformed/edited metadata. The standard
  comparison fields remain the same; RGB/count are not added to its scope.
- **VCF:** a missing ID is generated differently; a deletion with REF `AT` has
  a point end in VariantTrack but a two-base range in VariantAnnotation; the
  typed PL array pads a haploid likelihood pair with a missing third slot.
  The other tested INFO/FORMAT values, sample identities and raw phased/missing
  genotypes agree after the explicitly described projection. Header declarations,
  arbitrary Number=A/R/G fields and genotype-preserving writes are not solved.
- **Windows paths:** a direct absolute `E:/...` BigWig path was interpreted as
  a UCSC URL protocol in this environment. Using the fixture directory and its
  basename succeeds. The helper uses that path form for standard imports and
  restores the working directory even after errors. No package parser changes.

All 14 format/reader executions succeed in the recorded run. The native BigWig
full-memory warning is retained in `0.99.8-reader-status.csv`. It is not an error
or a reason to suppress warnings in package code.

## Resource observations and remaining qualification

Three raw reads per format/reader are recorded with warm package namespaces.
Garbage collection runs before timing; returned R-object bytes are recorded
with `object.size()`. RDS serialization and canonical projection are outside the
read timing. These tiny fixtures, warm OS caches and different S3/S4 object
contents cannot establish speed or memory superiority. Returned object size is
not peak process memory, allocation volume or a large-file benchmark.

Dependency versions and raw schema inventories are recorded. This run is Windows
only; Linux/macOS standard-reader qualification, large-file peak memory,
compressed/indexed region queries, complete attribute escaping/header fidelity
and GenomicFeatures hierarchy mapping remain open. Until those contracts are
settled, keep the existing defaults and evaluate opt-in adapters separately.

The release preserves identifier tokens rather than globally decoding/reencoding
all GFF3 attributes. BED header/comment/extra columns and general VCF header
fidelity are outside the structural contract. See the
[GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
and [UCSC BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).

Sources: [rtracklayer](https://bioconductor.org/packages/rtracklayer/),
[Bioc.gff](https://bioconductor.org/packages/Bioc.gff/),
[VariantAnnotation](https://bioconductor.org/packages/VariantAnnotation/),
[Bioconductor reuse guidance](https://contributions.bioconductor.org/reusebioc.html).
