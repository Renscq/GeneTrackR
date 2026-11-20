# Preliminary-review response draft: GeneTrackR 0.99.8

Draft for the submission issue. BiocStaging synchronization and reviewer
acceptance are pending.
Review: https://github.com/Bioconductor/BiocContributions/issues/177#issuecomment-6853763197

## Metadata and recommended public entry points

Removed roxygen2 from Suggests. The README/reference organization recommends
Feature() and read_vcf(), while keeping FeatureTrack()/read_vcf_track() as
compatibility aliases. The optional fnd role and ORCID are not invented without
verified maintainer information. Acceptance of the alias strategy remains open.

## Vignettes, comparison and Bioconductor entry workflow

Fifteen vignettes were consolidated into six numbered executable BiocStyle
workflows; navigation and old article redirects follow the same organization.
The introduction compares GeneTrackR with Gviz/rtracklayer/VariantAnnotation.
The annotation workflow demonstrates entry from existing GRanges, metadata
edits, bidirectional conversion and now GFF3/BED12 structural round trips.

## Reuse of standard classes and readers

Feature now has a bidirectional GRanges boundary. The class/I/O discussion
explains retained workflow state and remaining GenomicFeatures/VCF mappings.
Existing S3 classes, public signatures and analysis algorithms remain intact;
convenience alone is not offered as justification for duplicate infrastructure.

The shared-fixture developer assessment uses rtracklayer, Bioc.gff and
VariantAnnotation. In 0.99.7 it exposed a missing normalized GFF3 parent link
and BED12 block/thick metadata. The 0.99.8 candidate resolves those structural
metadata gaps: parent_ids retains the complete Parent list, and flat BED12
objects retain structural fields through filtering/GRanges/read-write paths.
The legacy parent_id and primary hierarchy remain unchanged. Regenerated GFF
attributes keep all parents, inconsistent edits fail, and affected lossy
retrieve_feature model rebuilds are guarded. This is not a multi-parent graph.
The guard includes selected gene-level and secondary-parent links. Merge rename
maps update every parent token. Internal statistics, plots and haplotype location
retain the established primary-model projection, with explicit equality and
input-immutability regressions; those projections are not metadata exports.

Eight reader comparisons retain all input records. GFF3/GTF/BED6/BED12/BigWig
agree on the explicitly recorded fields; three VCF representation differences
remain visible. Old evidence is not rewritten. Resource observations are small
warm-cache reads and object.size, not performance rankings or peak-memory proof.

Rsamtools is already reused for indexed queries. GenomicFeatures hierarchy
mapping, a complete VariantAnnotation::VCF adapter, genotype-preserving export,
large-data/region-query/other-platform qualification and an eventual standard
backend decision remain open. Defaults and mandatory dependencies are unchanged.
We request guidance on the acceptable next interoperability milestone.

## Coverage, checks and remaining advisories

The previous same-method covr 4 tests-only baseline is 9400/11456 = 82.05% for
0.99.7. New tests verify complete parent lists, regenerated attributes, exact
BED12 coordinates/geometry/RGB, repeated names, empty output, malformed data,
metadata edits and compatibility with the primary model. The 80% floor and
no-exclusion method remain unchanged. Current 0.99.8 evidence is stored under
tools/coverage/0.99.8-*; use its actual final measurements/checker reports when
posting, rather than the reviewer's differently measured 74.34% baseline.

After the independent-review compatibility repairs, the candidate passes 1232
package assertions (88 additional assertions) and 115 developer-tool assertions
with no failures or skips. The same-method tests-only coverage is
9577/11618 = 82.43%, versus the preserved 0.99.7 82.05% baseline. No exclusions
are requested. Primary-model tables/legacy columns for four input formats and
six legacy GenePred export byte streams match the 0.99.7 source snapshot;
six analysis files differ only by the internal primary-retrieval delegate.

The repaired candidate's final full R CMD check reports Status: OK, including
six vignette builds/rebuilds and the PDF manual. BiocCheckGitClone reports 0/0/0;
new-package BiocCheck 1.49.30 reports 0 ERROR / 0 WARNING / 6 NOTE after a complete
retry of a Support Site HTTP 429 warning. All original attempt reports are
retained. The final Watched Tags profile fetch still returned an informational
HTTP 429, so that online profile was not reconfirmed in this run. No diagnostic
severity was downgraded or removed. The checked archive/source proof is recorded
in the local 0.99.8 delivery report; official builds still await authorization.

The [devel workflow history](https://github.com/Renscq/GeneTrackR/actions?query=branch%3Adevel)
records remote checks for each commit. The local checks above cover the exact
0.99.8 source archive. Official candidate builds follow staging synchronization.

Existing BiocCheck advisories concern ORCID, warning suppression, long functions,
line length, indentation and unverifiable Bioc-Devel subscription. Missing
identity/subscription information requires maintainer verification. Broad style,
function decomposition and diagnostic refactoring remain focused future tasks;
these explanations do not imply the reviewer has waived the advisories.
