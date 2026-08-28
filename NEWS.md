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
