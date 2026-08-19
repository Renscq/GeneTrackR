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
