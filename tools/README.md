# README figures and the package website

These developer tools do not install dependencies or publish anything locally.
Run all commands from the GeneTrackR package root with R 4.6.0 or later.

## Regenerate the README

Install the package dependencies and the developer packages `pkgload`, `ragg`,
`systemfonts` and `png` in your development environment. Install Arial before
regenerating the gallery; the generator rejects a missing Arial font rather
than substituting another font.

```sh
Rscript tools/render_readme.R
Rscript tools/render_readme.R --check
```

Edit `README.qmd`; do not edit the generated `README.md`. The QMD is a Markdown
source, with display-only R code fences. `render_readme.R` strips its YAML
front matter and separately executes the three figure workflows against this
checkout via `pkgload::load_all()`. It does **not** invoke Quarto or execute the
installation examples.

The generator uses only `inst/extdata/gtr_demo_*` files, checks sample IDs,
missingness, group sizes and variant coordinates, and saves the three gallery
PNGs to `man/figures/` at 500 dpi. New gallery text is regular black Arial;
continuous plots use classic styling. The haplotype diagram retains the
package's specialized table and connector layout. The random seed is 1014.

Commit the source, generated README and figures together when publication is
authorized. Check mode is read-only and fails when the Markdown is stale or a
PNG is missing or cannot be decoded, including a missing PNG end marker.
CI installs `png` through `Config/Needs/website` and checks existing figures
rather than regenerating them, so build runners do not need an Arial installation.

## Build a local website

Install this checkout and the developer dependency `pkgdown` in your chosen
development library, then run the following in R:

```r
Sys.setenv(GENETRACKR_PKGDOWN = "true")
pkgdown::check_pkgdown()
pkgdown::build_site_github_pages(
  dest_dir = "pkgdown-site",
  new_process = FALSE,
  install = TRUE
)
```

The `install = TRUE` setting installs this source into a temporary library
before running article and reference examples. The generated homepage is
`pkgdown-site/index.html`. Six BiocStyle workflows are installed from `vignettes/`;
`_pkgdown.yml` controls their navigation and the nine retired article redirects.
Retired local vignette names are no longer installed. Generated site files are
ignored by Git and excluded from the source package.

The tutorials apply explicit Arial, regular black text to returned plots and
saved figures. Continuous plots use classic styling; haplotype/LD diagrams
retain their specialized layouts. Install Arial for matching text metrics on
another machine; successful CI execution does not verify installed font faces.
All six workflows use knitr's `ragg_png` device, and the PNG export example
uses `ragg::agg_png`. The suggested dependency `ragg` removes the tutorials'
requirement for Cairo/XQuartz while retaining sizes, aspect ratios and 500 dpi.
This changes documentation rendering only, not package plotting APIs or I/O.

The homepage uses `pkgdown/sidebar.html` for its project links. This avoids
making a local build depend on automatic package-repository availability probes.

Use the navigation menus to open local article pages during a preview. The
README tutorial table uses production URLs so that it also works on GitHub.

## Publishing configuration

`.github/workflows/pkgdown.yaml` builds pushes and pull requests targeting
`devel` or `development`, and supports manual builds. Every successful build
uploads a `pkgdown-site-preview` artifact. Only a non-pull-request run on the
repository's default branch uploads a Pages artifact and enters the deployment
job. The website build installs this checkout into a temporary library before
running its article and reference examples.

Before the first authorized deployment, select **Settings > Pages > Build and
deployment > Source > GitHub Actions** in the GitHub repository. A branch-based
Jekyll build does not execute pkgdown. If the `github-pages` environment has
deployment branch restrictions, they must allow the repository's default
branch. Local validation does not change either remote setting.

After GeneTrackR appears in a Bioconductor release, update the installation
section in `README.qmd` to lead with `BiocManager::install("GeneTrackR")`, then
regenerate the README. Until then, retain the submission-candidate wording.

## Measure tests-only coverage

Install the development dependency `covr >= 4.0.0` in your chosen development
library. It is declared under `Config/Needs/coverage`, not as a runtime import.
From the package root:

```sh
Rscript tools/check_coverage.R --output-dir temp/coverage --minimum 80 --baseline tools/coverage/0.99.7-files.csv
```

The script runs package tests only. It requests no line/function exclusions and
fails when tests fail or measured coverage is below the selected floor. The
80 percent floor is a project milestone, not a threshold specified by the
preliminary review. Use `--minimum 0` for an exploratory baseline run.

Reports include raw RDS coverage, eligible/covered line counts, per-file CSV,
optional before/after comparison, R/covr versions, scope, platform, source/test
checksums and session information. Code-free documentation files do not have
instrumented lines. The aggregation is checked against covr's own line metric.
Covr 4 changed line counting; compare results under the recorded engine and
scope. Tests-only coverage is not equivalent to coverage that also runs
vignettes/examples, and assertion counts are not coverage percentages.

The untouched 0.99.5 snapshot is recorded under `tools/coverage/0.99.5-*`.
Historical 0.99.6 measurements remain under `tools/coverage/0.99.6-*`, with
their original `0.99.5-to-0.99.6.csv` comparison. The 0.99.7 evidence remains
under `tools/coverage/0.99.7-*`, with `0.99.6-to-0.99.7.csv`. Current 0.99.8
evidence is stored separately under `tools/coverage/0.99.8-*`, with
`0.99.7-to-0.99.8.csv` and exact source/test checksums. Historical evidence
is not overwritten.
The Windows R 4.6.1 / covr 4.0.0 tests-only measurement is:

| Source | Covered / eligible lines | Coverage |
|---|---|---|
| Untouched 0.99.5 | 8556 / 11455 | 74.69% |
| 0.99.6 with contract regressions and length-summary correction | 9184 / 11456 | 80.17% |
| 0.99.7 with signal, export and variant-plot regressions | 9400 / 11456 | 82.05% |
| 0.99.8 with annotation metadata and primary-model compatibility regressions | 9577 / 11618 | 82.43% |

The three focused 0.99.7 files improve `bwg_ops.R` from 48.58% to 55.49%,
`plot_variant.R` from 55.47% to 95.31%, and `write_feature.R` from 61.27% to
93.92%. The signal module still has substantial untested branches; reaching the
overall floor is not a claim of complete coverage. These tests preserve the
existing normalization, aggregation and plotting algorithms.

CI uploads full reports;
it does not send coverage to an external reporting service. Regression tests
cover existing contracts and the length-summary correction described below.

## Length-summary correction and compatibility

Transcript-level exon/CDS/UTR length summaries now contain one row per transcript
and feature. A scalar label prevents the input `feature` column from masking
the function argument and recycling the total once per segment. This changes
the previously incorrect row counts and plot weights, including CDS/UTR with
`unit = "auto"`; distinct transcripts with the same length remain separate.
Segment records, inclusive length formulas, coordinate/strand rules, schemas,
aliases and other analysis algorithms are unchanged. The annotation workflow
checks direct summaries against explicit aggregation of segment records.

The existing VCF writer is also explicitly characterized: it exports eight
site fields with a minimal metadata header, omitting FORMAT/sample genotypes.
The export workflow retains the original VCF and complete objects separately.
Versions 0.99.6 through 0.99.8 retain this backend and do not claim a lossless
VCF round trip.

## Annotation metadata contracts in 0.99.8

GFF3 now adds `parent_ids`, a character-list column containing every original
Parent token. `parent_id` remains the first parent; primary gene-model tables
and record counts are unchanged. Regenerated GFF attributes preserve the full
list, require consistent primary/companion edits and accept legacy objects
without the companion column. Original text still wins by default with
`keep_attributes = TRUE`. This is structured metadata preservation, not a new
multi-parent analysis graph. A lossy `retrieve_feature()` primary-model rebuild
now fails for affected shared-parent records; filter feature-level GRanges and
then call `as_feature()` to retain the original records and metadata.
Merge renaming applies the existing source/ID map to every Parent token.
Internal statistics/plots/haplotype calculations use a temporary primary-model
projection with unchanged algorithms; the original annotation is not mutated.

BED12 remains flat Feature intervals with additive `block_count`, integer-list
`block_sizes`/`block_starts`, `thick_start`/`thick_end` and `item_rgb` columns.
Offsets are zero-based relative bp; sizes are bp widths; thick spans are closed
1-based coordinates, with end == start - 1 for empty spans. The new writer path
preserves structure and repeated names; existing GenePred output is unchanged.
Counts, geometry, exact integers and RGB are validated before export. BED3-11
behavior is retained. Invalid BED12 records/edits now error rather than silently
flattening their structure. Examples in the annotation workflow also demonstrate
the optional nonstandard granges mode and its explicit empty-span convention.

## Compare standard readers

The standalone `compare_standard_io.R` helper evaluates the current readers
against `rtracklayer`, `Bioc.gff` and `VariantAnnotation`. These are optional
developer dependencies, not new package Imports/Suggests. Install them in a
separate development library and install this same GeneTrackR version before
running the helper; the helper refuses a mismatched installed version.

```sh
Rscript tools/compare_standard_io.R --path . --output-dir temp/standard-io --iterations 3
Rscript -e 'testthat::test_file("tools/tests/test-compare-standard-io.R", stop_on_failure = TRUE)'
Rscript -e 'testthat::test_file("tools/tests/test-standard-io-readers.R", stop_on_failure = TRUE)'
```

The second standalone test file needs the optional reader packages; it explicitly
skips when those tools are absent. All optional-reader tests were run with the
tools present for the 0.99.7 and 0.99.8 assessments. They are not package tests and do not
contribute to the tests-only runtime coverage percentage.

Use a fresh output directory. Reports retain raw reader objects, input hashes,
projected tables, unmatched rows, field-level differences, reader warnings/errors,
source hashes, package versions and per-read timing/object-size observations.
The helper reports observed differences without converting them into successful
equivalence checks. The CSV records exactly which fields were compared.

See [the recorded assessment](io-evaluation/README.md),
[the unposted review response](review-response-0.99.8.md) and
[the progress record](review-progress-0.99.8.md). The assessment does not activate
an alternative backend or establish a performance ranking. Large-data peak
memory, indexed-query comparisons and other operating systems remain future work.

## Complete package-check CI

```sh
Rscript tools/check_package.R --path . --output-dir ../GeneTrackR-package-check
```

Use a fresh output directory outside the actual package source. The default
is a persistent sibling named `<source-directory>-package-check`; source-local
paths, including `temp/package-check` when checking the repository root, are
rejected before any checker/build log is opened. `.Rbuildignore` does not
prevent Windows from trying to copy an active source-local build log.
The helper runs the source-repository checker,
full R CMD build/check (including all vignettes and the manual), and new-package
BiocCheck, retaining logs and structured severity results. It fails on any
BiocCheck ERROR/WARNING or an R CMD check result other than Status: OK.
If all error/warning diagnostics are the recognized Support Site HTTP 429
messages, it
repeats the complete new-package BiocCheck after 30 and then 60 seconds
(at most three attempts). Each attempt retains its checker object, native
reports and severity row. This includes the ERROR used by BiocCheck 1.48 and
the WARNING used by 1.49; original severities are never downgraded.
Other errors/warnings fail without retry; exhausted
rate limiting still fails. Optional local BiocVersion metadata is recorded as
`not installed` when absent rather than adding a developer dependency.
It requires the development dependency BiocCheck and an available Pandoc/LaTeX
toolchain. For local release staging, pass the verified staging directory as
`--path`; automatically generated source-check reports can be adjacent to that
source directory.

Run the standalone report-path regressions from the package root:

```sh
Rscript -e 'testthat::test_file("tools/tests/test-check-package.R", stop_on_failure = TRUE)'
```

`check-package.yaml` and `test-coverage.yaml` keep the official Bioconductor
devel container's R with `setup-r`'s `install-r: false`. Bioconductor devel
does not always use R-devel. The jobs independently resolve the expected
devel repository, use `BiocManager::install(version = "devel")` to validate
the R/Bioconductor pair, and compare the installed version and actual
repository with that expectation before persisting it. Windows/macOS use R release and
its matching default Bioconductor release. All check/coverage jobs explicitly
use current CRAN instead of a dated container snapshot, persist the selected
Bioconductor repositories for dependency installation, and record actual R,
Bioconductor, repositories and libraries before dependencies are resolved.
Coverage additionally records/asserts the installed covr version before
measuring; `covr >= 4.0.0`, tests-only scope and the 80% floor are unchanged.

Linux package checks refresh apt metadata and record/install the texinfo
candidate before TinyTeX; PDF checks are retained. All package-check runners
use the official `TinyTeX` prebuilt bundle containing Courier, then verify
`pcrr8t.tfm` before building. This avoids a live CTAN/self-update phase;
font verification output is retained in `tex-fonts.log`. Check/coverage outputs live
under `${{ runner.temp }}` outside the checkout. Runner diagnostics are created
before toolchain setup using a step-level runner context and `GITHUB_ENV`,
and reports upload even on failure. An absent artifact
warns instead of masking the original failed check; check failures still fail
the job. Artifact uploads use the verified Node 24 action `v7.0.2`.

These jobs upload evidence only and
have read-only repository permissions; they contain no publishing, pushing or
BiocStaging steps. The pkgdown preview upload uses the same supported artifact
action; website deployment conditions are unchanged.
Configuration and local commands can be checked without executing
GitHub jobs; remote CI results require a later authorized push.
