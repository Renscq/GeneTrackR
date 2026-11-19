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
Rscript tools/check_coverage.R --output-dir temp/coverage --minimum 80 --baseline tools/coverage/0.99.5-files.csv
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
Current release measurements are stored under `tools/coverage/0.99.6-*`, with
a per-file comparison in `0.99.5-to-0.99.6.csv` and source/test checksums.
The Windows R 4.6.1 / covr 4.0.0 tests-only measurement is:

| Source | Covered / eligible lines | Coverage |
|---|---|---|
| Untouched 0.99.5 | 8556 / 11455 | 74.69% |
| 0.99.6 with contract regressions and length-summary correction | 9184 / 11456 | 80.17% |

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
0.99.6 does not replace this backend or claim a lossless VCF round trip.

## Complete package-check CI

```sh
Rscript tools/check_package.R --path . --output-dir temp/package-check
```

Use a fresh output directory. The helper runs the source-repository checker,
full R CMD build/check (including all vignettes and the manual), and new-package
BiocCheck, retaining logs and structured severity results. It fails on any
BiocCheck ERROR/WARNING or an R CMD check result other than Status: OK.
It requires the development dependency BiocCheck and an available Pandoc/LaTeX
toolchain. For local release staging, pass the verified staging directory as
`--path`; automatically generated source-check reports can be adjacent to that
source directory.

`check-package.yaml` configures Linux in the official Bioconductor devel image
so the R/Bioconductor pair is supplied by that image. Windows/macOS jobs use R
release with matching Bioconductor dependencies. `test-coverage.yaml` uses the
same devel image for tests-only coverage. These jobs upload evidence only and
have read-only repository permissions; they contain no publishing, pushing or
BiocStaging steps. The separate existing pkgdown publication workflow is
unchanged. Configuration and local commands can be checked without executing
GitHub jobs; remote CI results require a later authorized push.
