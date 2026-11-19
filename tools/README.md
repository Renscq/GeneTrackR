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
`pkgdown-site/index.html`. All 15 source tutorials remain in `vignettes/`;
`_pkgdown.yml` controls their website navigation. Generated site files are
ignored by Git and excluded from the source package.

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
