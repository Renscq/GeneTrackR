# GeneTrackR Bioconductor submission layout

GeneTrackR development and formal Bioconductor submission use two different Git trees. The 2026 submission service reads the GitHub repository **default branch only**, so that branch must contain package code and not development infrastructure.

## Development workspace

Keep the full development workspace on a non-default development branch. It may contain:

- `.github/` pkgdown workflows;
- `README.qmd`, `docs/`, and `_pkgdown.yml`;
- `tools/` preflight and README-generation scripts;
- `dev/` development notes;
- local IDE project files that are ignored by Git.

## Package-only submission branch

GeneTrackR 0.99.0 is the formal submission candidate. Prepare the package-only tree with:

```bash
Rscript tools/bioc_preflight.011.R .
Rscript tools/prepare_bioc_submission.005.R . ../GeneTrackR-bioc-submission
```

The package-only tree contains `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `README.md`, `R/`, `man/`, `inst/`, `tests/`, and `vignettes/`, plus one top-level `.gitignore`. It must not contain `.github/`, `dev/`, `docs/`, `tools/`, `README.qmd`, `_pkgdown.yml`, `.Rbuildignore`, `pkgdown-site/`, or an IDE project file.

## Current Bioconductor devel environment

The formal candidate targets Bioconductor 3.24 devel and R 4.6.0. Use a clean R 4.6 installation and then:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "devel")
BiocManager::valid()
BiocManager::install(c("BiocCheck", "devtools"))
```

Do the final checks in the package-only Git clone:

```r
devtools::document()
devtools::test()
devtools::check()

BiocCheck::BiocCheckGitClone(
  ".",
  `quit-with-status` = FALSE
)
```

Also build a source tarball and run BiocCheck on the tarball so package-size checks are included:

```bash
R CMD build .
R CMD check --as-cran GeneTrackR_0.99.0.tar.gz
```

```r
BiocCheck::BiocCheck(
  "GeneTrackR_0.99.0.tar.gz",
  `new-package` = TRUE,
  `quit-with-status` = FALSE
)
```

## Accounts and repository requirements

Before opening the submission issue:

1. register on the Bioconductor Support Site using the maintainer identity/email used by GeneTrackR;
2. subscribe to the `bioc-devel` mailing list or follow the Bioconductor package-submission community channel;
3. ensure the submitting GitHub user is also the maintainer listed in `DESCRIPTION`;
4. add at least one SSH public key to that GitHub account;
5. make the package-only branch the GitHub repository default branch;
6. ensure the default branch contains GeneTrackR 0.99.0 and only package code.

## Formal submission

Starting in 2026, new packages are submitted through the `Bioconductor/BiocContributions` GitHub repository. Open a new submission issue, use `GeneTrackR` as the issue title, provide the GitHub repository URL, and confirm the submission/maintenance checkboxes in the issue template.

Opening the issue triggers the 2026 automated validation. If validation passes, read the policy comment and reply exactly `/accept-policies`. The submission service then clones the package into the Bioconductor staging organization and registers it with the submission R-universe. Follow the issue comment to add/push to that staging repository; subsequent review builds are triggered from the staging location, not from the original personal GitHub repository. Every build-triggering revision must increment only the patch component (`0.99.1`, `0.99.2`, and so on).

Only after acceptance is the package moved to the canonical `git.bioconductor.org` devel repository and a BiocCredentials account created for SSH-based maintenance.

No GitHub branch changes or remote writes are performed by GeneTrackR preparation scripts.
