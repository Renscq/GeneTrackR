# GeneTrackR Bioconductor submission layout

GeneTrackR development and Bioconductor submission use two different repository trees.

## Development workspace

Keep the full development workspace on a non-submission branch. It may contain:

- `.github/` pkgdown workflows;
- `README.qmd`, `docs/`, and `_pkgdown.yml`;
- `tools/` preflight and README-generation scripts;
- `dev/` development notes;
- local IDE project files that are ignored by Git.

## Package-only submission branch

At formal submission, the GitHub default branch must contain only the package tree used by Bioconductor. Prepare it from the development checkout with:

```bash
Rscript tools/prepare_bioc_submission.001.R . ../GeneTrackR-bioc-submission
```

The package-only tree contains `DESCRIPTION`, `NAMESPACE`, `README.md`, `R/`, `man/`, `inst/`, `tests/`, and `vignettes/`, plus one top-level `.gitignore`. GeneTrackR uses the standard `GPL-3` identifier in `DESCRIPTION` and does not ship separate `NEWS.md` or `LICENSE` files.

Do not copy `.github/`, `dev/`, `docs/`, `tools/`, `README.qmd`, `_pkgdown.yml`, `.Rbuildignore`, `pkgdown-site/`, or an RStudio project file into the submission branch.

## Suggested branch transition

Before opening the Bioconductor submission issue:

1. preserve the current full workspace on a development branch;
2. create a clean package-only branch from the prepared submission tree;
3. make the package-only branch the GitHub default branch used by the submission issue;
4. run `R CMD build`, `R CMD check`, `BiocCheck::BiocCheckGitClone()`, and `BiocCheck::BiocCheck(..., new-package = TRUE)` against that clean Git clone;
5. change the package version to `0.99.0` only for the formal submission candidate.

No GitHub branch changes are performed by the preparation script.
