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
Rscript tools/bioc_preflight.008.R .
Rscript tools/prepare_bioc_submission.004.R . ../GeneTrackR-bioc-submission
```

The package-only tree contains `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `README.md`, `R/`, `man/`, `inst/`, `tests/`, and `vignettes/`, plus one top-level `.gitignore`. GeneTrackR uses the standard `GPL-3` identifier in `DESCRIPTION`, retains `NEWS.md` as the package release history, and does not ship a separate `LICENSE` file.

Do not copy `.github/`, `dev/`, `docs/`, `tools/`, `README.qmd`, `_pkgdown.yml`, `.Rbuildignore`, `pkgdown-site/`, or an RStudio project file into the submission branch.

## Suggested branch transition

Before opening the Bioconductor submission issue:

1. preserve the current full workspace on a development branch;
2. create a clean package-only branch from the prepared submission tree;
3. make the package-only branch the GitHub default branch used by the submission issue;
4. run `R CMD build`, `R CMD check`, `BiocCheck::BiocCheckGitClone()`, and `BiocCheck::BiocCheck(..., new-package = TRUE)` against that clean Git clone;
5. change the package version to `0.99.0` only for the formal submission candidate.

No GitHub branch changes are performed by the preparation script.

## Windows Git ownership troubleshooting

`BiocCheckGitClone()` uses `gert::git_ls()` to inspect tracked files when the
source directory is a Git clone. On Windows, libgit2 can reject a repository
whose filesystem owner differs from the user running R, even when the repository
is the developer's own project.

First reproduce the Git-layer check directly:

```r
gert::git_ls(".")
```

For a trusted local checkout, add only the GeneTrackR repository to Git's global
safe-directory list:

```bash
git config --global --add safe.directory "E:/rensc/programme/Rscripts/GeneTrackR"
```

If the libgit2 error explicitly names the `.git` directory and the repository
root entry alone is insufficient, also trust that exact path:

```bash
git config --global --add safe.directory "E:/rensc/programme/Rscripts/GeneTrackR/.git"
```

Avoid `safe.directory=*`; trusting one explicit project preserves Git's ownership
protection for other repositories. Re-run `gert::git_ls(".")` before re-running
`BiocCheck::BiocCheckGitClone()`.
