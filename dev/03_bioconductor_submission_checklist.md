# GeneTrackR 0.99.0 submission checklist

## Before submission

- [ ] Use R 4.6.0 and Bioconductor 3.24 devel.
- [ ] `BiocManager::valid()` reports a consistent devel installation.
- [ ] Bioconductor Support Site account is registered with the GeneTrackR maintainer identity/email.
- [ ] The submitter is the maintainer listed in `DESCRIPTION`.
- [ ] The submitting GitHub account has an SSH public key.
- [ ] Subscribe to `bioc-devel` (recommended for package maintainers).
- [ ] Preserve the full GeneTrackR workspace on a non-default development branch.
- [ ] Make the package-only GeneTrackR 0.99.0 branch the GitHub default branch.
- [ ] Confirm the GitHub default branch contains no `.github/`, `dev/`, `docs/`, `tools/`, `README.qmd`, `_pkgdown.yml`, `.Rbuildignore`, `pkgdown-site/`, or IDE project files.

## Safe GitHub branch preparation

A low-risk way to avoid rewriting the existing development branch is to create an independent package-only branch from the prepared submission tree:

```bash
cd ../GeneTrackR-bioc-submission
git init
git switch -c bioc-submission
git remote add origin git@github.com:Renscq/GeneTrackR.git
git add .
git commit -m "chore(bioc): prepare 0.99.0 submission"
git push -u origin bioc-submission
```

Then, in GitHub repository settings, change the repository **default branch** to `bioc-submission` before opening the Bioconductor issue. Do not delete the existing full development branch. The submission service uses only the current default branch.

## Local validation

- [ ] `devtools::document()` completes without skipped/nameless topics.
- [ ] `devtools::test()` passes.
- [ ] `devtools::check()` reports 0 errors and 0 warnings; resolve avoidable notes.
- [ ] `BiocCheck::BiocCheckGitClone()` reports 0 errors and 0 warnings.
- [ ] `R CMD build .` produces `GeneTrackR_0.99.0.tar.gz`.
- [ ] `R CMD check --as-cran GeneTrackR_0.99.0.tar.gz` completes successfully.
- [ ] `BiocCheck::BiocCheck("GeneTrackR_0.99.0.tar.gz", new-package = TRUE)` reports 0 errors and 0 warnings.
- [ ] Network-dependent repository checks are rerun if CRAN/Bioconductor URLs timed out.

## Expected advisory notes

These do not justify large refactors unless the reviewer requests them:

- optional ORCID metadata if no ORCID is supplied;
- optional `fnd` role when no grant/funder should be declared;
- intentionally scoped `suppressWarnings()` calls;
- long-function and formatting recommendations;
- optional `CITATION` when no associated package publication/DOI exists.

## Submission

- [ ] Open a new issue in `Bioconductor/BiocContributions`.
- [ ] Use `GeneTrackR` as the issue title.
- [ ] Provide the GitHub repository URL.
- [ ] Confirm all public-review, maintainer-email, naming-policy, code-of-conduct, and devel/release maintenance checkboxes.
- [ ] Monitor the issue for moderation and build reports.

## After validation / policy acceptance

- [ ] If validation succeeds, reply to the issue exactly `/accept-policies`.
- [ ] Follow the issue instructions to link/push to the Bioconductor staging GitHub repository created for the submission.
- [ ] Trigger review builds from the staging repository, not from the original personal GitHub repository.
- [ ] Increment the patch version for each build-triggering revision: `0.99.1`, `0.99.2`, ...
- [ ] Reply to reviewer comments in the submission issue with a concise point-by-point summary of changes.

## After acceptance

- [ ] Follow the acceptance instructions to add the canonical `git.bioconductor.org` remote.
- [ ] Configure the BiocCredentials/SSH access created during acceptance.
- [ ] Maintain future devel/release changes through the canonical Bioconductor Git repository and keep GitHub synchronized as desired.
