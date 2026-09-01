#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-09-01
# Version: dev006
# Function: Prepare a package-only GeneTrackR tree for Bioconductor submission
# Input: Development source directory and optional output directory
# Output: Package-only source tree with development files excluded

args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args) >= 1L) args[[1L]] else "."
output_root <- if (length(args) >= 2L) args[[2L]] else "GeneTrackR-bioc-submission-0.99.1"

source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

required_top_level <- c(
  "DESCRIPTION",
  "NAMESPACE",
  "NEWS.md",
  "README.md",
  "R",
  "man",
  "inst",
  "tests",
  "vignettes"
)

missing <- required_top_level[!file.exists(file.path(source_root, required_top_level))]
if (length(missing) > 0L) {
  stop("Source tree is missing required package paths: ", paste(missing, collapse = ", "))
}

desc <- read.dcf(file.path(source_root, "DESCRIPTION"))[1L, , drop = TRUE]
if (!grepl("^0\\.99\\.[0-9]+$", desc[["Version"]])) {
  stop(
    "Formal Bioconductor submission tree requires a 0.99.x package version; found ",
    desc[["Version"]]
  )
}

if (dir.exists(output_root) || file.exists(output_root)) {
  unlink(output_root, recursive = TRUE, force = TRUE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

copy_path <- function(name) {
  src <- file.path(source_root, name)
  if (!file.exists(src)) {
    return(invisible(FALSE))
  }
  ok <- file.copy(
    src,
    output_root,
    recursive = dir.exists(src),
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(ok)) {
    stop("Failed to copy submission path: ", name)
  }
  invisible(TRUE)
}

submission_paths <- c(
  ".gitignore",
  "DESCRIPTION",
  "NAMESPACE",
  "NEWS.md",
  "README.md",
  "R",
  "man",
  "inst",
  "tests",
  "vignettes"
)
invisible(vapply(submission_paths, copy_path, logical(1L)))

readme_file <- file.path(output_root, "README.md")
if (file.exists(readme_file)) {
  readme_lines <- readLines(readme_file, warn = FALSE, encoding = "UTF-8")
  readme_lines <- readme_lines[!grepl(
    "README.md is assembled from README.qmd",
    readme_lines,
    fixed = TRUE
  )]
  while (length(readme_lines) > 0L && !nzchar(trimws(readme_lines[[1L]]))) {
    readme_lines <- readme_lines[-1L]
  }
  writeLines(readme_lines, readme_file, useBytes = TRUE)
}

# The submission clone must not contain IDE, workflow, pkgdown-source, or
# development-tool files. Keep one top-level .gitignore with local-only files.
gitignore <- c(
  ".Rproj.user",
  ".Rhistory",
  ".RData",
  ".Ruserdata",
  "*.Rproj"
)
writeLines(gitignore, file.path(output_root, ".gitignore"), useBytes = TRUE)

forbidden <- c(
  ".github",
  "dev",
  "docs",
  "tools",
  "README.qmd",
  "_pkgdown.yml",
  "GeneTrackR.Rproj",
  ".Rbuildignore",
  "pkgdown-site"
)
present_forbidden <- forbidden[file.exists(file.path(output_root, forbidden))]
if (length(present_forbidden) > 0L) {
  stop(
    "Submission tree still contains development-only paths: ",
    paste(present_forbidden, collapse = ", ")
  )
}

message("Prepared package-only Bioconductor source tree: ", output_root)
message("Next run R CMD build/check and BiocCheck on a Git clone of this tree.")
