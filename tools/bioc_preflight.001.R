#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-08-30
# Version: dev001
# Function: Static Bioconductor-readiness preflight for the GeneTrackR source tree
# Input: GeneTrackR package source directory
# Output: PASS/WARN/FAIL readiness report and process exit status

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) > 0L && !startsWith(args[[1L]], "--")) args[[1L]] else "."
root <- normalizePath(root, winslash = "/", mustWork = TRUE)

failures <- character()
warnings <- character()
passes <- character()

add_fail <- function(x) failures <<- c(failures, x)
add_warn <- function(x) warnings <<- c(warnings, x)
add_pass <- function(x) passes <<- c(passes, x)

path <- function(...) file.path(root, ...)

description_file <- path("DESCRIPTION")
if (!file.exists(description_file)) {
  stop("DESCRIPTION was not found under: ", root)
}

desc <- read.dcf(description_file)[1L, , drop = TRUE]
required_fields <- c(
  "Package", "Type", "Title", "Version", "Authors@R", "Description",
  "License", "URL", "BugReports", "Encoding", "Depends", "Imports",
  "Suggests", "VignetteBuilder", "biocViews", "NeedsCompilation"
)
missing_fields <- setdiff(required_fields, names(desc))
if (length(missing_fields) > 0L) {
  add_fail(paste0("DESCRIPTION is missing fields: ", paste(missing_fields, collapse = ", ")))
} else {
  add_pass("DESCRIPTION contains the 0.8.x preflight metadata fields.")
}

if (identical(unname(desc[["NeedsCompilation"]]), "no")) {
  add_pass("NeedsCompilation is no.")
} else {
  add_fail("NeedsCompilation must remain 'no' for the pure-R package.")
}

forbidden_fields <- c("Remotes", "LinkingTo")
present_forbidden <- intersect(forbidden_fields, names(desc))
if (length(present_forbidden) > 0L) {
  add_fail(paste0("Remove unsupported or compiled dependency fields: ", paste(present_forbidden, collapse = ", ")))
} else {
  add_pass("DESCRIPTION has no Remotes or LinkingTo field.")
}

imports <- if ("Imports" %in% names(desc)) trimws(strsplit(desc[["Imports"]], ",", fixed = TRUE)[[1L]]) else character()
if ("Rcpp" %in% imports) {
  add_fail("Rcpp remains in Imports.")
} else {
  add_pass("Rcpp is absent from Imports.")
}

views <- if ("biocViews" %in% names(desc)) trimws(strsplit(desc[["biocViews"]], ",", fixed = TRUE)[[1L]]) else character()
if (length(setdiff(views, "Software")) < 2L) {
  add_fail("biocViews should include at least two software vocabulary terms in addition to Software.")
} else {
  add_pass(paste0("biocViews: ", paste(views, collapse = ", ")))
}

sentence_parts <- strsplit(desc[["Description"]], "[.!?]+")[[1L]]
sentence_parts <- sentence_parts[nzchar(trimws(sentence_parts))]
if (length(sentence_parts) < 3L) {
  add_fail("DESCRIPTION Description must contain at least three complete sentences.")
} else {
  add_pass("DESCRIPTION has at least three complete sentences.")
}

version <- package_version(desc[["Version"]])
if (version < package_version("0.99.0")) {
  add_warn(paste0(
    "Development version is ", desc[["Version"]],
    "; change to 0.99.0 only when creating the formal Bioconductor submission candidate."
  ))
}

compiled_paths <- list.files(
  root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE,
  pattern = "\\.(o|obj|so|dll|a|lib)$",
  ignore.case = TRUE
)
if (dir.exists(path("src")) || file.exists(path("R", "bw_cpp_backend.R")) || length(compiled_paths) > 0L) {
  add_fail("Compiled backend files or build artifacts remain in the source tree.")
} else {
  add_pass("No src directory, compiled BigWig wrapper, or compiled build artifact was found.")
}

r_files <- list.files(path("R"), pattern = "\\.[Rr]$", full.names = TRUE)
r_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
if (any(grepl("Rcpp::|\\.Call\\s*\\(", r_text, perl = TRUE))) {
  add_fail("R production sources still reference Rcpp or .Call().")
} else {
  add_pass("R production sources contain no Rcpp:: or .Call() references.")
}

vignette_files <- list.files(path("vignettes"), pattern = "\\.Rmd$", full.names = TRUE)
bad_vignettes <- character()
for (file in vignette_files) {
  x <- readLines(file, warn = FALSE, n = 40L)
  ok <- length(x) > 0L && trimws(x[[1L]]) == "---" &&
    any(grepl("%\\\\VignetteIndexEntry\\{", x)) &&
    any(grepl("%\\\\VignetteEngine\\{knitr::rmarkdown\\}", x)) &&
    any(grepl("%\\\\VignetteEncoding\\{UTF-8\\}", x))
  if (!ok) bad_vignettes <- c(bad_vignettes, basename(file))
}
if (length(vignette_files) == 0L) {
  add_fail("No Rmd vignettes were found.")
} else if (length(bad_vignettes) > 0L) {
  add_fail(paste0("Formal vignette metadata is missing from: ", paste(bad_vignettes, collapse = ", ")))
} else {
  add_pass(paste0(length(vignette_files), " Rmd vignettes contain formal knitr vignette metadata."))
}

readme_qmd <- path("README.qmd")
if (file.exists(readme_qmd)) {
  readme <- paste(readLines(readme_qmd, warn = FALSE), collapse = "\n")
  required_links <- c(
    "https://renscq.github.io/GeneTrackR/",
    "https://github.com/Renscq/GeneTrackR/issues",
    "#citation",
    "#license"
  )
  missing_links <- required_links[!vapply(required_links, grepl, logical(1), x = readme, fixed = TRUE)]
  if (length(missing_links) > 0L) {
    add_fail(paste0("README navigation is missing: ", paste(missing_links, collapse = ", ")))
  } else {
    add_pass("README navigation exposes documentation, issues, citation, and license links.")
  }
}

if (!grepl("export\\(as_granges\\)", paste(readLines(path("NAMESPACE"), warn = FALSE), collapse = "\n"))) {
  add_fail("NAMESPACE does not export as_granges(), weakening GenomicRanges interoperability.")
} else {
  add_pass("as_granges() is exported for GenomicRanges interoperability.")
}

submission_only_dirs <- c(".github", "dev", "docs", "tools")
present_submission_dirs <- submission_only_dirs[dir.exists(file.path(root, submission_only_dirs))]
if (length(present_submission_dirs) > 0L) {
  add_warn(paste0(
    "Default-branch submission hygiene remains for a later stage: move non-package development directories before formal submission: ",
    paste(present_submission_dirs, collapse = ", "), "."
  ))
}

add_warn(
  "Before submission, document why GeneTrackR keeps native bedGraph/WIG/BigWig I/O instead of delegating common genomic-file import to rtracklayer."
)
add_warn(
  "Run R CMD build, R CMD check, BiocCheck::BiocCheckGitClone(), and BiocCheck::BiocCheck(new-package = TRUE) in the Bioconductor devel environment."
)

cat("GeneTrackR Bioconductor preflight\n")
cat("Source: ", root, "\n\n", sep = "")
for (x in passes) cat("[PASS] ", x, "\n", sep = "")
for (x in warnings) cat("[WARN] ", x, "\n", sep = "")
for (x in failures) cat("[FAIL] ", x, "\n", sep = "")
cat("\nSummary: ", length(passes), " PASS, ", length(warnings), " WARN, ", length(failures), " FAIL\n", sep = "")

quit(status = if (length(failures) > 0L) 1L else 0L)
