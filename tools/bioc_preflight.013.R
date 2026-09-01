#!/usr/bin/env Rscript

# Author: Shuchao Ren
# Date: 2026-09-01
# Version: dev013
# Function: Formal 0.99.2 Bioconductor submission preflight for the GeneTrackR source tree
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

read_text_file <- function(file) {
  if (!file.exists(file)) return("")
  paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

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
  add_pass("DESCRIPTION contains the formal submission metadata fields.")
}

if (identical(unname(desc[["NeedsCompilation"]]), "no")) {
  add_pass("NeedsCompilation is no.")
} else {
  add_fail("NeedsCompilation must remain 'no' for the pure-R package.")
}

forbidden_fields <- c("Remotes", "Additional_repositories", "LinkingTo")
present_forbidden <- intersect(forbidden_fields, names(desc))
if (length(present_forbidden) > 0L) {
  add_fail(paste0("Remove unsupported or compiled dependency fields: ", paste(present_forbidden, collapse = ", ")))
} else {
  add_pass("DESCRIPTION has no Remotes, Additional_repositories, or LinkingTo field.")
}

split_dependencies <- function(field) {
  if (!field %in% names(desc)) return(character())
  out <- trimws(strsplit(desc[[field]], ",", fixed = TRUE)[[1L]])
  sub("\\s*\\(.*$", "", out)
}

imports <- split_dependencies("Imports")
suggests <- split_dependencies("Suggests")
depends <- split_dependencies("Depends")

if ("Rcpp" %in% imports) {
  add_fail("Rcpp remains in Imports.")
} else {
  add_pass("Rcpp is absent from Imports.")
}

duplicate_dependencies <- unique(c(
  intersect(imports, suggests),
  intersect(imports, depends),
  intersect(suggests, depends)
))
if (length(duplicate_dependencies) > 0L) {
  add_fail(paste0(
    "Packages appear in more than one dependency field: ",
    paste(duplicate_dependencies, collapse = ", ")
  ))
} else {
  add_pass("Depends, Imports, and Suggests contain no duplicate package entries.")
}

required_imports <- c(
  "data.table", "ggplot2", "patchwork", "rlang", "RColorBrewer",
  "GenomicRanges", "IRanges", "S4Vectors", "grDevices", "grid",
  "stats", "tools", "utils"
)
missing_imports <- setdiff(required_imports, imports)
if (length(missing_imports) > 0L) {
  add_fail(paste0(
    "Runtime namespace dependencies are missing from Imports: ",
    paste(missing_imports, collapse = ", ")
  ))
} else {
  add_pass("Runtime namespace dependencies are declared in Imports.")
}

if ("Rsamtools" %in% suggests && !"Rsamtools" %in% imports) {
  add_pass("Rsamtools remains an optional Suggests dependency for indexed query acceleration.")
} else {
  add_fail("Rsamtools should remain in Suggests, not Imports, while tabix acceleration is optional.")
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

version <- unname(desc[["Version"]])
if (!grepl("^0\\.99\\.[0-9]+$", version)) {
  add_fail(paste0(
    "Formal Bioconductor submission version must be 0.99.x; found ",
    version, "."
  ))
} else {
  add_pass(paste0("Formal Bioconductor pre-release version: ", version))
}

if (!identical(version, "0.99.2")) {
  add_fail(paste0("Current GeneTrackR submission candidate must be 0.99.2; found ", version, "."))
} else {
  add_pass("Current submission candidate version is 0.99.2.")
}

if (!identical(unname(desc[["Title"]]), "Genomic Track Visualization and Haplotype Analysis")) {
  add_fail("DESCRIPTION Title is not the canonical GeneTrackR title.")
} else {
  add_pass("DESCRIPTION Title matches the canonical GeneTrackR title.")
}

authors_r <- unname(desc[["Authors@R"]])
if (!grepl('given = "Shuchao"', authors_r, fixed = TRUE) ||
    !grepl('family = "Ren"', authors_r, fixed = TRUE) ||
    !grepl('email = "rensc0718@163.com"', authors_r, fixed = TRUE)) {
  add_fail("DESCRIPTION Authors@R must identify Shuchao Ren as author/maintainer with the canonical maintainer email.")
} else {
  add_pass("DESCRIPTION Authors@R identifies Shuchao Ren as author/maintainer.")
}

if (!identical(unname(desc[["Date"]]), "2026-09-01")) {
  add_fail("DESCRIPTION Date must be 2026-09-01 for the 0.99.2 candidate.")
} else {
  add_pass("DESCRIPTION Date matches the 0.99.2 candidate date.")
}

if (file.exists(path("inst", "CITATION"))) {
  add_warn("inst/CITATION is present; keep it only when GeneTrackR has an associated publication/DOI. Package-only citation is generated from DESCRIPTION.")
} else {
  add_pass("Package citation is generated from DESCRIPTION; no custom inst/CITATION is present.")
}

gitattributes_file <- path(".gitattributes")
if (file.exists(gitattributes_file)) {
  gitattributes_text <- read_text_file(gitattributes_file)
  if (grepl("filter=lfs", gitattributes_text, fixed = TRUE)) {
    add_fail("Git LFS is not allowed for Bioconductor new-package submissions.")
  } else {
    add_pass(".gitattributes does not configure Git LFS.")
  }
} else {
  add_pass("No .gitattributes/Git LFS configuration is present.")
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

namespace_text <- paste(readLines(path("NAMESPACE"), warn = FALSE), collapse = "\n")
if (!grepl("export\\(as_granges\\)", namespace_text)) {
  add_fail("NAMESPACE does not export as_granges(), weakening GenomicRanges interoperability.")
} else {
  add_pass("as_granges() is exported for GenomicRanges interoperability.")
}

as_granges_file <- path("R", "as_granges.R")
if (!file.exists(as_granges_file)) {
  add_fail("R/as_granges.R is missing.")
} else {
  as_granges_text <- paste(readLines(as_granges_file, warn = FALSE), collapse = "\n")
  required_classes <- c("GenePred", "Feature", "FeatureTrack", "BwgTrack", "VariantTrack")
  missing_classes <- required_classes[!vapply(
    required_classes,
    function(class_name) grepl(class_name, as_granges_text, fixed = TRUE),
    logical(1L)
  )]
  if (length(missing_classes) > 0L) {
    add_fail(paste0(
      "as_granges() does not expose all core genomic track classes: ",
      paste(missing_classes, collapse = ", ")
    ))
  } else {
    add_pass("Annotation, signal, and variant track objects expose a GRanges conversion path.")
  }
}

retrieve_vcf_text <- paste(readLines(path("R", "retrieve_vcf.R"), warn = FALSE), collapse = "\n")
if (grepl('"GRanges"', retrieve_vcf_text, fixed = TRUE)) {
  add_pass("retrieve_vcf() exposes a direct GRanges return path.")
} else {
  add_fail("retrieve_vcf() does not expose a GRanges return path.")
}

development_only_paths <- c(
  ".github", "dev", "docs", "tools", "README.qmd", "_pkgdown.yml",
  "GeneTrackR.Rproj", "pkgdown-site"
)
present_development_paths <- development_only_paths[
  file.exists(file.path(root, development_only_paths))
]
if (length(present_development_paths) > 0L) {
  add_warn(paste0(
    "This is a development workspace, not a package-only submission clone. ",
    "Prepare the default-branch submission tree without: ",
    paste(present_development_paths, collapse = ", "), "."
  ))
} else {
  add_pass("The source root contains no known development-only submission paths.")
}

if (file.exists(path("GeneTrackR.Rproj"))) {
  add_warn("GeneTrackR.Rproj is suitable for local development only and must not be tracked in the Bioconductor submission clone.")
}

signal_vignette <- path("vignettes", "signal-tracks.Rmd")
if (file.exists(signal_vignette)) {
  signal_text <- paste(readLines(signal_vignette, warn = FALSE), collapse = "\n")
  if (grepl("rtracklayer", signal_text, fixed = TRUE) && grepl("GRanges", signal_text, fixed = TRUE)) {
    add_pass("Signal documentation explains the native-I/O/rtracklayer boundary and GRanges exchange path.")
  } else {
    add_fail("Signal documentation must explain the native-I/O/rtracklayer boundary and GRanges exchange path.")
  }
}

# Bioconductor documentation/runtime checks added for 0.8.2.
example_sources <- c(r_files, list.files(path("man"), pattern = "\\.Rd$", full.names = TRUE))
dontrun_hits <- character()
for (file in example_sources) {
  x <- readLines(file, warn = FALSE)
  if (any(grepl("\\\\dontrun\\{", x))) {
    dontrun_hits <- c(dontrun_hits, sub(paste0("^", root, "/?"), "", file))
  }
}
if (length(dontrun_hits) > 0L) {
  add_fail(paste0(
    "Runnable documentation still contains \\dontrun{} blocks: ",
    paste(unique(dontrun_hits), collapse = ", ")
  ))
} else {
  add_pass("Roxygen and Rd examples contain no \\dontrun{} blocks.")
}

required_vignette_sections <- c("## Introduction", "## Installation", "## Session information")
missing_vignette_sections <- character()
missing_vignette_metadata <- character()
missing_landing_links <- character()
missing_dpi_policy <- character()
vignette_index_numbers <- character()
for (file in vignette_files) {
  x <- readLines(file, warn = FALSE)
  text <- paste(x, collapse = "\n")
  absent <- required_vignette_sections[!vapply(
    required_vignette_sections, grepl, logical(1L), x = text, fixed = TRUE
  )]
  if (length(absent) > 0L) {
    missing_vignette_sections <- c(
      missing_vignette_sections,
      paste0(basename(file), " [", paste(absent, collapse = ", "), "]")
    )
  }
  if (!grepl('author: "Shuchao Ren"', text, fixed = TRUE) ||
      !grepl('date: "2026-08-31"', text, fixed = TRUE)) {
    missing_vignette_metadata <- c(missing_vignette_metadata, basename(file))
  }
  if (!grepl("https://renscq.github.io/GeneTrackR/", text, fixed = TRUE)) {
    missing_landing_links <- c(missing_landing_links, basename(file))
  }
  if (!grepl("GENETRACKR_PKGDOWN", text, fixed = TRUE) ||
      !grepl("dpi = gtr_vignette_dpi", text, fixed = TRUE)) {
    missing_dpi_policy <- c(missing_dpi_policy, basename(file))
  }
  index_line <- grep("%\\VignetteIndexEntry{", x, value = TRUE, fixed = TRUE)
  if (length(index_line) == 1L) {
    entry_text <- strsplit(index_line, "{", fixed = TRUE)[[1L]][[2L]]
    index_number <- substr(entry_text, 1L, 2L)
    if (grepl("^[0-9]{2}$", index_number)) {
      vignette_index_numbers <- c(vignette_index_numbers, index_number)
    }
  }
}
if (length(missing_vignette_sections) > 0L) {
  add_fail(paste0(
    "Bioconductor vignette sections are incomplete: ",
    paste(missing_vignette_sections, collapse = "; ")
  ))
} else {
  add_pass("All formal vignettes include Introduction, Installation, and Session information sections.")
}
if (length(missing_vignette_metadata) > 0L) {
  add_fail(paste0(
    "Vignette author/date metadata is missing from: ",
    paste(missing_vignette_metadata, collapse = ", ")
  ))
} else {
  add_pass("All formal vignettes include author and last-modification metadata.")
}
if (length(missing_landing_links) > 0L) {
  add_fail(paste0(
    "Package landing-page links are missing from: ",
    paste(missing_landing_links, collapse = ", ")
  ))
} else {
  add_pass("All formal vignettes link to the GeneTrackR package website.")
}
if (length(missing_dpi_policy) > 0L) {
  add_fail(paste0(
    "Build-aware vignette DPI policy is missing from: ",
    paste(missing_dpi_policy, collapse = ", ")
  ))
} else {
  add_pass("All formal vignettes use the build-aware 144/300 dpi policy.")
}
expected_vignette_indices <- sprintf("%02d", seq_along(vignette_files))
if (!identical(sort(vignette_index_numbers), expected_vignette_indices)) {
  add_fail("Formal vignette index entries must provide one unique 01-15 display-order prefix.")
} else {
  add_pass("Formal vignette index entries define one unique 01-15 display order.")
}

pkgdown_workflow <- path(".github", "workflows", "pkgdown.yaml")
if (file.exists(pkgdown_workflow)) {
  workflow_text <- paste(readLines(pkgdown_workflow, warn = FALSE), collapse = "\n")
  if (grepl('GENETRACKR_PKGDOWN: "true"', workflow_text, fixed = TRUE)) {
    add_pass("pkgdown explicitly requests 300-dpi documentation rendering.")
  } else {
    add_fail("pkgdown workflow does not set GENETRACKR_PKGDOWN=true.")
  }
}

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
file_sizes <- file.info(all_files)$size
names(file_sizes) <- all_files
large_files <- names(file_sizes)[is.finite(file_sizes) & file_sizes > 5 * 1024^2]
if (length(large_files) > 0L) {
  add_fail(paste0(
    "Individual files exceed the 5 MB Bioconductor software-package guideline: ",
    paste(basename(large_files), collapse = ", ")
  ))
} else {
  add_pass("No individual source-tree file exceeds 5 MB.")
}

# Estimate the R-build payload by excluding development-only paths already covered by .Rbuildignore.
relative_files <- sub(paste0("^", root, "/?"), "", all_files)
excluded_prefixes <- c(".github/", "dev/", "docs/", "tools/", "pkgdown-site/", ".Rproj.user/")
keep <- !vapply(relative_files, function(x) {
  x %in% c("GeneTrackR.Rproj", "README.qmd", "_pkgdown.yml") ||
    any(startsWith(x, excluded_prefixes))
}, logical(1L))
estimated_build_size <- sum(file_sizes[keep], na.rm = TRUE)
if (estimated_build_size > 10 * 1024^2) {
  add_fail(sprintf(
    "Estimated R-build payload is %.2f MB, above the 10 MB Bioconductor guideline.",
    estimated_build_size / 1024^2
  ))
} else {
  add_pass(sprintf(
    "Estimated R-build payload is %.2f MB, below the 10 MB guideline.",
    estimated_build_size / 1024^2
  ))
}

rbuildignore_file <- path(".Rbuildignore")
if (file.exists(rbuildignore_file)) {
  ignore_text <- paste(readLines(rbuildignore_file, warn = FALSE), collapse = "\n")
  required_ignores <- c("GeneTrackR\\.Rproj", "\\.github", "dev", "docs", "tools", "README\\.qmd", "_pkgdown\\.yml")
  missing_ignores <- required_ignores[!vapply(
    required_ignores, grepl, logical(1L), x = ignore_text
  )]
  if (length(missing_ignores) > 0L) {
    add_fail(paste0(
      ".Rbuildignore is missing development-only exclusions: ",
      paste(missing_ignores, collapse = ", ")
    ))
  } else {
    add_pass(".Rbuildignore excludes development-only project, docs, tools, and workflow sources from R builds.")
  }
}

# Current BiocCheck checks require every R vignette chunk to have a label.
unlabeled_chunks <- character()
duplicate_chunk_labels <- character()
for (file in vignette_files) {
  x <- readLines(file, warn = FALSE)
  chunk_lines <- grep("^```\\{r", x, value = TRUE)
  labels <- character()
  for (line in chunk_lines) {
    header <- sub("^```\\{r", "", line)
    header <- sub("\\}.*$", "", header)
    header <- trimws(header)
    if (!nzchar(header) || startsWith(header, ",")) {
      unlabeled_chunks <- c(unlabeled_chunks, basename(file))
      next
    }
    label <- trimws(strsplit(header, ",", fixed = TRUE)[[1L]][[1L]])
    if (!nzchar(label) || grepl("=", label, fixed = TRUE)) {
      unlabeled_chunks <- c(unlabeled_chunks, basename(file))
    } else {
      labels <- c(labels, label)
    }
  }
  dup <- unique(labels[duplicated(labels)])
  if (length(dup) > 0L) {
    duplicate_chunk_labels <- c(
      duplicate_chunk_labels,
      paste0(basename(file), " [", paste(dup, collapse = ", "), "]")
    )
  }
}
if (length(unlabeled_chunks) > 0L) {
  add_fail(paste0(
    "All vignette R chunks need labels; anonymous chunks remain in: ",
    paste(unique(unlabeled_chunks), collapse = ", ")
  ))
} else {
  add_pass("Every R code chunk in every formal vignette has a label.")
}
if (length(duplicate_chunk_labels) > 0L) {
  add_fail(paste0(
    "Duplicate vignette chunk labels remain: ",
    paste(duplicate_chunk_labels, collapse = "; ")
  ))
} else {
  add_pass("Formal vignettes contain no duplicate R chunk labels.")
}

# Flag common BiocCheck R-code style issues without rewriting analysis code.
r_code_only <- unlist(lapply(r_files, function(file) {
  x <- readLines(file, warn = FALSE)
  x[!grepl("^[[:space:]]*#", x)]
}), use.names = FALSE)
if (any(grepl("<<-", r_code_only, fixed = TRUE))) {
  add_warn("R production code still contains <<-; Bioconductor recommends avoiding super-assignment.")
} else {
  add_pass("R production code contains no <<- super-assignment.")
}
if (any(grepl("\\bsapply\\s*\\(", r_code_only, perl = TRUE))) {
  add_warn("R production code contains sapply(); Bioconductor recommends vapply() where possible.")
} else {
  add_pass("R production code contains no sapply() calls.")
}
if (any(grepl("\\b1\\s*:\\s*(length|nrow|ncol)\\s*\\(", r_code_only, perl = TRUE))) {
  add_warn("R production code contains 1:length()/1:nrow()/1:ncol() style iteration.")
} else {
  add_pass("R production code contains no 1:length()/1:nrow()/1:ncol() iteration patterns.")
}

# Formatting and submission-facing metadata checks.
long_line_count <- 0L
for (file in c(r_files, vignette_files)) {
  x <- readLines(file, warn = FALSE)
  long_line_count <- long_line_count + sum(nchar(x, type = "width") > 80L)
}
if (long_line_count > 0L) {
  add_warn(paste0(
    long_line_count,
    " R/vignette source lines exceed 80 characters; BiocCheck may report a formatting NOTE."
  ))
} else {
  add_pass("R and vignette source lines are at most 80 characters wide.")
}

readme_file <- path("README.md")
if (file.exists(readme_file)) {
  readme_text <- paste(readLines(readme_file, warn = FALSE), collapse = "\n")
  if (grepl('BiocManager::install("GeneTrackR")', readme_text, fixed = TRUE)) {
    add_pass("README documents the Bioconductor installation path.")
  } else {
    add_fail("README should document BiocManager::install(\"GeneTrackR\").")
  }
}

readme_file <- path("README.md")
if (file.exists(readme_file)) {
  readme_text <- read_text_file(readme_file)
  citation_ok <-
    grepl("Shuchao Ren (2026)", readme_text, fixed = TRUE) &&
    grepl("R package version 0.99.2", readme_text, fixed = TRUE) &&
    grepl("author = {{Shuchao Ren}}", readme_text, fixed = TRUE)
  if (citation_ok) {
    add_pass("README citation metadata identifies Shuchao Ren and GeneTrackR 0.99.2.")
  } else {
    add_fail("README citation metadata is not synchronized to Shuchao Ren / GeneTrackR 0.99.2.")
  }
}

if (identical(unname(desc[["License"]]), "GPL-3")) {
  add_pass("DESCRIPTION uses the standard GPL-3 license identifier.")
} else {
  add_warn("Use a specific standard license identifier before formal submission.")
}

news_file <- path("NEWS.md")
if (!file.exists(news_file)) {
  add_fail("NEWS.md is required for GeneTrackR Bioconductor release notes.")
} else {
  news_text <- paste(readLines(news_file, warn = FALSE), collapse = "\n")
  if (!grepl("# GeneTrackR 0.99.2", news_text, fixed = TRUE)) {
    add_fail("NEWS.md does not contain the current GeneTrackR 0.99.2 submission entry.")
  } else {
    add_pass("NEWS.md is present and contains the current submission entry.")
  }
}

# Bioconductor requires bundled extdata provenance documentation.
demo_doc <- path("man", "GeneTrackR-demo-data.Rd")
demo_source_doc <- path("inst", "scripts", "README.md")
if (file.exists(demo_doc) && file.exists(demo_source_doc)) {
  add_pass("Bundled extdata files have a package help topic and inst/scripts provenance documentation.")
} else {
  add_fail("Bundled extdata needs both a package help topic and inst/scripts provenance documentation.")
}


# 0.8.8 runnable writer example hotfix checks
r_files <- list.files(path("R"), pattern = "\\.R$", full.names = TRUE)
r_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)

if (any(grepl(".Deprecated(", r_text, fixed = TRUE))) {
  add_fail("BiocCheck .Deprecated usage remains in R source.")
} else {
  add_pass("BiocCheck .Deprecated usage is absent from R source.")
}

if (any(grepl("<<-", r_text, fixed = TRUE))) {
  add_fail("BiocCheck <<- usage remains in R source.")
} else {
  add_pass("BiocCheck <<- usage is absent from R source.")
}

if (grepl("R (>= 4.6.0)", desc[["Depends"]], fixed = TRUE)) {
  add_pass("DESCRIPTION uses the Bioconductor R dependency floor R >= 4.6.0.")
} else {
  add_fail("DESCRIPTION should depend on R >= 4.6.0 for the Bioconductor 3.24 devel baseline.")
}

if (grepl("VariantAnnotation", desc[["biocViews"]], fixed = TRUE)) {
  add_pass("DESCRIPTION includes the suggested VariantAnnotation biocView.")
} else {
  add_fail("DESCRIPTION is missing the suggested VariantAnnotation biocView.")
}

# 0.8.8 runnable writer example hotfix checks
class_feature_file <- path("R", "class_feature_track.R")
if (file.exists(class_feature_file)) {
  class_feature_text <- paste(readLines(class_feature_file, warn = FALSE), collapse = "\n")
  method_count <- lengths(regmatches(
    class_feature_text,
    gregexpr("#' @rdname as_feature", class_feature_text, fixed = TRUE)
  ))
  if (method_count >= 4L) {
    add_pass("All as_feature() S3 methods share the documented as_feature topic.")
  } else {
    add_fail("All as_feature() S3 methods must use @rdname as_feature.")
  }
}

merge_test_file <- path("tests", "testthat", "test-merge-feature.R")
if (file.exists(merge_test_file)) {
  merge_test_text <- paste(readLines(merge_test_file, warn = FALSE), collapse = "\n")
  if (grepl("expect_warning(\n    expect_error(\n      merge_feature(first, second, conflict = \"error\")", merge_test_text, fixed = TRUE)) {
    add_fail("The merge_feature conflict='error' regression test still expects a preliminary warning.")
  } else {
    add_pass("The merge_feature conflict='error' regression test expects the direct error path.")
  }
}


# 0.8.8 runnable writer example hotfix checks
write_feature_source <- read_text_file(path("R", "write_feature.R"))
if (grepl("gtr_demo_features.bed", write_feature_source, fixed = TRUE) &&
    grepl("write_feature_track(features, outfile)", write_feature_source, fixed = TRUE)) {
  add_fail("write_feature_track() example still uses generic BED features with the default BED12 writer.")
} else if (all(vapply(
  c("gtr_demo.genePredExt", "features <- as_feature(gp)", "write_feature_track(features, outfile)"),
  function(token) grepl(token, write_feature_source, fixed = TRUE),
  logical(1L)
))) {
  add_pass("write_feature_track() example uses a gene-model Feature compatible with default BED12 output.")
}


# 0.8.9 final BiocCheck source cleanup checks
advanced_source <- read_text_file(path("R", "advanced_parameters.R"))
demo_source <- read_text_file(path("R", "demo_data.R"))
class_feature_source <- read_text_file(path("R", "class_feature_track.R"))
write_bwg_source <- read_text_file(path("R", "write_bwg.R"))
merge_feature_source <- read_text_file(path("R", "merge_feature.R"))

if (grepl("@return No value is returned; this is a documentation-only help topic.", advanced_source, fixed = TRUE) &&
    grepl("@return No value is returned; this is a documentation-only help topic.", demo_source, fixed = TRUE)) {
  add_pass("Documentation-only help topics declare explicit return values for BiocCheck.")
} else {
  add_fail("Documentation-only help topics must declare explicit return values.")
}

if (all(vapply(
  c("#' @examples", "GenePred(gp$transcripts, gp$exons, gp$genes)"),
  function(token) grepl(token, class_feature_source, fixed = TRUE),
  logical(1L)
))) {
  add_pass("GenePred() runnable example is preserved in roxygen source.")
} else {
  add_fail("GenePred() needs a runnable example in the roxygen source block.")
}

if (grepl('stop(paste0("File exists: ", file)', write_bwg_source, fixed = TRUE)) {
  add_fail("write_bwg.R still uses paste inside the file-exists condition signal.")
} else if (grepl('stop(sprintf("File exists: %s", file)', write_bwg_source, fixed = TRUE)) {
  add_pass("write_bwg.R avoids paste inside the file-exists condition signal.")
}

if (grepl('cannot be merged with `conflict = \"error\"`', merge_feature_source, fixed = TRUE)) {
  add_fail("merge_feature strict conflict message can trigger BiocCheck signal-in-signaler detection.")
} else if (grepl("cannot be merged under strict conflict handling", merge_feature_source, fixed = TRUE)) {
  add_pass("merge_feature strict conflict message avoids BiocCheck signal-in-signaler false positives.")
}

current_r_version <- package_version(as.character(getRversion()))
if (current_r_version < package_version("4.6.0")) {
  add_fail(paste0(
    "Formal submission checks require R >= 4.6.0; current R is ",
    current_r_version, "."
  ))
} else {
  add_pass(paste0("Submission preflight is running under R ", current_r_version, "."))
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
  bioc_version <- as.character(BiocManager::version())
  if (identical(bioc_version, "3.24")) {
    add_pass("Bioconductor devel version is 3.24.")
  } else {
    add_fail(paste0(
      "Formal submission checks target Bioconductor 3.24 devel; current version is ",
      bioc_version, "."
    ))
  }
} else {
  add_fail("BiocManager is required to verify the Bioconductor 3.24 devel environment.")
}

add_warn(
  "Confirm the maintainer is registered on the Bioconductor Support Site before submission."
)
add_warn(
  "Confirm the submitting GitHub account has an SSH public key and is the DESCRIPTION maintainer."
)
add_warn(
  "Before formal review, retain the documented native-I/O/rtracklayer interoperability rationale for reviewer discussion."
)
add_warn(
  "Run R CMD build/check, BiocCheckGitClone(), and BiocCheck(new-package = TRUE) on the package-only Git default branch."
)

cat("GeneTrackR Bioconductor preflight\n")
cat("Source: ", root, "\n\n", sep = "")
for (x in passes) cat("[PASS] ", x, "\n", sep = "")
for (x in warnings) cat("[WARN] ", x, "\n", sep = "")
for (x in failures) cat("[FAIL] ", x, "\n", sep = "")
cat("\nSummary: ", length(passes), " PASS, ", length(warnings), " WARN, ", length(failures), " FAIL\n", sep = "")

quit(status = if (length(failures) > 0L) 1L else 0L)
