#!/usr/bin/env Rscript

# Generate README.md from README.qmd and docs/*.qmd without executing R code.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
source_file <- file.path(root, "README.qmd")
output_file <- file.path(root, "README.md")

read_utf8 <- function(path) {
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

strip_yaml_front_matter <- function(lines) {
  if (length(lines) == 0L || trimws(lines[[1L]]) != "---") {
    return(lines)
  }

  end_idx <- which(trimws(lines[-1L]) == "---")
  if (length(end_idx) == 0L) {
    stop("README.qmd has an unterminated YAML front matter block.")
  }

  end_idx <- end_idx[[1L]] + 1L
  if (end_idx >= length(lines)) {
    return(character())
  }

  lines[(end_idx + 1L):length(lines)]
}

expand_includes <- function(lines, base_dir, stack = character()) {
  include_pattern <- "^\\s*\\{\\{<\\s*include\\s+([^>]+?)\\s*>\\}\\}\\s*$"
  out <- character()

  for (line in lines) {
    match <- regexec(include_pattern, line, perl = TRUE)
    fields <- regmatches(line, match)[[1L]]

    if (length(fields) == 0L) {
      out <- c(out, line)
      next
    }

    include_path <- trimws(fields[[2L]])
    include_path <- sub("^[\"']", "", include_path)
    include_path <- sub("[\"']$", "", include_path)
    include_file <- normalizePath(
      file.path(base_dir, include_path),
      winslash = "/",
      mustWork = TRUE
    )

    if (include_file %in% stack) {
      stop("Recursive README include detected: ", include_file)
    }

    include_lines <- read_utf8(include_file)
    include_lines <- expand_includes(
      include_lines,
      dirname(include_file),
      c(stack, include_file)
    )
    out <- c(out, include_lines)
  }

  out
}

normalize_code_fences <- function(lines) {
  # README is display-only: keep R examples as syntax-highlighted Markdown code.
  sub("^```\\{r\\}[[:space:]]*$", "```r", lines)
}

render_readme <- function() {
  lines <- read_utf8(source_file)
  lines <- strip_yaml_front_matter(lines)
  lines <- expand_includes(lines, root, normalizePath(source_file, winslash = "/"))
  lines <- normalize_code_fences(lines)
  lines
}

rendered <- render_readme()

if (check_only) {
  if (!file.exists(output_file)) {
    message("README.md is missing. Run: Rscript tools/render_readme.R")
    quit(status = 1L)
  }

  current <- read_utf8(output_file)
  if (!identical(current, rendered)) {
    message("README.md is out of sync with README.qmd/docs/*.qmd.")
    message("Run: Rscript tools/render_readme.R")
    quit(status = 1L)
  }

  message("README.md is synchronized with README.qmd/docs/*.qmd.")
  quit(status = 0L)
}

writeLines(rendered, output_file, useBytes = TRUE)
message("Rendered README.md from README.qmd/docs/*.qmd without executing code chunks.")
