#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-11-19
# Version: dev001
# Function: Run complete source, archive and new-package validation for CI.
# Input: A package source directory and a fresh report directory.
# Output: Source archive, R CMD check, BiocCheck objects and severity summaries.

run_r_command <- function(arguments, logfile) {
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
  status <- system2(executable, arguments, stdout = logfile, stderr = logfile)
  if (status != 0L) {
    stop(sprintf("R command exited %d; see %s.", status, logfile), call. = FALSE)
  }
  invisible(status)
}

checker_counts <- function(result, stage) {
  data.frame(stage = stage, ERROR = length(result$error),
               WARNING = length(result$warning), NOTE = length(result$note))
}

package_check_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- paste(
    "Usage: Rscript tools/check_package.R [options]",
    "  --path DIR          Package source (default: current directory)",
    "  --output-dir DIR    Fresh reports/archive directory (default: temp/package-check)",
    "  --help              Show this help",
    "Runs full build/check, BiocCheckGitClone and new-package BiocCheck.",
    "No checks are skipped. ERROR or WARNING results fail the run.", sep = "\n")
  if (any(args %in% c("--help", "-h"))) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  values <- list(path = ".", output_dir = "temp/package-check")
  mapping <- c("--path" = "path", "--output-dir" = "output_dir")
  index <- 1L
  while (index <= length(args)) {
    flag <- args[[index]]
    if (!flag %in% names(mapping) || index == length(args)) {
      stop(paste("Invalid or incomplete option:", flag, "\n", help), call. = FALSE)
    }
    values[[mapping[[flag]]]] <- args[[index + 1L]]
    index <- index + 2L
  }
  if (!requireNamespace("BiocCheck", quietly = TRUE)) {
    stop("Install the developer dependency BiocCheck to run package checks.", call. = FALSE)
  }
  source_path <- normalizePath(values$path, winslash = "/", mustWork = TRUE)
  desc <- read.dcf(file.path(source_path, "DESCRIPTION"))
  package <- unname(desc[1L, "Package"])
  version <- unname(desc[1L, "Version"])
  if (!grepl("^[A-Za-z][A-Za-z0-9.]*$", package) || !grepl("^[0-9.]+$", version)) {
    stop("Package name or version cannot be used as a source archive name.", call. = FALSE)
  }
  dir.create(values$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(values$output_dir, winslash = "/", mustWork = TRUE)
  archive <- file.path(output_dir, paste0(package, "_", version, ".tar.gz"))
  if (file.exists(archive) || dir.exists(file.path(output_dir, paste0(package, ".Rcheck")))) {
    stop("Choose a fresh output directory to retain earlier validation results.", call. = FALSE)
  }
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(output_dir)
  cat("Checking source repository layout.\n")
  clone <- BiocCheck::BiocCheckGitClone(source_path, `quit-with-status` = FALSE)
  saveRDS(clone, file.path(output_dir, "bioccheck-gitclone.rds"))
  clone_counts <- checker_counts(clone, "BiocCheckGitClone")
  utils::write.csv(clone_counts, file.path(output_dir, "bioccheck-summary.csv"), row.names = FALSE)
  if (clone_counts$ERROR > 0L || clone_counts$WARNING > 0L) {
    stop("Source repository checks reported ERROR or WARNING results.", call. = FALSE)
  }
  cat("Building all package vignettes.\n")
  run_r_command(c("CMD", "build", shQuote(source_path)), file.path(output_dir, "build.log"))
  if (!file.exists(archive)) stop("The expected source archive was not built.", call. = FALSE)
  cat("Running full R CMD check, including tests, vignettes and manual.\n")
  run_r_command(c("CMD", "check", shQuote(archive)), file.path(output_dir, "r-cmd-check.log"))
  check_log <- file.path(output_dir, paste0(package, ".Rcheck"), "00check.log")
  if (!file.exists(check_log) || !any(grepl("^Status: OK", readLines(check_log, warn = FALSE)))) {
    stop("R CMD check did not report Status: OK; inspect the check log.", call. = FALSE)
  }
  cat("Running full new-package BiocCheck.\n")
  result <- BiocCheck::BiocCheck(archive, checkDir = output_dir,
    install = TRUE, `new-package` = TRUE, `quit-with-status` = FALSE)
  saveRDS(result, file.path(output_dir, "bioccheck-new-package.rds"))
  counts <- rbind(clone_counts, checker_counts(result, "BiocCheck-new-package"))
  utils::write.csv(counts, file.path(output_dir, "bioccheck-summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
  print(counts, row.names = FALSE)
  if (any(counts$ERROR > 0L) || any(counts$WARNING > 0L)) {
    stop("BiocCheck reported ERROR or WARNING results.", call. = FALSE)
  }
  cat("Complete package validation passed. Archive:", archive, "\n")
  invisible(counts)
}

if (sys.nframe() == 0L) package_check_main()
