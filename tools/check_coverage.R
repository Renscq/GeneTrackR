#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-11-19
# Version: dev001
# Function: Measure and report reproducible tests-only package coverage.
# Input: A package source path, output directory and optional coverage floor.
# Output: RDS, CSV, provenance and optional before/after comparison reports.

coverage_summary <- function(coverage) {
  tally <- covr::tally_coverage(coverage, by = "line")
  if (!nrow(tally)) stop("No instrumented package lines were found.", call. = FALSE)
  lines <- stats::aggregate(tally$value,
    by = list(filename = tally$filename, line = tally$line), FUN = min)
  names(lines)[3L] <- "value"
  lines <- lines[order(lines$filename, lines$line), , drop = FALSE]
  files <- do.call(rbind, lapply(split(lines, lines$filename), function(rows) {
    eligible <- !is.na(rows$value)
    total <- sum(eligible)
    covered <- sum(rows$value[eligible] > 0)
    data.frame(filename = rows$filename[1L], lines = total,
               covered = covered, coverage = 100 * covered / total)
  }))
  rownames(files) <- NULL
  overall <- data.frame(lines = sum(files$lines), covered = sum(files$covered),
                        coverage = 100 * sum(files$covered) / sum(files$lines))
  native <- covr::percent_coverage(coverage, by = "line")
  if (!isTRUE(all.equal(overall$coverage, native, tolerance = 1e-10))) {
    stop("Coverage aggregation disagrees with covr's line metric.", call. = FALSE)
  }
  list(lines = lines, files = files, overall = overall)
}

write_coverage_reports <- function(coverage, package_path, output_dir,
                                   baseline = NULL, source_ref = "working-tree") {
  package_path <- normalizePath(package_path, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- coverage_summary(coverage)
  desc <- read.dcf(file.path(package_path, "DESCRIPTION"))
  saveRDS(coverage, file.path(output_dir, "coverage.rds"))
  utils::write.csv(summary$lines, file.path(output_dir, "coverage-lines.csv"), row.names = FALSE)
  utils::write.csv(summary$files, file.path(output_dir, "coverage-files.csv"), row.names = FALSE)
  utils::write.csv(summary$overall, file.path(output_dir, "coverage-summary.csv"), row.names = FALSE)
  metadata <- c(Package = unname(desc[1L, "Package"]), Version = unname(desc[1L, "Version"]),
    Scope = "tests", Metric = "line", Exclusions = "none requested",
    Engine = paste("covr", utils::packageVersion("covr")), R = R.version.string,
    Platform = R.version$platform, SourceRef = source_ref,
    CreatedUTC = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"))
  write.dcf(as.data.frame(as.list(metadata)), file.path(output_dir, "coverage-metadata.dcf"))
  source_files <- c(list.files(file.path(package_path, "R"), "[.]R$", full.names = TRUE),
                    list.files(file.path(package_path, "tests"), full.names = TRUE, recursive = TRUE))
  relative <- substring(source_files, nchar(package_path) + 2L)
  utils::write.csv(data.frame(filename = relative, md5 = unname(tools::md5sum(source_files))),
                   file.path(output_dir, "source-checksums.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
  if (!is.null(baseline)) {
    old <- utils::read.csv(baseline, stringsAsFactors = FALSE)
    required <- c("filename", "lines", "covered", "coverage")
    if (!all(required %in% names(old)) || anyDuplicated(old$filename)) {
      stop("Baseline CSV has an invalid coverage schema.", call. = FALSE)
    }
    comparison <- merge(old[required], summary$files, by = "filename", all = TRUE,
                         suffixes = c("_before", "_after"), sort = TRUE)
    comparison$delta <- comparison$coverage_after - comparison$coverage_before
    utils::write.csv(comparison, file.path(output_dir, "coverage-comparison.csv"), row.names = FALSE)
  }
  summary
}

coverage_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- paste(
    "Usage: Rscript tools/check_coverage.R [options]",
    "  --path DIR          Package source (default: current directory)",
    "  --output-dir DIR    Reports (default: temp/coverage)",
    "  --minimum NUMBER    Minimum line coverage, 0 to 100 (default: 0)",
    "  --baseline FILE     Optional previous coverage-files.csv for comparison",
    "  --help              Show this help",
    "Scope is tests only; all runtime R code is measured without exclusions.", sep = "\n")
  if (any(args %in% c("--help", "-h"))) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  values <- list(path = ".", output_dir = "temp/coverage", minimum = "0", baseline = NULL)
  index <- 1L
  mapping <- c("--path" = "path", "--output-dir" = "output_dir",
                "--minimum" = "minimum", "--baseline" = "baseline")
  while (index <= length(args)) {
    flag <- args[[index]]
    if (!flag %in% names(mapping) || index == length(args)) {
      stop(paste("Invalid or incomplete option:", flag, "\n", help), call. = FALSE)
    }
    values[[mapping[[flag]]]] <- args[[index + 1L]]
    index <- index + 2L
  }
  minimum <- suppressWarnings(as.numeric(values$minimum))
  if (length(minimum) != 1L || !is.finite(minimum) || minimum < 0 || minimum > 100) {
    stop("--minimum must be a finite number from 0 to 100.", call. = FALSE)
  }
  if (!requireNamespace("covr", quietly = TRUE) || utils::packageVersion("covr") < "4.0.0") {
    stop("Install the developer dependency covr >= 4.0.0 to measure coverage.", call. = FALSE)
  }
  if (!is.null(values$baseline) && !file.exists(values$baseline)) {
    stop("The baseline CSV does not exist.", call. = FALSE)
  }
  dir.create(values$output_dir, recursive = TRUE, showWarnings = FALSE)
  values$path <- normalizePath(values$path, winslash = "/", mustWork = TRUE)
  values$output_dir <- normalizePath(values$output_dir, winslash = "/", mustWork = TRUE)
  coverage <- covr::package_coverage(path = values$path, type = "tests",
    relative_path = TRUE, quiet = TRUE, line_exclusions = NULL,
    function_exclusions = NULL, allow_failures = FALSE,
    outDir = file.path(values$output_dir, "test-output"))
  source_ref <- Sys.getenv("GITHUB_SHA", unset = "working-tree")
  summary <- write_coverage_reports(coverage, values$path, values$output_dir,
                                     values$baseline, source_ref)
  percent <- summary$overall$coverage
  cat(sprintf("Tests-only line coverage: %.2f%% (%d/%d lines)\n",
               percent, summary$overall$covered, summary$overall$lines))
  if (percent + 1e-10 < minimum) {
    stop(sprintf("Line coverage %.2f%% is below the %.2f%% floor.", percent, minimum), call. = FALSE)
  }
  invisible(summary)
}

if (sys.nframe() == 0L) coverage_main()
