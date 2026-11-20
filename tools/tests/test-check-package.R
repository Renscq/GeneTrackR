# Author: Rensc
# Date: 2026-11-20
# Version: dev002
# Function: Regress safe package-check report paths outside the source tree.
# Input: tools/check_package.R sourced by the standalone test runner.
# Output: testthat assertions; temporary fixture directories only.

source(file.path("..", "check_package.R"), local = TRUE)

testthat::test_that("package checks reject report directories inside the source", {
  fixture <- tempfile("gtr-check-path-")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  source_path <- file.path(fixture, "GeneTrackR")
  dir.create(source_path)
  writeLines(c("Package: GeneTrackR", "Version: 0.99.8"),
    file.path(source_path, "DESCRIPTION"))

  testthat::expect_error(package_check_main(c("--path", source_path,
    "--output-dir", source_path)), "outside the package source")
  testthat::expect_error(package_check_main(c("--path", source_path,
    "--output-dir", file.path(source_path, "temp", "checks"))),
    "outside the package source")
  testthat::expect_false(file.exists(file.path(source_path, "build.log")))
})

testthat::test_that("default and explicit report paths are persistent siblings", {
  fixture <- tempfile("gtr-check-path-")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  source_path <- file.path(fixture, "GeneTrackR")
  dir.create(source_path)

  output <- package_check_output(source_path)
  testthat::expect_identical(basename(output), "GeneTrackR-package-check")
  testthat::expect_identical(dirname(output),
    normalizePath(fixture, winslash = "/", mustWork = TRUE))
  prefixed_sibling <- file.path(fixture, "GeneTrackR-reports")
  testthat::expect_identical(package_check_output(source_path, prefixed_sibling),
    normalizePath(prefixed_sibling, winslash = "/", mustWork = TRUE))
  testthat::expect_error(package_check_output(source_path,
    file.path(output, "..", "GeneTrackR", "nested")),
    "outside the package source")
  if (.Platform$OS.type == "windows") {
    testthat::expect_error(package_check_output(source_path, toupper(source_path)),
      "outside the package source")
  }
})

testthat::test_that("Support Site rate limiting retries complete checks and retains attempts", {
  output <- tempfile("gtr-bioc-retry-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  warning_result <- list(error = list(), warning = list(checkSupportReg = list(
    list(checkSupportReg = "Unable to retrieve email info from the Support Site:\n  HTTP status 429"))),
    note = list("Existing advisory"))
  clean_result <- list(error = list(), warning = list(), note = list("Existing advisory"))
  calls <- 0L
  checker <- function(package, checkDir, install, `new-package`, `quit-with-status`) {
    calls <<- calls + 1L
    testthat::expect_identical(package, "fixture.tar.gz")
    testthat::expect_true(install)
    testthat::expect_true(`new-package`)
    testthat::expect_false(`quit-with-status`)
    if (calls == 1L) warning_result else clean_result
  }
  result <- run_new_package_check("fixture.tar.gz", output,
    checker = checker, retry_delays = c(0, 0))
  testthat::expect_identical(calls, 2L)
  testthat::expect_identical(result, clean_result)
  attempts <- utils::read.csv(file.path(output, "bioccheck-attempts.csv"))
  testthat::expect_identical(attempts$WARNING, c(1L, 0L))
  testthat::expect_identical(readRDS(file.path(output,
    "bioccheck-attempt-1", "result.rds")), warning_result)
})

testthat::test_that("unrelated checker errors and warnings are never retried", {
  support <- "Unable to retrieve email info from the Support Site:\n  HTTP status 429"
  cases <- list(
    list(error = list("Invalid source"), warning = list(support), note = list()),
    list(error = list(), warning = list(support, "Invalid namespace"), note = list()),
    list(error = list(), warning = list("Other service: HTTP status 429"), note = list()),
    list(error = list(), warning = list(), note = list("Existing advisory")))
  for (fixture in cases) {
    output <- tempfile("gtr-bioc-no-retry-")
    dir.create(output)
    calls <- 0L
    checker <- function(...) { calls <<- calls + 1L; fixture }
    observed <- run_new_package_check("fixture.tar.gz", output,
      checker = checker, retry_delays = c(0, 0))
    testthat::expect_identical(calls, 1L)
    testthat::expect_identical(observed, fixture)
    unlink(output, recursive = TRUE)
  }
})

testthat::test_that("persistent rate limiting exhausts a bounded number of full attempts", {
  output <- tempfile("gtr-bioc-retry-bound-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  fixture <- list(error = list(), warning = list(
    "Unable to retrieve email info from the Support Site:\n  HTTP status 429"), note = list())
  calls <- 0L
  checker <- function(...) { calls <<- calls + 1L; fixture }
  observed <- run_new_package_check("fixture.tar.gz", output,
    checker = checker, retry_delays = c(0, 0))
  testthat::expect_identical(calls, 3L)
  testthat::expect_identical(observed, fixture)
  testthat::expect_equal(nrow(utils::read.csv(file.path(output, "bioccheck-attempts.csv"))), 3L)
})

testthat::test_that("legacy Support Site HTTP 429 errors retry without hiding severity", {
  output <- tempfile("gtr-bioc-legacy-retry-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  fixture <- list(error = list(checkSupportReg = list(list(checkSupportReg =
    "Unable to find your email in the Support Site:\n  HTTP 429 Too Many Requests."))),
    warning = list(), note = list())
  clean <- list(error = list(), warning = list(), note = list())
  calls <- 0L
  checker <- function(...) { calls <<- calls + 1L; if (calls == 1L) fixture else clean }
  result <- run_new_package_check("fixture.tar.gz", output,
    checker = checker, retry_delays = c(0, 0))
  testthat::expect_identical(calls, 2L)
  testthat::expect_identical(result, clean)
  attempts <- utils::read.csv(file.path(output, "bioccheck-attempts.csv"))
  testthat::expect_identical(attempts$ERROR, c(1L, 0L))
  testthat::expect_identical(readRDS(file.path(output,
    "bioccheck-attempt-1", "result.rds")), fixture)
})

testthat::test_that("legacy rate limits stay errors after retry exhaustion or real failures", {
  rate_limit <- "Unable to find your email in the Support Site: HTTP 429 Too Many Requests."
  cases <- list(
    list(result = list(error = list(rate_limit), warning = list(), note = list()), attempts = 3L),
    list(result = list(error = list(rate_limit, "Invalid source"), warning = list(), note = list()), attempts = 1L),
    list(result = list(error = list("Unable to find your email in the Support Site: HTTP 404 Not Found."),
      warning = list(), note = list()), attempts = 1L),
    list(result = list(error = list(rate_limit), warning = list("Invalid namespace"), note = list()), attempts = 1L))
  for (fixture in cases) {
    output <- tempfile("gtr-bioc-legacy-failure-")
    dir.create(output)
    calls <- 0L
    checker <- function(...) { calls <<- calls + 1L; fixture$result }
    result <- run_new_package_check("fixture.tar.gz", output,
      checker = checker, retry_delays = c(0, 0))
    testthat::expect_identical(calls, fixture$attempts)
    testthat::expect_identical(result, fixture$result)
    unlink(output, recursive = TRUE)
  }
})
