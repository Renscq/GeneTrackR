# Author: Rensc
# Date: 2026-11-20
# Version: dev002
# Function: Exercise the complete optional shared-fixture reader assessment.
# Input: Installed evaluation dependencies and the matching package checkout.
# Output: Full harness regressions with record and report accounting.
source(file.path("..", "compare_standard_io.R"), local = TRUE)

testthat::test_that("shared-fixture assessment retains all readers and records", {
    for (package in c("GeneTrackR", "Bioc.gff", "rtracklayer", "VariantAnnotation")) {
        testthat::skip_if_not_installed(package)
    }
    output <- tempfile("gtr-io-assessment-")
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    source <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
    result <- standard_io_main(c("--path", source, "--output-dir", output,
        "--iterations", "1"))
    status <- utils::read.csv(file.path(output, "reader-status.csv"))
    testthat::expect_equal(nrow(status), 14L)
    testthat::expect_true(all(status$status == "ok"))
    testthat::expect_equal(nrow(result), 8L)
    testthat::expect_true(all(result$native_only_rows == 0L))
    testthat::expect_true(all(result$standard_only_rows == 0L))
    expected <- c(gff3 = 8L, gtf = 8L, bed6 = 3L, bed12 = 2L, vcf = 4L, bigwig = 4L)
    testthat::expect_equal(result$matched_rows, unname(expected[result$format]))
    testthat::expect_true(all(result$status[result$format %in% c("gff3", "gtf", "bed6", "bed12", "bigwig")] ==
        "equivalent"))
    difference <- utils::read.csv(file.path(output, "bed12-rtracklayer-differences.csv"))
    testthat::expect_equal(nrow(difference), 0L)
    order <- utils::read.csv(file.path(output, "vcf-sample-order.csv"))
    testthat::expect_identical(order$native, c("S2", "S1"))
    testthat::expect_identical(order$standard, order$native)
    timing <- utils::read.csv(file.path(output, "timing.csv"))
    testthat::expect_equal(nrow(timing), 14L)
    testthat::expect_true(all(timing$object_bytes > 0))
    testthat::expect_true(file.exists(file.path(output, "source-checksums.csv")))
    testthat::expect_error(standard_io_main(c("--path", source,
        "--output-dir", output)), "fresh")
})
