# Author: Rensc
# Date: 2026-11-20
# Version: dev002
# Function: Test record accounting and failure visibility in reader evaluation.
# Input: The standalone comparison helper and deterministic in-memory fixtures.
# Output: testthat results; no package or backend changes.
tool <- file.path("..", "compare_standard_io.R")
if (file.exists(tool)) source(tool, local = TRUE)

testthat::test_that("reader comparison aligns keys and retains missing values", {
    native <- data.frame(id = c("a", "b"), start = c(1L, 5L),
        score = c(NA_real_, 2.5), strand = c("+", "*"))
    standard <- native[c(2L, 1L), ]
    result <- compare_io_tables(native, standard, "id", c("start", "score", "strand"))
    testthat::expect_identical(result$summary$status, "equivalent")
    testthat::expect_equal(result$summary$matched_rows, 2L)
    testthat::expect_equal(nrow(result$differences), 0L)
    testthat::expect_equal(nrow(result$native_only), 0L)
    standard$score[2L] <- 0
    result <- compare_io_tables(native, standard, "id", c("start", "score", "strand"))
    testthat::expect_identical(result$summary$status, "different")
    testthat::expect_identical(result$differences$field, "score")
    testthat::expect_identical(result$differences$native, NA_character_)
    testthat::expect_identical(result$differences$standard, "0")
})

testthat::test_that("reader comparison reports every unmatched record", {
    native <- data.frame(id = c("a", "b"), value = c(1, 2))
    standard <- data.frame(id = c("b", "c", "d"), value = c(2, 3, 4))
    result <- compare_io_tables(native, standard, "id", "value")
    testthat::expect_equal(result$summary$native_rows, 2L)
    testthat::expect_equal(result$summary$standard_rows, 3L)
    testthat::expect_equal(result$summary$matched_rows, 1L)
    testthat::expect_equal(result$summary$native_only_rows, 1L)
    testthat::expect_equal(result$summary$standard_only_rows, 2L)
    testthat::expect_identical(result$native_only$id, "a")
    testthat::expect_identical(result$standard_only$id, c("c", "d"))
    testthat::expect_equal(nrow(result$differences), 0L)
    testthat::expect_identical(result$summary$status, "different")
})

testthat::test_that("ambiguous keys and invalid schemas fail instead of joining", {
    input <- data.frame(id = c("a", "a"), value = c(1, 2))
    testthat::expect_error(compare_io_tables(input, input, "id", "value"),
        "Duplicate")
    input$id <- c("a", NA)
    testthat::expect_error(compare_io_tables(input, input, "id", "value"),
        "missing")
    input$id <- c("a", "b")
    testthat::expect_error(compare_io_tables(input, input, "id", "absent"),
        "columns")
    testthat::expect_error(compare_io_tables(input, input, "id", "value",
        tolerance = -1), "tolerance")
})

testthat::test_that("numeric tolerance does not hide missing or infinite values", {
    native <- data.frame(id = letters[seq_len(3L)], value = c(1, NA, Inf))
    standard <- data.frame(id = native$id, value = c(1 + 1e-8, NA, Inf))
    result <- compare_io_tables(native, standard, "id", "value", tolerance = 1e-6)
    testthat::expect_identical(result$summary$status, "equivalent")
    standard$value[3L] <- -Inf
    result <- compare_io_tables(native, standard, "id", "value", tolerance = 1e-6)
    testthat::expect_identical(result$differences$field, "value")
    testthat::expect_equal(nrow(result$differences), 1L)
})

testthat::test_that("reader measurements expose warnings and errors", {
    measured <- measure_io_reader(function() data.frame(value = seq_len(3L)),
        iterations = 2L)
    testthat::expect_identical(measured$status, "ok")
    testthat::expect_equal(nrow(measured$timing), 2L)
    testthat::expect_true(all(measured$timing$object_bytes > 0))
    warned <- measure_io_reader(function() {
        warning("Fixture advisory")
        data.frame(value = 1)
    }, iterations = 1L)
    testthat::expect_match(warned$warnings, "Fixture advisory")
    failed <- measure_io_reader(function() stop("Fixture failure"), iterations = 1L)
    testthat::expect_identical(failed$status, "error")
    testthat::expect_match(failed$error, "Fixture failure")
    testthat::expect_equal(nrow(failed$timing), 0L)
    testthat::expect_error(measure_io_reader(identity, iterations = 0L),
        "iterations")
})

testthat::test_that("range projection retains list metadata without S3 coercion", {
    ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 5L),
        strand = "+")
    S4Vectors::mcols(ranges)$Parent <- IRanges::CharacterList(c("t1", "t2"))
    table <- io_ranges_frame(ranges)
    testthat::expect_identical(table$chrom, "chr1")
    testthat::expect_identical(table$start, 1L)
    testthat::expect_identical(io_column(table, "Parent"), "t1,t2")
})

testthat::test_that("local-path reader calls restore the working directory on error", {
    file <- tempfile(fileext = ".bw")
    writeLines("fixture", file)
    on.exit(unlink(file), add = TRUE)
    before <- getwd()
    value <- io_with_local_path(file, function(path) {
        testthat::expect_identical(path, basename(file))
        testthat::expect_true(file.exists(path))
        "read"
    })
    testthat::expect_identical(value, "read")
    testthat::expect_identical(getwd(), before)
    testthat::expect_error(io_with_local_path(file,
        function(path) stop("Fixture local error")), "Fixture local error")
    testthat::expect_identical(getwd(), before)
})

testthat::test_that("BED12 projection exposes absent native block metadata", {
    ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 20L),
        strand = "+")
    S4Vectors::mcols(ranges)$name <- "t1"
    S4Vectors::mcols(ranges)$score <- 0
    S4Vectors::mcols(ranges)$thick <- IRanges::IRanges(2L, 18L)
    S4Vectors::mcols(ranges)$blocks <- IRanges::IRangesList(
        IRanges::IRanges(c(1L, 10L), c(5L, 20L)))
    standard <- io_annotation_table(ranges, "bed12", "rtracklayer")
    testthat::expect_identical(standard$block_sizes, "5,11")
    testthat::expect_identical(standard$block_starts, "0,9")
    testthat::expect_identical(standard$thick_start, 1L)
    native <- io_annotation_table(list(data = data.frame(chrom = "chr1",
        start = 1L, end = 20L, strand = "+", name = "t1", score = 0)),
        "bed12", "native")
    testthat::expect_true(is.na(native$block_sizes))
    result <- compare_io_tables(native, standard,
        c("chrom", "start", "end", "strand", "name"),
        c("score", "block_sizes", "block_starts", "thick_start", "thick_end"))
    testthat::expect_equal(result$summary$different_values, 4L)
})

testthat::test_that("native projection uses preserved BED12 and complete parent metadata", {
    data <- data.frame(chrom = "chr1", start = 1L, end = 20L,
        strand = "+", name = "t1", score = 0, thick_start = 2L, thick_end = 18L)
    data$block_sizes <- list(c(5L, 11L))
    data$block_starts <- list(c(0L, 9L))
    projected <- io_annotation_table(list(data = data), "bed12", "native")
    testthat::expect_identical(projected$block_sizes, "5,11")
    testthat::expect_identical(projected$block_starts, "0,9")
    testthat::expect_identical(projected$thick_start, 1L)
    testthat::expect_identical(projected$thick_end, 18L)
    data$type <- "exon"
    data$source <- "fixture"
    data$phase <- NA_character_
    data$feature_id <- "e1"
    data$parent_id <- "t1"
    data$parent_ids <- list(c("t1", "t2"))
    projected <- io_annotation_table(list(data = data), "gff3", "native")
    testthat::expect_identical(projected$parent, "t1,t2")
})

testthat::test_that("INFO projection distinguishes flags from missing numeric fields", {
    input <- c("DP=12;AF=0.25,0.5;FLAG", "DP=.;AF=.", ".")
    testthat::expect_identical(io_info_values(input, "DP"), c("12", NA, NA))
    testthat::expect_identical(io_info_values(input, "AF"), c("0.25,0.5", NA, NA))
    testthat::expect_identical(io_info_values(input, "FLAG"), c("TRUE", "FALSE", "FALSE"))
})

testthat::test_that("comparison CLI rejects incomplete and destructive report requests", {
    testthat::expect_output(standard_io_main("--help"), "Usage:")
    testthat::expect_error(standard_io_main("--unknown"), "Usage:")
    testthat::expect_error(standard_io_main(c("--iterations", "0")), "positive integer")
    testthat::expect_error(standard_io_main(c("--output-dir", tempdir())), "fresh")
})
