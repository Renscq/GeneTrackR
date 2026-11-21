test_that("lazy indexed bedGraph clips hits and reads the full file only after index failure", {
  skip_if_not_installed("Rsamtools")

  source <- tempfile(fileext = ".bedgraph")
  compressed <- paste0(source, ".gz")
  index <- paste0(compressed, ".tbi")
  on.exit(unlink(c(source, compressed, index)), add = TRUE)
  writeLines(c(
    "chr1\t0\t4\t9",
    "chr1\t4\t10\t1.5",
    "chr1\t15\t25\t2.5",
    "chr1\t25\t29\t-4",
    "chr2\t0\t20\t7"
  ), source)
  Rsamtools::bgzip(source, dest = compressed)
  Rsamtools::indexTabix(compressed, format = "bed")

  full_file_reads <- 0L
  suppressMessages(trace(
    "read_signal_file_memory",
    tracer = function() full_file_reads <<- full_file_reads + 1L,
    print = FALSE,
    where = asNamespace("GeneTrackR")
  ))
  on.exit(suppressMessages(untrace(
    "read_signal_file_memory", where = asNamespace("GeneTrackR")
  )), add = TRUE)

  track <- read_bwg(compressed, mode = "lazy", sample_names = "signal",
                    strand = "+", verbose = FALSE)
  expect_s3_class(track, "BwgTrack")
  expect_null(track$data)
  expect_true(track$samples$has_tabix)
  expect_true(track$samples$use_tabix)
  # Exercise Rsamtools even when a system tabix executable is available.
  track$samples$tabix_backend <- "Rsamtools"

  expected <- list(
    sample_id = rep("signal", 2L),
    chrom = rep("chr1", 2L),
    start = c(7L, 16L),
    end = c(10L, 20L),
    value = c(1.5, 2.5),
    strand = rep("+", 2L)
  )
  hits <- retrieve_bwg(track, "chr1", 7L, 20L)
  expect_s3_class(hits, "data.table")
  expect_identical(as.list(hits), expected)
  expect_identical(full_file_reads, 0L)

  full_file_reads <- 0L
  gap <- retrieve_bwg(track, "chr1", 11L, 15L)
  expect_identical(as.list(gap), lapply(expected, function(column) column[0L]))
  expect_identical(full_file_reads, 0L)

  full_file_reads <- 0L
  expect_identical(unlink(index), 0L)
  expect_false(file.exists(index))
  recovered <- retrieve_bwg(track, "chr1", 7L, 20L)
  expect_identical(as.list(recovered), expected)
  expect_identical(full_file_reads, 1L)
})
