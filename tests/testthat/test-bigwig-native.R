test_that("native BigWig binary helpers decode integer and float values", {
  raw_u32 <- as.raw(c(0x26, 0xfc, 0x8f, 0x88))
  expect_equal(
    GeneTrackR:::.gtr_bw_u32(raw_u32, 0L, "little"),
    GeneTrackR:::.GTR_BIGWIG_MAGIC
  )

  raw_float <- as.raw(c(0x00, 0x00, 0xc0, 0x3f))
  expect_equal(
    GeneTrackR:::.gtr_bw_float32_vector(raw_float, 0L, "little"),
    1.5,
    tolerance = 1e-7
  )
})


test_that("native BigWig metadata reader resolves chromosome lengths", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  metadata <- GeneTrackR:::bigwig_metadata_native(file, use_cache = FALSE)
  observed <- GeneTrackR:::bigwig_seqinfo_native(file)
  expected <- data.table::data.table(
    chrom = c("chr1", "chr2"),
    length = c(1000L, 500L)
  )

  expect_identical(metadata$header$endian, "little")
  expect_equal(metadata$header$version, 4L)
  expect_equal(observed, expected)
})


test_that("native BigWig query returns 1-based closed intervals", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  observed <- GeneTrackR:::bigwig_query_native(
    file = file,
    chrom = "chr1",
    start = 1L,
    end = 200L
  )

  expect_equal(observed$start, c(10L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 110L))
  expect_equal(observed$value, c(1.5, 2.5, -3.25), tolerance = 1e-6)
})


test_that("native BigWig query clips intervals to requested boundaries", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  observed <- GeneTrackR:::bigwig_query_native(
    file = file,
    chrom = "chr1",
    start = 11L,
    end = 105L
  )

  expect_equal(observed$start, c(11L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 105L))
})


test_that("native BigWig query returns an empty table outside chromosome range", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  observed <- GeneTrackR:::bigwig_query_native(
    file = file,
    chrom = "chr2",
    start = 501L,
    end = 600L
  )

  expect_s3_class(observed, "data.table")
  expect_equal(nrow(observed), 0L)
  expect_identical(
    names(observed),
    c("chrom", "start", "end", "value")
  )
})


test_that("native BigWig reader rejects invalid binary input", {
  file <- tempfile(fileext = ".bigwig")
  writeBin(as.raw(rep(0L, 64L)), file)

  expect_error(
    GeneTrackR:::bigwig_metadata_native(file, use_cache = FALSE),
    "not a valid BigWig"
  )
})


test_that("native BigWig adapters preserve the low-level reader result", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  low_level <- GeneTrackR:::bigwig_query_native(
    file,
    chrom = "chr1",
    start = 1L,
    end = 200L
  )
  adapted <- GeneTrackR:::query_bigwig_native(
    file = file,
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 200L,
    strand = "*"
  )

  expect_equal(
    adapted[, .(chrom, start, end, value)],
    low_level[, .(chrom, start, end, value)],
    tolerance = 1e-6
  )
  expect_true(all(adapted$sample_id == "sampleA"))
  expect_true(all(adapted$strand == "*"))
})
