test_that("native BigWig wrappers return GeneTrackR signal schema", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  observed <- GeneTrackR:::query_bigwig_native(
    file = file,
    sample_id = "sampleA",
    chrom = "chr1",
    start = 11L,
    end = 105L,
    strand = "*"
  )

  expect_identical(
    names(observed),
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  expect_true(all(observed$sample_id == "sampleA"))
  expect_true(all(observed$strand == "*"))
  expect_equal(observed$start, c(11L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 105L))
})


test_that("native full-memory BigWig wrapper reads all chromosomes", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  observed <- GeneTrackR:::read_bigwig_whole_native(
    file = file,
    sample_id = "sampleA",
    strand = "*"
  )

  expect_identical(
    names(observed),
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  expect_equal(unique(observed$chrom), c("chr1", "chr2"))
  expect_true(all(observed$sample_id == "sampleA"))
})


test_that("read_bwg uses native R BigWig metadata and memory reader", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  expect_warning(
    observed <- read_bwg(
      file,
      format = "bigwig",
      sample_names = "sampleA",
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )

  expect_identical(observed$meta$backend, "GeneTrackR-native-R")
  expect_equal(
    observed$seqinfo,
    data.table::data.table(
      sample_id = c("sampleA", "sampleA"),
      chrom = c("chr1", "chr2"),
      length = c(1000L, 500L)
    )
  )
  expect_equal(
    observed$data[chrom == "chr1", .(start, end, value)],
    data.table::data.table(
      start = c(10L, 20L, 100L),
      end = c(12L, 20L, 110L),
      value = c(1.5, 2.5, -3.25)
    ),
    tolerance = 1e-6
  )
})


test_that("retrieve_bwg lazy BigWig queries use the native R path", {
  file <- test_path("fixtures", "native_reader_test.bigwig")
  track <- read_bwg(
    file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy",
    verbose = FALSE
  )

  observed <- retrieve_bwg(
    track,
    chrom = "chr1",
    start = 11L,
    end = 105L
  )

  expect_equal(observed$start, c(11L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 105L))
  expect_equal(observed$value, c(1.5, 2.5, -3.25), tolerance = 1e-6)
})


test_that("legacy BwgTrack seqinfo fallback uses native R metadata", {
  file <- test_path("fixtures", "native_reader_test.bigwig")
  legacy <- BwgTrack(
    samples = data.table::data.table(
      sample_id = "sampleA",
      file = normalizePath(file, winslash = "/", mustWork = TRUE),
      format = "bigwig",
      strand = "*"
    ),
    data = NULL,
    meta = list(mode = "lazy", coordinate = "1-based closed"),
    seqinfo = NULL
  )

  observed <- seqinfo_bwg(legacy)
  expect_equal(
    observed,
    data.table::data.table(
      sample_id = c("sampleA", "sampleA"),
      chrom = c("chr1", "chr2"),
      length = c(1000L, 500L)
    )
  )
})


test_that("public BigWig read paths no longer call the compiled reader", {
  read_body <- paste(deparse(body(read_bwg)), collapse = "\n")
  query_body <- paste(deparse(body(GeneTrackR:::.query_bwg_internal)), collapse = "\n")
  seqinfo_body <- paste(deparse(body(seqinfo_bwg)), collapse = "\n")

  expect_false(grepl("read_bigwig_whole_cpp", read_body, fixed = TRUE))
  expect_false(grepl("bw_seqinfo_cpp", read_body, fixed = TRUE))
  expect_false(grepl("query_bigwig_cpp", query_body, fixed = TRUE))
  expect_false(grepl("bw_seqinfo_cpp", seqinfo_body, fixed = TRUE))
})
