test_that("unified memory reader returns the canonical signal schema", {
  bedgraph <- tempfile(fileext = ".bedgraph")
  writeLines(c("chr1\t9\t12\t1.5", "chr1\t19\t20\t2.5"), bedgraph)

  wig <- tempfile(fileext = ".wig")
  writeLines(c(
    "variableStep chrom=chr1 span=3",
    "10\t1.5",
    "variableStep chrom=chr1 span=1",
    "20\t2.5"
  ), wig)

  bigwig <- test_path("fixtures", "native_reader_test.bigwig")
  expected_columns <- c("sample_id", "chrom", "start", "end", "value", "strand")

  observed_bg <- GeneTrackR:::read_signal_file_memory(
    bedgraph, "bedgraph", "sampleA", "+"
  )
  observed_wig <- GeneTrackR:::read_signal_file_memory(
    wig, "wig", "sampleA", "*"
  )
  observed_bw <- GeneTrackR:::read_signal_file_memory(
    bigwig, "bigwig", "sampleA", "*"
  )

  expect_identical(names(observed_bg), expected_columns)
  expect_identical(names(observed_wig), expected_columns)
  expect_identical(names(observed_bw), expected_columns)
  expect_equal(observed_bg[, .(chrom, start, end, value)], observed_wig[, .(chrom, start, end, value)])
})


test_that("unified text readers support compressed bedGraph and WIG files", {
  bedgraph <- tempfile(fileext = ".bedgraph.gz")
  con <- gzfile(bedgraph, open = "wt")
  writeLines(c("track type=bedGraph", "chr1\t9\t12\t1.5", "chr2\t4\t7\t-2"), con)
  close(con)

  wig <- tempfile(fileext = ".wig.gz")
  con <- gzfile(wig, open = "wt")
  writeLines(c(
    "fixedStep chrom=chr1 start=10 step=10 span=3",
    "1.5",
    "2.5"
  ), con)
  close(con)

  bg <- GeneTrackR:::read_signal_file_memory(bedgraph, "bedgraph", "bg", "+")
  wg <- GeneTrackR:::read_signal_file_memory(wig, "wig", "wg", "*")

  expect_equal(bg$start, c(10L, 5L))
  expect_equal(bg$end, c(12L, 7L))
  expect_equal(wg$start, c(10L, 20L))
  expect_equal(wg$end, c(12L, 22L))
})


test_that("unified regional dispatcher clips text and BigWig queries consistently", {
  bedgraph <- tempfile(fileext = ".bedgraph")
  writeLines(c("chr1\t9\t12\t1.5", "chr1\t19\t25\t2.5"), bedgraph)

  bg_sample <- data.table::data.table(
    sample_id = "sampleA",
    file = bedgraph,
    format = "bedgraph",
    strand = "+",
    use_tabix = FALSE,
    tabix_backend = NA_character_,
    tabix_empty_fallback = FALSE
  )
  bg <- GeneTrackR:::query_signal_file_region(bg_sample, "chr1", 11L, 22L)
  expect_equal(bg$start, c(11L, 20L))
  expect_equal(bg$end, c(12L, 22L))

  bw_sample <- data.table::data.table(
    sample_id = "sampleA",
    file = test_path("fixtures", "native_reader_test.bigwig"),
    format = "bigwig",
    strand = "*",
    use_tabix = FALSE,
    tabix_backend = NA_character_,
    tabix_empty_fallback = FALSE
  )
  bw <- GeneTrackR:::query_signal_file_region(bw_sample, "chr1", 11L, 105L)
  expect_equal(bw$start, c(11L, 20L, 100L))
  expect_equal(bw$end, c(12L, 20L, 105L))
})


test_that("public signal I/O paths use unified dispatchers", {
  read_body <- paste(deparse(body(read_bwg)), collapse = "\n")
  query_body <- paste(deparse(body(GeneTrackR:::.query_bwg_internal)), collapse = "\n")
  write_body <- paste(deparse(body(write_bwg)), collapse = "\n")

  expect_match(read_body, "read_signal_file_memory", fixed = TRUE)
  expect_match(read_body, "build_signal_seqinfo", fixed = TRUE)
  expect_false(grepl("read_bedgraph_file", read_body, fixed = TRUE))
  expect_false(grepl("read_wig_file", read_body, fixed = TRUE))
  expect_false(grepl("read_bigwig_whole_native", read_body, fixed = TRUE))

  expect_match(query_body, "query_signal_file_region", fixed = TRUE)
  expect_false(grepl("query_bigwig_native", query_body, fixed = TRUE))
  expect_false(grepl("query_bedgraph_tabix", query_body, fixed = TRUE))

  expect_match(write_body, "write_signal_file_memory", fixed = TRUE)
  expect_match(write_body, "signal_output_path", fixed = TRUE)
  expect_false(grepl("write_bedgraph_table", write_body, fixed = TRUE))
  expect_false(grepl("write_wig_table", write_body, fixed = TRUE))
  expect_false(grepl("write_bigwig_native", write_body, fixed = TRUE))
})


test_that("unified writer round-trips bedGraph and WIG without coordinate drift", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(10L, 20L, 5L),
    end = c(12L, 22L, 7L),
    value = c(1.5, 2.5, -3.25),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(
      sample_id = "sampleA",
      strand = "*",
      norm_method = "none"
    ),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )

  for (format in c("bedgraph", "wig")) {
    outdir <- tempfile(pattern = paste0("gtr_io_", format, "_"))
    dir.create(outdir)
    written <- write_bwg(
      track,
      outdir = outdir,
      format = format,
      overwrite = TRUE
    )
    observed <- read_bwg(
      written$file,
      format = format,
      sample_names = "sampleA",
      mode = "memory",
      verbose = FALSE
    )$data
    data.table::setorderv(observed, c("sample_id", "chrom", "start", "end"))
    expected <- data.table::copy(signal)
    data.table::setorderv(expected, c("sample_id", "chrom", "start", "end"))
    expect_equal(observed[, .(sample_id, chrom, start, end, value)], expected[, .(sample_id, chrom, start, end, value)])
  }
})
