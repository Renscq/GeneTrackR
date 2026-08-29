test_that("BwgTrack exposes the schema-v2 core slots without breaking positional construction", {
  samples <- data.table::data.table(sample_id = "sampleA")
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 10L,
    end = 20L,
    value = 1,
    strand = "*"
  )

  x <- BwgTrack(
    samples,
    signal,
    list(mode = "memory"),
    make_empty_validation()
  )

  expect_named(x, c("samples", "data", "seqinfo", "meta", "validation"))
  expect_null(x$seqinfo)
  expect_identical(x$meta$coordinate, "1-based closed")
  expect_identical(x$meta$schema_version, "2")
  expect_identical(x$meta$mode, "memory")
})


test_that("BwgTrack accepts and normalizes sample-specific seqinfo", {
  samples <- data.table::data.table(
    sample_id = c("sampleA", "sampleB")
  )
  seqinfo <- data.frame(
    sample_id = c("sampleB", "sampleA", "sampleA"),
    chrom = c("chr2", "chr2", "chr1"),
    length = c(2000, 2000, 1000)
  )

  x <- BwgTrack(samples = samples, seqinfo = seqinfo)

  expect_s3_class(x$seqinfo, "data.table")
  expect_identical(names(x$seqinfo)[1:3], c("sample_id", "chrom", "length"))
  expect_type(x$seqinfo$sample_id, "character")
  expect_type(x$seqinfo$chrom, "character")
  expect_type(x$seqinfo$length, "integer")
  expect_equal(
    x$seqinfo[, .(sample_id, chrom)],
    data.table::data.table(
      sample_id = c("sampleA", "sampleA", "sampleB"),
      chrom = c("chr1", "chr2", "chr2")
    )
  )
})


test_that("BwgTrack rejects inconsistent seqinfo", {
  samples <- data.table::data.table(sample_id = "sampleA")

  expect_error(
    BwgTrack(
      samples = samples,
      seqinfo = data.frame(sample_id = "sampleB", chrom = "chr1", length = 1000)
    ),
    "unknown sample IDs"
  )

  expect_error(
    BwgTrack(
      samples = samples,
      seqinfo = data.frame(
        sample_id = c("sampleA", "sampleA"),
        chrom = c("chr1", "chr1"),
        length = c(1000, 2000)
      )
    ),
    "conflicting chromosome lengths"
  )
})


test_that("memory bedGraph reading records schema-v2 seqinfo with unknown lengths", {
  file <- tempfile(fileext = ".bedgraph")
  writeLines(c(
    "chr1\t0\t10\t1",
    "chr1\t10\t20\t2",
    "chr2\t5\t15\t3"
  ), file)

  x <- read_bwg(
    file,
    format = "bedgraph",
    sample_names = "sampleA",
    mode = "memory",
    verbose = FALSE
  )

  expect_identical(x$meta$schema_version, "2")
  expect_equal(
    x$seqinfo[, .(sample_id, chrom)],
    data.table::data.table(
      sample_id = c("sampleA", "sampleA"),
      chrom = c("chr1", "chr2")
    )
  )
  expect_true(all(is.na(x$seqinfo$length)))
  expect_equal(seqinfo_bwg(x), x$seqinfo)
})


test_that("BwgTrack subsets preserve relevant seqinfo", {
  samples <- data.table::data.table(
    sample_id = c("sampleA", "sampleB"),
    file = NA_character_,
    format = "memory",
    strand = "*",
    has_strand = FALSE,
    scale_factor = 1
  )
  signal <- data.table::data.table(
    sample_id = c("sampleA", "sampleA", "sampleB"),
    chrom = c("chr1", "chr2", "chr1"),
    start = c(10L, 30L, 15L),
    end = c(20L, 40L, 25L),
    value = c(1, 2, 3),
    strand = "*"
  )
  seqinfo <- data.table::data.table(
    sample_id = c("sampleA", "sampleA", "sampleB", "sampleB"),
    chrom = c("chr1", "chr2", "chr1", "chr2"),
    length = c(1000L, 2000L, 1000L, 2000L)
  )
  x <- BwgTrack(samples = samples, data = signal, seqinfo = seqinfo)

  retrieved <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    samples = "sampleA",
    as = "BwgTrack"
  )
  sliced <- slice_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    samples = "sampleA",
    as = "BwgTrack"
  )

  expected_seqinfo <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    length = 1000L
  )
  expect_equal(retrieved$seqinfo, expected_seqinfo)
  expect_equal(sliced$seqinfo, expected_seqinfo)
})


test_that("merge_bwg preserves seqinfo for independent samples", {
  x1 <- BwgTrack(
    samples = data.table::data.table(
      sample_id = "sampleA",
      norm_method = "none"
    ),
    data = data.table::data.table(
      sample_id = "sampleA", chrom = "chr1", start = 1L,
      end = 10L, value = 1, strand = "*"
    ),
    seqinfo = data.table::data.table(
      sample_id = "sampleA", chrom = "chr1", length = 1000L
    )
  )
  x2 <- BwgTrack(
    samples = data.table::data.table(
      sample_id = "sampleB",
      norm_method = "none"
    ),
    data = data.table::data.table(
      sample_id = "sampleB", chrom = "chr2", start = 1L,
      end = 10L, value = 2, strand = "*"
    ),
    seqinfo = data.table::data.table(
      sample_id = "sampleB", chrom = "chr2", length = 2000L
    )
  )

  merged <- merge_bwg(x1, x2)

  expect_equal(
    merged$seqinfo,
    data.table::data.table(
      sample_id = c("sampleA", "sampleB"),
      chrom = c("chr1", "chr2"),
      length = c(1000L, 2000L)
    )
  )
  expect_identical(merged$meta$schema_version, "2")
})
