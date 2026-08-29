test_that("write_bwg routes BigWig output through the native R dispatcher", {
  write_body <- paste(deparse(body(write_bwg)), collapse = "\n")
  dispatcher_body <- paste(
    deparse(body(GeneTrackR:::write_signal_file_memory)),
    collapse = "\n"
  )

  expect_match(write_body, "write_signal_file_memory", fixed = TRUE)
  expect_false(grepl("write_bigwig_libbigwig", write_body, fixed = TRUE))
  expect_match(dispatcher_body, "write_bigwig_native", fixed = TRUE)
  expect_false(grepl("write_bigwig_libbigwig", dispatcher_body, fixed = TRUE))
})


test_that("native R BigWig writer round-trips signal intervals", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(10L, 20L, 100L),
    end = c(12L, 22L, 101L),
    value = c(1.5, 2.5, -3.5),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    seqinfo = data.table::data.table(
      sample_id = "sampleA",
      chrom = c("chr1", "chr2"),
      length = c(1000L, 2000L)
    ),
    meta = list(mode = "memory")
  )

  outdir <- tempfile(pattern = "gtr_native_writer_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(1000L, 2000L)),
    overwrite = TRUE
  )

  expect_true(file.exists(written$file))
  expect_gt(file.info(written$file)$size, 0)
  metadata <- GeneTrackR:::bigwig_metadata_native(written$file, use_cache = FALSE)
  expect_gt(metadata$header$n_levels, 0L)

  observed_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy",
    verbose = FALSE
  )
  expect_equal(
    seqinfo_bwg(observed_track),
    data.table::data.table(
      sample_id = "sampleA",
      chrom = c("chr1", "chr2"),
      length = c(1000L, 2000L)
    )
  )

  observed <- data.table::rbindlist(list(
    retrieve_bwg(observed_track, chrom = "chr1", start = 1L, end = 1000L),
    retrieve_bwg(observed_track, chrom = "chr2", start = 1L, end = 2000L)
  ))
  data.table::setorderv(observed, c("chrom", "start", "end"))
  expected <- data.table::copy(signal)
  data.table::setorderv(expected, c("chrom", "start", "end"))

  expect_equal(
    observed[, .(sample_id, chrom, start, end)],
    expected[, .(sample_id, chrom, start, end)]
  )
  expect_equal(observed$value, expected$value, tolerance = 1e-4)
})


test_that("native R BigWig writer supports multiple samples", {
  signal <- data.table::data.table(
    sample_id = c("plus", "plus", "minus", "minus"),
    chrom = "chr1",
    start = c(10L, 20L, 40L, 50L),
    end = c(12L, 22L, 42L, 52L),
    value = c(10, 20, 15, 25),
    strand = c("+", "+", "-", "-")
  )
  track <- BwgTrack(
    samples = data.table::data.table(
      sample_id = c("plus", "minus"),
      strand = c("+", "-")
    ),
    data = signal,
    meta = list(mode = "memory")
  )

  outdir <- tempfile(pattern = "gtr_native_multi_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = "chr1", size = 1000L),
    overwrite = TRUE
  )

  expect_equal(written$sample_id, c("plus", "minus"))
  expect_true(all(file.exists(written$file)))

  for (sid in c("plus", "minus")) {
    observed_track <- read_bwg(
      written$file[written$sample_id == sid],
      format = "bigwig",
      sample_names = sid,
      mode = "lazy",
      verbose = FALSE
    )
    observed <- retrieve_bwg(
      observed_track,
      chrom = "chr1",
      start = 1L,
      end = 1000L
    )
    expected <- signal[sample_id == sid]
    expect_equal(
      observed[, .(chrom, start, end)],
      expected[, .(chrom, start, end)]
    )
    expect_equal(observed$value, expected$value, tolerance = 1e-4)
  }
})


test_that("native R BigWig writer rejects overlaps and chromosome overflow", {
  overlap <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA"),
    data = data.table::data.table(
      sample_id = "sampleA",
      chrom = "chr1",
      start = c(10L, 12L),
      end = c(15L, 20L),
      value = c(1, 2),
      strand = "*"
    ),
    meta = list(mode = "memory")
  )
  expect_error(
    write_bwg(
      overlap,
      outdir = tempfile(pattern = "gtr_overlap_"),
      format = "bigwig",
      chrom_sizes = data.frame(chrom = "chr1", size = 1000L),
      overwrite = TRUE
    ),
    "overlapping intervals"
  )

  overflow <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA"),
    data = data.table::data.table(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 95L,
      end = 110L,
      value = 1,
      strand = "*"
    ),
    meta = list(mode = "memory")
  )
  expect_error(
    write_bwg(
      overflow,
      outdir = tempfile(pattern = "gtr_overflow_"),
      format = "bigwig",
      chrom_sizes = data.frame(chrom = "chr1", size = 100L),
      overwrite = TRUE
    ),
    "beyond chromosome sizes"
  )
})


test_that("native R BigWig writer supports multi-level chromosome indexes", {
  n_chrom <- 300L
  chrom <- sprintf("chr%03d", seq_len(n_chrom))
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = chrom,
    start = 11L,
    end = 20L,
    value = as.numeric((seq_len(n_chrom) %% 7L) - 3L),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA"),
    data = signal,
    meta = list(mode = "memory")
  )
  outdir <- tempfile(pattern = "gtr_multilevel_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = chrom, size = rep(1000L, n_chrom)),
    overwrite = TRUE
  )

  observed <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy",
    verbose = FALSE
  )
  si <- seqinfo_bwg(observed)
  expect_equal(nrow(si), n_chrom)
  expect_equal(si$chrom[c(1L, n_chrom)], chrom[c(1L, n_chrom)])
  expect_equal(
    retrieve_bwg(observed, chrom = chrom[1L], start = 1L, end = 1000L)$value,
    signal$value[1L],
    tolerance = 1e-4
  )
  expect_equal(
    retrieve_bwg(observed, chrom = chrom[n_chrom], start = 1L, end = 1000L)$value,
    signal$value[n_chrom],
    tolerance = 1e-4
  )
})
