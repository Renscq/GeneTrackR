.signal_frame_for_test <- function(x, columns = c("sample_id", "chrom", "start", "end", "value", "strand")) {
  x <- data.table::copy(data.table::as.data.table(x))
  columns <- intersect(columns, names(x))
  if (nrow(x) > 0L) {
    order_columns <- intersect(c("sample_id", "chrom", "start", "end"), names(x))
    if (length(order_columns) > 0L) {
      data.table::setorderv(x, order_columns)
    }
  }
  as.data.frame(x[, columns, with = FALSE], stringsAsFactors = FALSE)
}


test_that("lazy and memory BigWig retrieval are equivalent", {
  file <- test_path("fixtures", "native_reader_test.bigwig")

  lazy <- read_bwg(
    file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy",
    verbose = FALSE
  )
  expect_warning(
    memory <- read_bwg(
      file,
      format = "bigwig",
      sample_names = "sampleA",
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )

  lazy_region <- retrieve_bwg(lazy, chrom = "chr1", start = 11L, end = 105L)
  memory_region <- retrieve_bwg(memory, chrom = "chr1", start = 11L, end = 105L)

  expect_equal(
    .signal_frame_for_test(lazy_region),
    .signal_frame_for_test(memory_region),
    tolerance = 1e-6
  )
})


test_that("native BigWig writer preserves chromosome-edge intervals", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 10L, 1000L),
    end = c(1L, 12L, 1000L),
    value = c(0, -2.5, 3.25),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )

  outdir <- tempfile(pattern = "gtr_boundary_bw_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = "chr1", size = 1000L),
    overwrite = TRUE
  )
  observed_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy",
    verbose = FALSE
  )
  observed <- retrieve_bwg(observed_track, chrom = "chr1", start = 1L, end = 1000L)

  expect_equal(
    .signal_frame_for_test(observed),
    .signal_frame_for_test(signal),
    tolerance = 1e-4
  )
  expect_equal(
    retrieve_bwg(observed_track, chrom = "chr1", start = 1000L, end = 1000L)$value,
    3.25,
    tolerance = 1e-4
  )
  expect_equal(
    nrow(retrieve_bwg(observed_track, chrom = "chr1", start = 1001L, end = 1100L)),
    0L
  )
})


test_that("strand-specific bedGraph and unstranded BigWig semantics remain distinct", {
  plus_file <- tempfile(fileext = ".bedgraph")
  minus_file <- tempfile(fileext = ".bedgraph")
  writeLines("chr1\t9\t12\t1", plus_file)
  writeLines("chr1\t9\t12\t2", minus_file)

  stranded <- read_bwg(
    c(plus_file, minus_file),
    format = "bedgraph",
    sample_names = c("plus", "minus"),
    strand = c("+", "-"),
    mode = "memory",
    verbose = FALSE
  )
  plus <- retrieve_bwg(stranded, chrom = "chr1", start = 1L, end = 20L, strand = "+")
  minus <- retrieve_bwg(stranded, chrom = "chr1", start = 1L, end = 20L, strand = "-")
  expect_identical(unique(plus$sample_id), "plus")
  expect_identical(unique(minus$sample_id), "minus")

  bigwig <- read_bwg(
    test_path("fixtures", "native_reader_test.bigwig"),
    format = "bigwig",
    sample_names = "bw",
    strand = "+",
    mode = "lazy",
    verbose = FALSE
  )
  expect_identical(bigwig$samples$strand, "*")
  expect_false(bigwig$samples$has_strand)
  expect_gt(nrow(retrieve_bwg(bigwig, chrom = "chr1", start = 1L, end = 200L, strand = "+")), 0L)
  expect_equal(
    nrow(retrieve_bwg(
      bigwig,
      chrom = "chr1",
      start = 1L,
      end = 200L,
      strand = "+",
      strand_policy = "strict"
    )),
    0L
  )
})


test_that("scale normalization is equivalent for lazy and memory BigWig tracks", {
  file <- test_path("fixtures", "native_reader_test.bigwig")
  lazy <- read_bwg(file, format = "bigwig", sample_names = "sampleA", mode = "lazy", verbose = FALSE)
  expect_warning(
    memory <- read_bwg(file, format = "bigwig", sample_names = "sampleA", mode = "memory", verbose = FALSE),
    "Full-memory bigWig loading"
  )

  lazy <- norm_bwg(lazy, method = "scale", scale_factor = c(sampleA = 2))
  memory <- norm_bwg(memory, method = "scale", scale_factor = c(sampleA = 2))
  lazy_region <- retrieve_bwg(lazy, chrom = "chr1", start = 1L, end = 200L)
  memory_region <- retrieve_bwg(memory, chrom = "chr1", start = 1L, end = 200L)

  expect_equal(
    .signal_frame_for_test(lazy_region),
    .signal_frame_for_test(memory_region),
    tolerance = 1e-6
  )
})


test_that("bedGraph WIG and BigWig round trips preserve the same canonical signal", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(1L, 20L, 5L),
    end = c(3L, 22L, 7L),
    value = c(0.5, -2.25, 3.75),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )
  chrom_sizes <- data.frame(chrom = c("chr1", "chr2"), size = c(100L, 100L))

  observed <- list()
  for (format in c("bedgraph", "wig", "bigwig")) {
    outdir <- tempfile(pattern = paste0("gtr_cross_format_", format, "_"))
    dir.create(outdir)
    written <- write_bwg(
      track,
      outdir = outdir,
      format = format,
      chrom_sizes = if (format == "bigwig") chrom_sizes else NULL,
      overwrite = TRUE
    )
    if (format == "bigwig") {
      reread <- read_bwg(written$file, format = format, sample_names = "sampleA", mode = "lazy", verbose = FALSE)
      observed[[format]] <- data.table::rbindlist(list(
        retrieve_bwg(reread, chrom = "chr1", start = 1L, end = 100L),
        retrieve_bwg(reread, chrom = "chr2", start = 1L, end = 100L)
      ))
    } else {
      reread <- read_bwg(written$file, format = format, sample_names = "sampleA", mode = "memory", verbose = FALSE)
      observed[[format]] <- reread$data
    }
  }

  expected <- .signal_frame_for_test(signal)
  expect_equal(.signal_frame_for_test(observed$bedgraph), expected, tolerance = 1e-10)
  expect_equal(.signal_frame_for_test(observed$wig), expected, tolerance = 1e-10)
  expect_equal(.signal_frame_for_test(observed$bigwig), expected, tolerance = 1e-4)
})


test_that("signal summaries and binning do not mutate caller-owned data tables", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(10L, 20L),
    end = c(12L, 22L),
    value = c(1, 2),
    strand = "*"
  )
  track <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA"),
    data = data.table::copy(signal),
    meta = list(mode = "memory")
  )
  original_track_data <- data.table::copy(track$data)
  summary_bwg(track)
  expect_identical(names(track$data), names(original_track_data))
  expect_equal(track$data, original_track_data)

  input <- data.table::copy(signal)
  original_input <- data.table::copy(input)
  bin_bwg(input, bin_size = 10L)
  expect_identical(names(input), names(original_input))
  expect_equal(input, original_input)
})


test_that("merge_bwg keep_first keeps records from the first duplicate sample only", {
  first <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", norm_method = "none"),
    data = data.table::data.table(
      sample_id = "sampleA", chrom = "chr1", start = 10L, end = 12L,
      value = 1, strand = "*"
    ),
    seqinfo = data.table::data.table(sample_id = "sampleA", chrom = "chr1", length = 100L),
    meta = list(mode = "memory")
  )
  second <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", norm_method = "none"),
    data = data.table::data.table(
      sample_id = "sampleA", chrom = "chr1", start = 20L, end = 22L,
      value = 9, strand = "*"
    ),
    seqinfo = data.table::data.table(sample_id = "sampleA", chrom = "chr1", length = 100L),
    meta = list(mode = "memory")
  )

  merged <- merge_bwg(first, second, sample_conflict = "keep_first")
  expect_equal(nrow(merged$samples), 1L)
  expect_equal(nrow(merged$data), 1L)
  expect_equal(merged$data$start, 10L)
  expect_equal(merged$data$value, 1)
  expect_equal(merged$seqinfo$length, 100L)
})
