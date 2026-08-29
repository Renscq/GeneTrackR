test_that("GeneTrackR package metadata no longer requires compiled code", {
  desc <- packageDescription("GeneTrackR")
  imports <- if (is.null(desc$Imports)) "" else desc$Imports
  linking_to <- if (is.null(desc$LinkingTo)) "" else desc$LinkingTo
  needs_compilation <- if (is.null(desc$NeedsCompilation)) "no" else desc$NeedsCompilation
  dependency_text <- paste(imports, linking_to)

  expect_false(grepl("Rcpp", dependency_text, fixed = TRUE))
  expect_identical(needs_compilation, "no")
})


test_that("compiled BigWig compatibility helpers are absent from the namespace", {
  ns <- asNamespace("GeneTrackR")
  compiled_helpers <- c(
    "bw_seqinfo_cpp",
    "bw_query_cpp",
    "query_bigwig_cpp",
    "read_bigwig_whole_cpp",
    "write_bigwig_libbigwig"
  )

  expect_false(any(vapply(
    compiled_helpers,
    exists,
    logical(1L),
    envir = ns,
    inherits = FALSE
  )))
})


test_that("pure R BigWig read and write remain functional without compiled helpers", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1"),
    start = c(10L, 30L),
    end = c(12L, 32L),
    value = c(1.5, -2.25),
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

  outdir <- tempfile(pattern = "gtr_pure_r_")
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
  observed <- retrieve_bwg(
    observed_track,
    chrom = "chr1",
    start = 1L,
    end = 1000L
  )

  expect_equal(
    observed[, .(sample_id, chrom, start, end)],
    signal[, .(sample_id, chrom, start, end)]
  )
  expect_equal(observed$value, signal$value, tolerance = 1e-4)
})
