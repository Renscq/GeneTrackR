test_that("read_bwg reads example bedGraph tracks", {
  files <- system.file(
    "extdata",
    c("example_signal_A.bedgraph", "example_signal_B.bedgraph"),
    package = "GeneTrackR"
  )
  bg <- read_bwg(files, format = "bedgraph", mode = "lazy", verbose = FALSE)

  expect_s3_class(bg, "BwgTrack")
  expect_equal(nrow(bg$samples), 2L)
  expect_equal(bg$samples$sample_id, c("example_signal_A", "example_signal_B"))
})

test_that("query_bwg returns regional bedGraph signal", {
  files <- system.file(
    "extdata",
    c("example_signal_A.bedgraph", "example_signal_B.bedgraph"),
    package = "GeneTrackR"
  )
  bg <- read_bwg(files, format = "bedgraph", mode = "lazy", verbose = FALSE)
  dt <- query_bwg(bg, chrom = "chr1", start = 101, end = 450)

  expect_true(nrow(dt) > 0L)
  expect_true(all(dt$start <= 450L))
  expect_true(all(dt$end >= 101L))
  expect_equal(sort(unique(dt$sample_id)), c("example_signal_A", "example_signal_B"))
})
