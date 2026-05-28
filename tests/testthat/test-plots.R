test_that("gene model plotting functions return ggplot objects", {
  gp <- read_genepred(
    system.file("extdata", "example.genePredExt", package = "GeneTrackR"),
    format = "genePredExt",
    verbose = FALSE
  )

  expect_s3_class(plot_gene(gp, gene_id = "GeneA"), "ggplot")
  expect_s3_class(plot_transcript(gp, transcript_id = "TxA1"), "ggplot")
  expect_s3_class(plot_region(gp, chrom = "chr1", start = 1, end = 1000), "ggplot")
})

test_that("signal and combined plotting functions return plot objects", {
  gp <- read_genepred(
    system.file("extdata", "example.genePredExt", package = "GeneTrackR"),
    format = "genePredExt",
    verbose = FALSE
  )
  bg <- read_bwg(
    system.file("extdata", c("example_signal_A.bedgraph", "example_signal_B.bedgraph"), package = "GeneTrackR"),
    format = "bedgraph",
    mode = "lazy",
    verbose = FALSE
  )

  expect_s3_class(plot_signal_region(bg, annotation = gp, chrom = "chr1", start = 101, end = 900), "patchwork")
  expect_true(inherits(plot_tracks(annotation = gp, signal = bg, gene_id = "GeneA"), "patchwork"))
  expect_true(inherits(plot_tracks(annotation = gp, signal = bg, transcript_id = "TxA1"), "patchwork"))
  expect_true(inherits(plot_tracks(annotation = gp, signal = bg, chrom = "chr1", start = 101, end = 900), "patchwork"))
})


test_that("signal plotting supports sample groups and replicate summaries", {
  gp <- read_genepred(
    system.file("extdata", "example.genePredExt", package = "GeneTrackR"),
    format = "genePredExt",
    verbose = FALSE
  )
  bg <- read_bwg(
    system.file("extdata", c("example_signal_A.bedgraph", "example_signal_B.bedgraph"), package = "GeneTrackR"),
    format = "bedgraph",
    mode = "lazy",
    verbose = FALSE
  )
  groups <- c(example_signal_A = "A", example_signal_B = "B")

  p <- plot_signal_region(
    bg,
    annotation = gp,
    chrom = "chr1",
    start = 101,
    end = 900,
    sample_groups = groups,
    signal_color_by = "group",
    signal_summary = "mean",
    bin_size = 100
  )
  expect_true(inherits(p, "patchwork"))
})
