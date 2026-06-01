test_that("read_bwg reads bedGraph tracks and retrieves regional signal", {
  files <- gtr_extdata(c("example_signal_A.bedgraph", "example_signal_B.bedgraph"))
  bg <- read_bwg(files, format = "bedgraph", mode = "lazy", verbose = FALSE)

  expect_s3_class(bg, "BwgTrack")
  expect_equal(nrow(bg$samples), 2L)
  dt <- retrieve_bwg(bg, chrom = "chr1", start = 101, end = 450)
  expect_gt(nrow(dt), 0L)
  expect_true(all(dt$start <= 450L))
  expect_true(all(dt$end >= 101L))
})

test_that("gene model, signal, and combined track plots are generated", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  bg <- read_bwg(gtr_extdata(c("example_signal_A.bedgraph", "example_signal_B.bedgraph")), format = "bedgraph", mode = "lazy", verbose = FALSE)
  bed <- read_bed(gtr_extdata("example_features.bed"), verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_variants.vcf"), mode = "memory", verbose = FALSE)

  expect_s3_class(plot_gene(gp, gene_id = "GeneA"), "ggplot")
  expect_s3_class(plot_transcript(gp, transcript_id = "TxA1"), "ggplot")
  expect_s3_class(plot_region(gp, chrom = "chr1", start = 1, end = 1200), "ggplot")

  expect_true(inherits(plot_signal_gene(bg, gp, gene_id = "GeneA"), "patchwork"))
  expect_true(inherits(plot_signal_transcript(bg, gp, transcript_id = "TxA1"), "patchwork"))
  expect_true(inherits(plot_signal_region(bg, chrom = "chr1", start = 101, end = 900, annotation = gp), "patchwork"))

  p <- plot_tracks(
    annotation = gp,
    signal = bg,
    features = bed,
    variants = vcf,
    chrom = "chr1",
    start = 1,
    end = 1200,
    gene_color_palette = "Set2"
  )
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("signal plotting supports sample groups and summaries", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  bg <- read_bwg(gtr_extdata(c("example_signal_A.bedgraph", "example_signal_B.bedgraph")), format = "bedgraph", mode = "lazy", verbose = FALSE)
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
