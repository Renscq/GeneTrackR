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
    gene_palette = "Set2"
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

test_that("track themes and signal style parameters are applied consistently", {
  dt <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 3L),
    end = c(2L, 4L),
    value = c(4, 25),
    strand = "*"
  )

  p <- plot_signal_core(
    dt,
    plot_type = "bar",
    signal_transform = "sqrt",
    signal_y_scale = "fixed",
    signal_y_limits = c(0, 4),
    signal_alpha = 0.60,
    signal_bar_width = 0.50,
    plot_theme = "classic",
    show_panel_border = FALSE
  )
  built <- ggplot2::ggplot_build(p)
  bar_layer <- built$data[[1L]]

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(bar_layer), 2L)
  expect_equal(bar_layer$xmax - bar_layer$xmin, rep(1, 2))
  expect_equal(unique(bar_layer$alpha), 0.60)
  expect_s3_class(p$theme$panel.border, "element_blank")
  expect_equal(p$scales$get_scales("y")$limits, c(0, 4))

  expect_warning(
    plot_signal_core(
      dt,
      signal_y_scale = "free",
      signal_y_limits = c(0, 10)
    ),
    "changed to 'fixed'"
  )
})

test_that("new plotting parameter validation is strict", {
  expect_error(normalize_signal_alpha(1.1), "between 0 and 1")
  expect_error(normalize_signal_bar_width(0), "greater than 0")
  expect_error(normalize_signal_y_limits(c(2, 1)), "increasing")
  expect_error(normalize_show_panel_border(1), "NULL, TRUE, or FALSE")
  expect_error(normalize_plot_theme("dark"), "arg")
})

test_that("gene model feature colors are stable across observed feature subsets", {
  utr_only <- make_gene_model_fill_colors(features = "UTR")
  cds_utr <- make_gene_model_fill_colors(features = c("CDS", "UTR"))
  exon_only <- make_gene_model_fill_colors(features = "exon")

  expect_identical(gene_model_feature_levels(), c("UTR", "CDS", "exon"))
  expect_true(all(c("UTR", "CDS", "exon") %in% names(utr_only)))
  expect_identical(names(utr_only)[seq_along(gene_model_feature_levels())], gene_model_feature_levels())
  expect_identical(unname(utr_only["UTR"]), unname(cds_utr["UTR"]))
  expect_identical(unname(exon_only["exon"]), unname(cds_utr["exon"]))

  alias_cols <- make_gene_model_fill_colors(features = c("five_prime_UTR", "three_prime_UTR"))
  expect_true("UTR" %in% names(alias_cols))
})


test_that("variant marker colors are stable and ordered", {
  cols <- make_variant_marker_fill_colors(variant_types = c("DEL", "SNP"))
  expect_identical(variant_marker_levels(), c("SNP", "Ind", "..."))
  expect_true(all(c("SNP", "Ind", "...") %in% names(cols)))
  expect_identical(normalize_variant_marker_type(c("SNP", "INS", "unknown")), c("SNP", "Ind", "..."))
})
