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
