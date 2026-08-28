test_that("compute_ld_block builds LDTrack from genotype-rich VCF", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  ld <- compute_ld_block(vcf, chrom = "chr1", start = 12342620, end = 12355500, variant_type = "snp", min_pair_samples = 3, verbose = FALSE)

  expect_s3_class(ld, "LDTrack")
  expect_gt(nrow(ld$variants), 1L)
  expect_gt(nrow(ld$data), 0L)
  expect_true(all(c("variant_i", "variant_j", "r2", "Dprime", "ld") %in% names(ld$data)))
  expect_equal(nrow(ld$matrix), nrow(ld$variants))
  expect_equal(ncol(ld$matrix), nrow(ld$variants))
})

test_that("plot_ld_block returns LDTrack with triangular LD figure", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  ld <- compute_ld_block(vcf, chrom = "chr1", start = 12342620, end = 12355500, variant_type = "snp", min_pair_samples = 3, verbose = FALSE)
  ld_with_figure <- plot_ld_block(ld, show_variant_labels = FALSE)
  expect_s3_class(ld_with_figure, "LDTrack")
  expect_true(inherits(ld_with_figure$figure, "ggplot") || inherits(ld_with_figure$figure, "patchwork"))
  expect_null(ld_with_figure$plot)

  p <- plot_ld_block(ld, show_region = TRUE, show_variant_labels = FALSE, return_object = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})


test_that("two-variant LD geometry is a single diamond-shaped heatmap cell", {
  pair_dt <- data.table::data.table(
    index_i = 1L,
    index_j = 2L,
    ld = 0.75
  )

  poly <- GeneTrackR:::make_ld_triangle_polygons(pair_dt, n_var = 2L)
  frame <- GeneTrackR:::make_ld_triangle_frame(2L)

  expect_equal(poly$x, c(1.5, 2, 1.5, 1))
  expect_equal(poly$y, c(0, 0.5, 1, 0.5))
  expect_equal(poly$ld, rep(0.75, 4L))
  expect_equal(frame$x, c(1.5, 2, 1.5, 1, 1.5))
  expect_equal(frame$y, c(0, 0.5, 1, 0.5, 0))
})

test_that("demo LD truth sets distinguish gradient, low, and two-variant regions", {
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE
  )

  high <- compute_ld_block(
    vcf,
    chrom = "chr1",
    start = 12342620,
    end = 12355500,
    variant_type = "snp",
    method = "r2",
    verbose = FALSE
  )
  expect_equal(nrow(high$variants), 12L)
  expect_equal(nrow(high$data), 66L)
  expect_equal(sum(high$data$r2 >= 0.8, na.rm = TRUE), 15L)
  expect_equal(sum(high$data$r2 < 0.2, na.rm = TRUE), 16L)
  expect_equal(min(high$data$r2, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(high$meta$sample_n, 36L)

  low <- compute_ld_block(
    vcf,
    chrom = "chr2",
    start = 1999000,
    end = 2008500,
    method = "r2",
    verbose = FALSE
  )
  expect_equal(nrow(low$variants), 4L)
  expect_equal(nrow(low$data), 6L)
  expect_equal(max(low$data$r2, na.rm = TRUE), 1 / 9, tolerance = 1e-12)
  expect_equal(sum(low$data$r2 < 1e-12, na.rm = TRUE), 5L)

  pair <- compute_ld_block(
    vcf,
    chrom = "chr2",
    start = 16995001,
    end = 17006000,
    method = "r2",
    verbose = FALSE
  )
  expect_equal(nrow(pair$variants), 2L)
  expect_equal(nrow(pair$data), 1L)
  expect_equal(pair$data$r2, 1, tolerance = 1e-12)
})

test_that("LD workflow can retain the genotype dosage matrix", {
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE
  )
  ld <- compute_ld_block(
    vcf,
    chrom = "chr1",
    start = 12342620,
    end = 12355500,
    variant_type = "snp",
    keep_genotype_matrix = TRUE,
    verbose = FALSE
  )

  expect_true(is.matrix(ld$genotype))
  expect_equal(dim(ld$genotype), c(12L, 36L))
  expect_identical(rownames(ld$genotype), as.character(ld$variants$variant_id))
})


test_that("plot_ld_block defaults to the Reds palette", {
  expect_identical(formals(plot_ld_block)$color_palette, "Reds")
  expect_identical(formals(GeneTrackR:::draw_ld_triangle_heatmap)$color_palette, "Reds")
  expect_identical(formals(GeneTrackR:::make_ld_continuous_palette)$color_palette, "Reds")
})
