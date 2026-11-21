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

ld_test_variants <- function() {
  VariantTrack(data.frame(
    chrom = "chr1",
    pos = c(40L, 20L, 10L, 30L),
    variant_id = c("sparse", "second", "first", "constant"),
    ref = c("A", "A", "AT", "A"),
    alt = c("G", "AT", "A", "G"),
    S1 = c("./.", "0/1", "0/0", "0/0"),
    S2 = c("0/1", "0/1", "0/0", "0/0"),
    S3 = c("1/1", "0/0", "0/1", "0/0"),
    S4 = c(".", "0|1", "0|1", "0|0"),
    S5 = c("0/.", "./.", "1/1", "0/0"),
    S6 = c("0/0", "1/1", "1/1", "0/0")
  ))
}

test_that("LDTrack sorts variant metadata and supports the legacy figure alias", {
  figure <- ggplot2::ggplot()
  ld <- LDTrack(
    data = data.frame(index_i = 1L, index_j = 2L, ld = 0.75),
    variants = data.frame(chrom = "chr1", pos = c(20L, 10L), variant_id = c("b", "a")),
    region = list(chrom = "chr1", start = 10L, end = 20L),
    meta = list(method = "r2", sample_n = 4L),
    plot = figure
  )

  expect_s3_class(ld$data, "data.table")
  expect_identical(ld$variants$variant_id, c("a", "b"))
  expect_identical(ld$variants$variant_index, 1:2)
  expect_identical(ld$figure, figure)
  expect_output(print(ld), "region    : chr1:10-20", fixed = TRUE)
  expect_error(
    LDTrack(variants = data.frame(chrom = "chr1", pos = 10L)),
    "`variants` must contain `chrom`, `pos`, and `variant_id` columns.",
    fixed = TRUE
  )
})

test_that("LD calculation aligns named samples and distinguishes uninformative pairs", {
  samples <- c("S4", "S2", "S5", "S1", "S3")
  ld <- compute_ld_block(
    ld_test_variants(), method = "Dprime", samples = samples, min_pair_samples = 4L,
    keep_genotype_matrix = TRUE, verbose = FALSE
  )

  ids <- c("first", "second", "constant", "sparse")
  expected_dosage <- rbind(
    c(1, 0, 2, 0, 1),
    c(1, 1, NA, 1, 0),
    c(0, 0, 0, 0, 0),
    c(NA, 1, NA, NA, 2)
  )
  dimnames(expected_dosage) <- list(ids, samples)
  expect_identical(ld$region, list(chrom = "chr1", start = 10L, end = 40L))
  expect_identical(ld$variants$variant_id, ids)
  expect_identical(ld$meta$sample_n, 5L)
  expect_identical(ld$genotype, expected_dosage)
  expect_identical(ld$data$variant_i, ids[ld$data$index_i])
  expect_identical(ld$data$variant_j, ids[ld$data$index_j])
  expect_identical(ld$data$n_samples, c(4L, 5L, 2L, 4L, 2L, 2L))

  # The first pair has dosage covariance -1/8 and Dmax of 3/32.
  missing <- rep(NA_real_, 5L)
  expect_equal(ld$data$r2, c(1 / 3, missing), tolerance = 1e-12)
  expect_equal(ld$data$D, c(-1 / 16, missing), tolerance = 1e-12)
  expect_equal(ld$data$Dprime_signed, c(-2 / 3, missing), tolerance = 1e-12)
  expect_equal(ld$data$ld, c(2 / 3, missing), tolerance = 1e-12)
  expect_identical(ld$data$ld, ld$data$Dprime)

  expected_matrix <- matrix(NA_real_, nrow = 4L, ncol = 4L, dimnames = list(ids, ids))
  diag(expected_matrix) <- 1
  expected_matrix[1L, 2L] <- expected_matrix[2L, 1L] <- 2 / 3
  expect_equal(ld$matrix, expected_matrix, tolerance = 1e-12)
})

test_that("LD calculation explains invalid samples, regions, and genotype-free input", {
  vcf <- ld_test_variants()
  no_genotypes <- VariantTrack(data.frame(chrom = "chr1", pos = c(10L, 20L)))

  expect_error(
    compute_ld_block(vcf, samples = "unknown", verbose = FALSE),
    "Samples not found in VCF: unknown", fixed = TRUE
  )
  expect_error(
    compute_ld_block(vcf, chrom = "chrMissing", verbose = FALSE),
    "No variants were found on the requested chromosome", fixed = TRUE
  )
  expect_error(
    compute_ld_block(vcf, chrom = "chr1", start = 10L, end = 10L, verbose = FALSE),
    "At least two variants are required for LD calculation", fixed = TRUE
  )
  expect_error(
    compute_ld_block(no_genotypes, verbose = FALSE),
    "No VCF sample genotype columns were found", fixed = TRUE
  )
})

test_that("LD plotting uses position and Reds fallbacks with optional layers disabled", {
  expect_warning(
    figure <- plot_ld_block(
      ld_test_variants(), variant_type = "ind", method = "Dprime",
      samples = c("S4", "S2", "S5", "S1", "S3"), min_pair_samples = 4L,
      label_by = "unknown", color_palette = "invalid-palette",
      show_region = TRUE, connect_region = FALSE, variant_marker_size = 0,
      return_object = FALSE, verbose = FALSE
    ),
    "Unknown `color_palette`: invalid-palette. Falling back to 'Reds'.", fixed = TRUE
  )

  expect_s3_class(figure, "patchwork")
  expect_length(figure[[1L]]$layers, 0L)
  variants <- attr(figure, "variant_data")
  expect_identical(variants$variant_id, c("first", "second"))
  expect_identical(variants$variant_label, c("10", "20"))
  expect_equal(attr(figure, "ld_data")$ld, 2 / 3, tolerance = 1e-12)

  heatmap <- figure[[2L]]
  expect_identical(heatmap$labels$title, "chr1:10-40 LD (Dprime)")
  fill_scale <- heatmap$scales$get_scales("fill")
  expect_identical(fill_scale$name, "D'")
  expect_identical(fill_scale$palette(c(0, 1)), RColorBrewer::brewer.pal(9L, "Reds")[c(1L, 9L)])
})
