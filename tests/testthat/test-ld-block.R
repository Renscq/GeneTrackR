test_that("compute_ld_block builds LDTrack from genotype-rich VCF", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  ld <- compute_ld_block(vcf, chrom = "chr1", start = 1, end = 1200, min_pair_samples = 3, verbose = FALSE)

  expect_s3_class(ld, "LDTrack")
  expect_gt(nrow(ld$variants), 1L)
  expect_gt(nrow(ld$data), 0L)
  expect_true(all(c("variant_i", "variant_j", "r2", "Dprime", "ld") %in% names(ld$data)))
  expect_equal(nrow(ld$matrix), nrow(ld$variants))
  expect_equal(ncol(ld$matrix), nrow(ld$variants))
})

test_that("plot_ld_block returns LDTrack with triangular LD figure", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  ld <- compute_ld_block(vcf, chrom = "chr1", start = 1, end = 1200, min_pair_samples = 3, verbose = FALSE)
  ld_with_figure <- plot_ld_block(ld, show_variant_labels = FALSE)
  expect_s3_class(ld_with_figure, "LDTrack")
  expect_true(inherits(ld_with_figure$figure, "ggplot") || inherits(ld_with_figure$figure, "patchwork"))
  expect_null(ld_with_figure$plot)

  p <- plot_ld_block(ld, show_region = TRUE, show_variant_labels = FALSE, return_object = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})
