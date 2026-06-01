test_that("VCF reader loads genotype-rich example VCF", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)

  expect_s3_class(vcf, "VariantTrack")
  expect_gte(nrow(vcf$data), 500L)
  sample_names <- if (is.null(vcf$meta$sample_names)) character() else vcf$meta$sample_names
  expect_gte(length(sample_names), 60L)
  expect_true(all(c("chrom", "pos", "variant_id", "ref", "alt") %in% names(vcf$data)))
  expect_true("S01" %in% sample_names)
})

test_that("retrieve_vcf subsets in-memory VariantTrack by chromosome and range", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  sub <- retrieve_vcf(vcf, chrom = "chr1", start = 1, end = 1200, as = "VariantTrack", verbose = FALSE)

  expect_s3_class(sub, "VariantTrack")
  expect_gt(nrow(sub$data), 0L)
  expect_true(all(sub$data$chrom == "chr1"))
  expect_true(all(sub$data$pos >= 1L & sub$data$pos <= 1200L))
})

test_that("read_vcf lazy mode falls back safely when no tabix index is available", {
  expect_warning(
    vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "lazy", verbose = FALSE),
    "requires a bgzip-compressed VCF"
  )
  expect_s3_class(vcf, "VariantTrack")
  expect_gt(nrow(vcf$data), 0L)
})

test_that("variant plotting returns ggplot objects", {
  vcf <- read_vcf(gtr_extdata("example_variants.vcf"), mode = "memory", verbose = FALSE)
  p <- plot_variant(vcf, chrom = "chr1", start = 1, end = 1200)
  expect_s3_class(p, "ggplot")
})
