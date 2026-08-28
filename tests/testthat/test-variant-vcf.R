test_that("VCF reader loads genotype-rich example VCF", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)

  expect_s3_class(vcf, "VariantTrack")
  expect_equal(nrow(vcf$data), 50L)
  sample_names <- if (is.null(vcf$meta$sample_names)) character() else vcf$meta$sample_names
  expect_equal(length(sample_names), 36L)
  expect_true(all(c("chrom", "pos", "variant_id", "ref", "alt") %in% names(vcf$data)))
  expect_true("S01" %in% sample_names)
})

test_that("retrieve_vcf subsets in-memory VariantTrack by chromosome and range", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  sub <- retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200, as = "VariantTrack", verbose = FALSE)

  expect_s3_class(sub, "VariantTrack")
  expect_gt(nrow(sub$data), 0L)
  expect_true(all(sub$data$chrom == "chr1"))
  expect_true(all(sub$data$pos >= 12339700L & sub$data$pos <= 12343200L))
})

test_that("read_vcf lazy mode falls back safely when no tabix index is available", {
  expect_warning(
    vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "lazy", verbose = FALSE),
    "requires a bgzip-compressed VCF"
  )
  expect_s3_class(vcf, "VariantTrack")
  expect_gt(nrow(vcf$data), 0L)
})

test_that("variant plotting returns ggplot objects", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  p <- plot_variant(vcf, chrom = "chr1", start = 12339700, end = 12343200)
  expect_s3_class(p, "ggplot")
})



test_that("retrieve_vcf supports whole-object and non-regional filters", {
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE
  )

  all_variants <- retrieve_vcf(vcf, verbose = FALSE)
  expect_equal(nrow(all_variants), nrow(vcf$data))

  snps <- retrieve_vcf(vcf, variant_type = "SNP", verbose = FALSE)
  expect_gt(nrow(snps), 0L)
  expect_true(all(snps$variant_type == "SNP"))

  selected <- retrieve_vcf(
    vcf,
    variant_id = c("varA03", "varA04"),
    verbose = FALSE
  )
  expect_setequal(selected$variant_id, c("varA03", "varA04"))

  expect_warning(
    patterned <- retrieve_vcf(
      vcf,
      pattern = "high_ld",
      fixed = TRUE,
      verbose = FALSE
    ),
    NA
  )
  expect_equal(nrow(patterned), 7L)

  expect_warning(
    patterned_upper <- retrieve_vcf(
      vcf,
      pattern = "HIGH_LD",
      fixed = TRUE,
      ignore_case = TRUE,
      verbose = FALSE
    ),
    NA
  )
  expect_setequal(patterned_upper$variant_id, patterned$variant_id)

  patterned_case_sensitive <- retrieve_vcf(
    vcf,
    pattern = "HIGH_LD",
    fixed = TRUE,
    ignore_case = FALSE,
    verbose = FALSE
  )
  expect_equal(nrow(patterned_case_sensitive), 0L)

  chr1 <- retrieve_vcf(vcf, chrom = "chr1", verbose = FALSE)
  expect_gt(nrow(chr1), 0L)
  expect_true(all(chr1$chrom == "chr1"))

  expect_error(
    retrieve_vcf(vcf, chrom = "chr1", start = 12340000L, verbose = FALSE),
    "both `start` and `end`"
  )
  expect_error(
    retrieve_vcf(vcf, upstream = 1000L, verbose = FALSE),
    "require `gene_id` or `transcript_id`"
  )
})


test_that("summary_vcf summarizes the full VariantTrack without a region", {
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE
  )

  whole <- summary_vcf(vcf)
  expect_gt(nrow(whole), 0L)
  expect_equal(sum(whole$n_variants), nrow(vcf$data))

  chr1 <- summary_vcf(vcf, chrom = "chr1")
  expect_gt(nrow(chr1), 0L)
  expect_equal(
    sum(chr1$n_variants),
    sum(vcf$data$chrom == "chr1")
  )

  expect_equal(summary(vcf), whole)
})
