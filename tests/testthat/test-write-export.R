test_that("feature and VCF writers create standard output files", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_variants.vcf"), mode = "memory", verbose = FALSE)

  gtf_file <- tempfile(fileext = ".gtf")
  bed6_file <- tempfile(fileext = ".bed")
  bed12_file <- tempfile(fileext = ".bed")
  vcf_file <- tempfile(fileext = ".vcf")

  write_feature(gp, gtf_file, format = "gtf", overwrite = TRUE)
  write_feature(gp, bed6_file, format = "bed6", overwrite = TRUE)
  write_feature(gp, bed12_file, format = "bed12", overwrite = TRUE)
  write_vcf(vcf, vcf_file, overwrite = TRUE)

  expect_true(file.exists(gtf_file))
  expect_true(file.exists(bed6_file))
  expect_true(file.exists(bed12_file))
  expect_true(file.exists(vcf_file))

  bed6_first <- strsplit(readLines(bed6_file, n = 1L), "\t", fixed = TRUE)[[1L]]
  bed12_first <- strsplit(readLines(bed12_file, n = 1L), "\t", fixed = TRUE)[[1L]]
  expect_equal(length(bed6_first), 6L)
  expect_equal(length(bed12_first), 12L)
})
