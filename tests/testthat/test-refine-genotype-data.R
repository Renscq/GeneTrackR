refine_gt_hap <- function(calls, format = rep("GT", ncol(calls)), alt = rep("G", ncol(calls))) {
  variants <- data.table::data.table(
    chrom = "chr1", pos = seq_len(ncol(calls)), variant_id = colnames(calls),
    ref = "A", alt = alt, FORMAT = format
  )
  for (sample in rownames(calls)) {
    data.table::set(variants, j = sample, value = calls[sample, ])
  }
  vcf <- VariantTrack(variants, meta = list(sample_names = rownames(calls)))
  hap_region_variant(vcf, chrom = "chr1", start = 1, end = ncol(calls), min_variant_number = 0)
}

test_that("the GeneA demo gives exact redundant groups of sizes 9, 1, and 1", {
  gp <- read_genepred(gtr_extdata("gtr_demo.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  hap <- hap_gene_variant(vcf, annotation = gp, gene_id = "GeneA", min_variant_number = 1)

  refined <- refine_genotype_data(hap)

  expect_equal(dim(refined$dosage), c(36L, 11L))
  expect_equal(sort(as.integer(table(refined$variant_map$block_id))), c(1L, 1L, 9L))
  expect_identical(rownames(refined$dosage), names(refined$pattern))
})

test_that("exact reversed columns recover every dosage and propagate protection", {
  calls <- cbind(
    first = c("0/0", "0/0", "0/1", "1/1", "1/1", "./."),
    reverse = c("1/1", "1/1", "0/1", "0/0", "0/0", "./."),
    recombinant = c("0/0", "0/1", "0/1", "1/1", "1/1", "./."),
    other_missing = c("0/0", "0/0", "0/1", "1/1", "1/1", "1/1"),
    constant = rep("0/0", 6),
    constant_reverse = rep("1/1", 6)
  )
  rownames(calls) <- paste0("s", 1:6)
  hap <- refine_gt_hap(calls)

  refined <- refine_genotype_data(hap, protect_variants = "reverse")
  map <- refined$variant_map

  expect_identical(map$block_id, c("rV001_2", "rV001_2", "rV002_1", "rV003_1", "rV004_1", "rV005_1"))
  expect_identical(map$representative_id[1:2], rep("first", 2))
  expect_identical(map$flipped, c(FALSE, TRUE, rep(FALSE, 4)))
  expect_identical(map$protected, c(TRUE, TRUE, rep(FALSE, 4)))
  expect_identical(unname(refined$block_supported), c(TRUE, TRUE, TRUE, FALSE, FALSE))
  restored <- refined$blocks[, map$block_id, drop = FALSE]
  restored[, map$flipped] <- 2 - restored[, map$flipped]
  colnames(restored) <- map$variant_id
  expect_identical(restored, refined$dosage)
})

test_that("raw GT separates heterozygotes, missing calls, and unsupported states", {
  calls <- cbind(v1 = c("5:0/1", "5:1/1", "5:1|0", "5:0/.", "5:2/2", "5:0/0/1"))
  rownames(calls) <- paste0("s", 1:6)
  hap <- refine_gt_hap(calls, format = "DP:GT")

  refined <- refine_genotype_data(hap)

  expect_identical(unname(refined$dosage[, "v1"]), c(1, 2, 1, NA, NA, NA))
  expect_identical(unname(refined$unsupported[, "v1"]), c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(refined$original_hap_id[["s1"]], refined$original_hap_id[["s2"]])
  expect_identical(refined$pattern[["s1"]], refined$pattern[["s3"]])
  expect_length(unique(refined$pattern[c("s1", "s2", "s4", "s5", "s6")]), 5L)
  expect_false(refined$block_supported[["rV001_1"]])
})

test_that("complex alleles retain distinct raw states as unsupported singletons", {
  calls <- cbind(v1 = c("0/0", "0/1", "0/2", "./."), v2 = c("0/0", "0/1", "0/2", "./."))
  rownames(calls) <- paste0("s", 1:4)
  hap <- refine_gt_hap(calls, alt = c("G,T", "G,T"))

  refined <- refine_genotype_data(hap)

  expect_true(all(is.na(refined$dosage)))
  expect_equal(colSums(refined$unsupported), c(v1 = 3, v2 = 3))
  expect_equal(ncol(refined$blocks), 2L)
  expect_length(unique(refined$pattern), 4L)
  matrix_input <- matrix(0, nrow(calls), ncol(calls), dimnames = dimnames(calls))
  expect_error(refine_genotype_data(hap, genotypes = matrix_input), "biallelic")
})

test_that("an explicit dosage matrix aligns fully by both ID axes", {
  calls <- cbind(a = c("0/0", "0/1", "1/1"), b = c("1/1", "0/1", "0/0"))
  rownames(calls) <- c("s1", "s2", "s3")
  hap <- refine_gt_hap(calls)
  raw <- refine_genotype_data(hap)
  dosages <- matrix(c(0, 1, 2, 2, 1, 0), ncol = 2, dimnames = dimnames(calls))
  hap$sample_haplotypes <- hap$sample_haplotypes[3:1]
  hap$variants <- hap$variants[2:1, c("chrom", "pos", "variant_id", "ref", "alt")]

  explicit <- refine_genotype_data(hap, genotypes = dosages[3:1, 2:1])

  expect_identical(explicit$pattern, raw$pattern)
  expect_identical(explicit$dosage, dosages)
  expect_identical(explicit$variant_map, raw$variant_map)
  expect_error(refine_genotype_data(hap), "raw GT")
  expect_error(refine_genotype_data(hap, genotypes = dosages[-1, ]), "sample IDs")
  expect_error(refine_genotype_data(hap, genotypes = dosages[, -1, drop = FALSE]), "variant IDs")
  dosages[1, 1] <- 0.5
  expect_error(refine_genotype_data(hap, genotypes = dosages), "0, 1, 2")
  expect_error(refine_genotype_data(hap, protect_variants = "absent"), "protect_variants")
})
