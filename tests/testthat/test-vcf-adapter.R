# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Regress the optional typed VCF interoperation contract.
# Input: A synthetic VCF with duplicate IDs, multiple alleles and mixed ploidy.
# Output: Assertions for standard and native variant workflows.

adapter_fixture <- function() {
  if (identical(Sys.getenv("GENETRACKR_REQUIRE_VARIANTANNOTATION"), "true")) {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
      stop("VariantAnnotation is required for the VCF adapter test suite.")
    }
  } else {
    testthat::skip_if_not_installed("VariantAnnotation")
  }
  VariantAnnotation::readVcf(
    testthat::test_path("fixtures", "standard-variants.vcf"), genome = NA_character_
  )
}

test_that("standard VCF objects round-trip with their typed fields intact", {
  source <- adapter_fixture()
  original <- serialize(source, NULL)
  track <- as_variant_track(source)

  expect_s3_class(track, "VariantTrack")
  expect_identical(track$data$pos, c(1L, 5L, 5L, 20L, 40L))
  expect_identical(track$data$start, track$data$end)
  expect_identical(track$data$alt, c("C,G", "A", "AG", "T", "<DEL>"))
  expect_identical(track$meta$sample_names, c("S2", "S1"))
  expect_identical(sub(":.*$", "", track$data$S2[c(1L, 3L, 4L)]), c("0|2", "./1", "1"))
  expect_identical(serialize(source, NULL), original)

  restored <- as_vcf(track)
  expect_s4_class(restored, "CollapsedVCF")
  expect_identical(VariantAnnotation::header(restored), VariantAnnotation::header(source))
  expect_identical(VariantAnnotation::fixed(restored), VariantAnnotation::fixed(source))
  expect_identical(VariantAnnotation::info(restored), VariantAnnotation::info(source))
  expect_identical(VariantAnnotation::geno(restored), VariantAnnotation::geno(source))
  expect_identical(rownames(restored), rownames(source))
  expect_identical(VariantAnnotation::info(restored)$END[[5L]], 50L)
  expect_true(VariantAnnotation::info(restored)$FLAG[[1L]])
})

test_that("region, row and sample selection maintain typed VCF associations", {
  source <- adapter_fixture()
  track <- as_variant_track(source[c(5L, 3L, 2L, 1L, 4L), ])
  selected <- retrieve_vcf(track, chrom = "chr1", start = 5, end = 5,
    as = "VariantTrack", verbose = FALSE)
  selected$data <- selected$data[c(2L, 1L)]
  checkpoint <- tempfile(fileext = ".rds")
  on.exit(unlink(checkpoint), add = TRUE)
  saveRDS(selected, checkpoint)

  result <- as_vcf(readRDS(checkpoint), samples = c("S1", "S2"))
  expected <- source[c(2L, 3L), c("S1", "S2")]
  expect_identical(selected$data$variant_id, rep("dup", 2L))
  expect_identical(colnames(result), c("S1", "S2"))
  expect_identical(VariantAnnotation::samples(VariantAnnotation::header(result)), colnames(result))
  expect_identical(VariantAnnotation::fixed(result), VariantAnnotation::fixed(expected))
  expect_identical(VariantAnnotation::info(result), VariantAnnotation::info(expected))
  expect_identical(VariantAnnotation::geno(result), VariantAnnotation::geno(expected))
  expect_identical(dim(as_vcf(selected, samples = "S2")), c(2L, 1L))

  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  VariantAnnotation::writeVcf(result, path)
  reread <- VariantAnnotation::readVcf(path, genome = NA_character_)
  expect_identical(dim(reread), dim(result))
  expect_identical(colnames(reread), colnames(result))
  original_header <- VariantAnnotation::header(VariantAnnotation::header(result))
  written_header <- VariantAnnotation::header(VariantAnnotation::header(reread))
  expect_identical(written_header[names(original_header)], original_header)
  expect_identical(VariantAnnotation::info(reread), VariantAnnotation::info(result))
  expect_identical(VariantAnnotation::geno(reread), VariantAnnotation::geno(result))
})

test_that("object conversion preserves high-precision genotype data and file dates", {
  adapter_fixture()
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2", "##fileDate=20210430",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t7\tv1\tA\tC\t10\tPASS\t.\tGT:GP\t1|0:0.123456789,0.345678912,0.530864299\t1:0.123456789,0.876543211"), path)
  source <- VariantAnnotation::readVcf(path, genome = NA_character_)
  selected <- as_vcf(as_variant_track(source), samples = c("S2", "S1"))
  expected <- source[, c("S2", "S1")]
  expect_identical(VariantAnnotation::geno(selected)$GP, VariantAnnotation::geno(expected)$GP)
  expect_identical(VariantAnnotation::geno(selected)$GT, VariantAnnotation::geno(expected)$GT)
  expect_identical(VariantAnnotation::meta(VariantAnnotation::header(selected))$fileDate,
    VariantAnnotation::meta(VariantAnnotation::header(source))$fileDate)
})

test_that("empty and sites-only standard variants have no phantom samples", {
  source <- adapter_fixture()
  track <- as_variant_track(source)
  empty <- as_vcf(retrieve_vcf(track, chrom = "chr1", start = 80, end = 90,
    as = "VariantTrack", verbose = FALSE))
  expect_identical(dim(empty), c(0L, 2L))
  expect_identical(colnames(empty), colnames(source))
  expect_identical(dim(as_vcf(as_variant_track(source[integer(), ]))), c(0L, 2L))

  sites <- as_vcf(track, samples = character())
  expect_identical(dim(sites), c(5L, 0L))
  expect_length(VariantAnnotation::samples(VariantAnnotation::header(sites)), 0L)
  expect_length(VariantAnnotation::geno(sites), 0L)
  sites_track <- as_variant_track(sites)
  expect_length(get_vcf_sample_columns(sites_track), 0L)
  expect_identical(VariantAnnotation::info(as_vcf(sites_track)), VariantAnnotation::info(sites))

  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  VariantAnnotation::writeVcf(sites, path)
  expect_identical(dim(VariantAnnotation::readVcf(path, genome = NA_character_)), c(5L, 0L))
  expect_false(any(grepl("^#CHROM.*FORMAT", readLines(path))))
})

test_that("unrecoverable variant edits and ambiguous samples fail explicitly", {
  source <- adapter_fixture()
  track <- as_variant_track(source)
  legacy <- read_vcf(testthat::test_path("fixtures", "standard-variants.vcf"), verbose = FALSE)
  expect_error(as_variant_track(VariantAnnotation::expand(source)), "CollapsedVCF")
  expect_error(as_vcf(legacy), "provenance")
  expect_error(as_vcf(merge_vcf(track, track)), "provenance")
  expect_error(as_vcf(track, samples = "unknown"), "Unknown VCF sample")
  expect_error(as_vcf(track, samples = c("S1", "S1")), "unique sample names")
  expect_error(as_vcf(track, samples = 1L), "unique sample names")

  edited <- track
  edited$data <- data.table::copy(edited$data)
  edited$data[1L, "alt" := "T"]
  expect_error(as_vcf(edited), "edited.*alt")
  edited$data <- data.table::copy(track$data)
  edited$data[1L, "S1" := "0/0"]
  expect_error(as_vcf(edited), "edited.*S1")
  edited$data <- track$data[c(1L, 1L)]
  expect_error(as_vcf(edited), "provenance")
  other <- as_variant_track(source)
  edited$data <- data.table::rbindlist(list(track$data[1L], other$data[2L]))
  expect_error(as_vcf(edited), "provenance")

  colnames(source) <- c("ID", "S1")
  expect_error(as_variant_track(source), "reserved.*column")

  incomplete <- adapter_fixture()
  info <- VariantAnnotation::info(incomplete)
  info$UNDECLARED <- seq_len(nrow(incomplete))
  suppressWarnings(VariantAnnotation::info(incomplete) <- info)
  expect_error(as_variant_track(incomplete), "declare every projected INFO and FORMAT field")
})

test_that("malformed standard VCF headers and ambiguous ALT elements are rejected", {
  source <- adapter_fixture()
  missing_format <- source
  header <- VariantAnnotation::header(source)
  fields <- VariantAnnotation::header(header)
  VariantAnnotation::header(missing_format) <- VariantAnnotation::VCFHeader(
    reference = VariantAnnotation::reference(header), samples = VariantAnnotation::samples(header),
    header = fields[names(fields) != "fileformat"]
  )
  expect_error(as_variant_track(missing_format), "fileformat")

  alternate <- VariantAnnotation::alt(source)
  alternate[[1L]] <- c("C", "")
  VariantAnnotation::alt(source) <- alternate
  expect_error(as_variant_track(source), "empty or missing ALT allele")
})

test_that("VCF metadata containing #CHROM is distinct from the column header", {
  adapter_fixture()
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth of a #CHROM record">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t7\tv1\tA\tC\t10\tPASS\tDP=3"), path)
  source <- VariantAnnotation::readVcf(path, genome = NA_character_)
  expect_identical(as_variant_track(source)$data$pos, 7L)
})

test_that("indexed standard VCF subsets retain range IDs and genotype arrays", {
  adapter_fixture()
  input <- testthat::test_path("fixtures", "standard-variants.vcf")
  indexed <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(indexed, paste0(indexed, ".tbi"))), add = TRUE)
  Rsamtools::bgzip(input, dest = indexed)
  Rsamtools::indexTabix(indexed, format = "vcf")
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 5L))
  source <- VariantAnnotation::readVcf(indexed, genome = NA_character_,
    param = VariantAnnotation::ScanVcfParam(which = region))
  track <- as_variant_track(source)
  native <- read_vcf(indexed, chrom = "chr1", start = 1L, end = 5L, verbose = FALSE)
  expect_identical(track$data$pos, native$data$pos)
  expect_identical(track$data$pos, c(1L, 5L, 5L))
  expect_identical(GenomicRanges::granges(as_vcf(track)), GenomicRanges::granges(source))
  expect_identical(VariantAnnotation::geno(as_vcf(track)), VariantAnnotation::geno(source))
})

test_that("the adapter leaves default VCF export and native analysis compatible", {
  source <- adapter_fixture()
  track <- as_variant_track(source)
  file <- testthat::test_path("fixtures", "standard-variants.vcf")
  native <- read_vcf(file, verbose = FALSE)
  expect_identical(get_vcf_sample_columns(track), get_vcf_sample_columns(native))

  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  write_vcf(track, path)
  columns <- strsplit(grep("^#CHROM", readLines(path), value = TRUE), "\t", fixed = TRUE)[[1L]]
  expect_length(columns, 8L)
  expect_false(any(grepl("gtr_vcf_row", names(S4Vectors::mcols(as_granges(track))))))

  native_hap <- hap_region_variant(native, "chr2", 20L, 40L)
  adapted_hap <- hap_region_variant(track, "chr2", 20L, 40L)
  expect_equal(adapted_hap$genotype_long, native_hap$genotype_long)
  expect_equal(adapted_hap$haplotypes, native_hap$haplotypes)
  expect_equal(adapted_hap$sample_haplotypes, native_hap$sample_haplotypes)
})
