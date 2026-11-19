coverage_vcf_file <- function() {
  file <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100>",
    "##contig=<ID=chr2,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1",
    "chr1\t5\t.\tA\tG\t10\tPASS\t.\tGT\t0|1\t1/1",
    "chr1\t10\tv2\tA\tAT\t.\tPASS\t.\tGT\t0/0\t./.",
    "chr2\t8\tv3\tC\tT\t20\tPASS\t.\tGT\t1/1\t0/1"
  ), file)
  file
}

test_that("plain VCF queries retain raw genotypes and explicit point coordinates", {
  file <- coverage_vcf_file()
  on.exit(unlink(file), add = TRUE)
  observed <- read_vcf(file, mode = "memory", verbose = FALSE, progress = FALSE)
  expect_identical(observed$meta$sample_names, c("S2", "S1"))
  expect_identical(observed$data$variant_id, c("chr1:5:A:G", "v2", "v3"))
  expect_identical(observed$data$S2, c("0|1", "0/0", "1/1"))
  expect_identical(observed$data$S1, c("1/1", "./.", "0/1"))
  expect_identical(observed$data$start, c(5L, 10L, 8L))
  expect_identical(observed$data$end, observed$data$pos)
  without <- read_vcf(file, keep_genotype = FALSE, mode = "memory", verbose = FALSE)
  expect_identical(without$meta$sample_names, character())
  expect_false(without$meta$has_format)
  expect_false(any(c("FORMAT", "S1", "S2") %in% names(without$data)))
  alias <- read_vcf_track(file, mode = "memory", verbose = FALSE)
  expect_equal(as.data.frame(alias$data), as.data.frame(observed$data))
  expect_identical(alias$meta, observed$meta)
  regional <- read_vcf(file, chrom = "chr1", start = 5L, end = 5L, verbose = FALSE)
  expect_identical(regional$data$pos, 5L)
  expect_identical(retrieve_vcf(file, chrom = "chr2", verbose = FALSE)$variant_id, "v3")
  ranges <- retrieve_vcf(observed, variant_id = "v2", as = "GRanges", verbose = FALSE)
  expect_identical(as.integer(IRanges::width(ranges)), 1L)
  expect_identical(S4Vectors::mcols(ranges)$alt, "AT")
  expect_equal(nrow(summary_vcf(observed)), 3L)
  expect_equal(nrow(summary_vcf(observed, by = "not_a_column")), 2L)
  expect_equal(nrow(summary_vcf(observed, chrom = "chr9")), 0L)
})

test_that("VCF logging and invalid location diagnostics follow documented modes", {
  file <- coverage_vcf_file()
  on.exit(unlink(file), add = TRUE)
  messages <- character()
  withCallingHandlers(read_vcf(file, verbose = TRUE, progress = TRUE),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_match(paste(messages, collapse = "\n"), "Progress")
  expect_match(paste(messages, collapse = "\n"), "samples: 2")
  expect_error(read_vcf(file, chrom = "chr1", start = 1L, verbose = FALSE),
                "must all be supplied")
  object <- read_vcf(file, mode = "memory", verbose = FALSE)
  expect_error(retrieve_vcf(object, start = 1L, end = 2L, verbose = FALSE), "chrom")
  expect_error(retrieve_vcf(object, chrom = c("chr1", "chr2"), start = 1L,
                             end = 2L, verbose = FALSE), "exactly one chromosome")
  expect_error(retrieve_vcf(object, chrom = "", verbose = FALSE), "non-empty")
  expect_error(retrieve_vcf(object, upstream = -1L, verbose = FALSE), "non-negative")
  expect_error(retrieve_vcf(object, gene_id = "plus_gene", chrom = "chr1",
                             annotation = coverage_gene_model(), verbose = FALSE), "Do not mix")
  located <- retrieve_vcf(object, gene_id = "plus_gene", annotation = coverage_gene_model(),
                           verbose = FALSE)
  expect_identical(located$pos, c(5L, 10L))
  bad <- tempfile(fileext = ".vcf")
  on.exit(unlink(bad), add = TRUE)
  writeLines("##fileformat=VCFv4.2", bad)
  expect_error(read_vcf(bad, verbose = FALSE), "#CHROM")
  expect_error(read_vcf(paste0(file, ".missing"), verbose = FALSE), "does not exist")
})

test_that("existing VCF export is explicitly sites-only and leaves input unchanged", {
  input <- coverage_vcf_file()
  output <- tempfile(fileext = ".vcf")
  on.exit(unlink(c(input, output)), add = TRUE)
  object <- read_vcf(input, mode = "memory", verbose = FALSE)
  before <- serialize(object, NULL)
  write_vcf(object, output, overwrite = TRUE)
  header <- readLines(output, warn = FALSE)
  header <- header[startsWith(header, "#CHROM")]
  expect_identical(strsplit(header, "\t", fixed = TRUE)[[1L]],
                   c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
  observed <- read_vcf(output, mode = "memory", verbose = FALSE)
  expect_identical(observed$meta$sample_names, character())
  expect_false(any(c("FORMAT", "S1", "S2") %in% names(observed$data)))
  expect_identical(observed$data$variant_id, c("chr1:5:A:G", "v2", "v3"))
  write_variant_track(object, output, overwrite = TRUE, include_header = FALSE)
  expect_true(startsWith(readLines(output, n = 1L), "#CHROM"))
  expect_identical(serialize(object, NULL), before)
})

test_that("real indexed VCF retrieval preserves sample schema in nonempty and empty regions", {
  skip_if_not_installed("Rsamtools")
  file <- coverage_vcf_file()
  compressed <- paste0(file, ".bgz")
  on.exit(unlink(c(file, compressed, paste0(compressed, ".tbi"))), add = TRUE)
  Rsamtools::bgzip(file, dest = compressed, overwrite = TRUE)
  Rsamtools::indexTabix(compressed, format = "vcf")
  messages <- character()
  lazy <- withCallingHandlers(read_vcf(compressed, verbose = TRUE, progress = TRUE),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_null(lazy$data)
  expect_true(lazy$meta$lazy)
  expect_true(lazy$meta$indexed)
  expect_identical(lazy$meta$sample_names, c("S2", "S1"))
  expect_match(paste(messages, collapse = "\n"), "Created lazy VariantTrack")
  expect_match(paste(capture.output(print(lazy)), collapse = "\n"), "lazy")
  selected <- retrieve_vcf(lazy, chrom = "chr1", start = 1L, end = 9L,
                            as = "VariantTrack", verbose = FALSE)
  expect_identical(selected$data$pos, 5L)
  expect_identical(selected$data$S2, "0|1")
  expect_false(selected$meta$lazy)
  messages <- character()
  direct <- withCallingHandlers(
    read_vcf(compressed, chrom = "chr1", start = 1L, end = 9L,
               keep_genotype = FALSE, verbose = TRUE, progress = TRUE),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_identical(direct$data$pos, 5L)
  expect_match(paste(messages, collapse = "\n"), "Scanning tabix index")
  expect_false(any(c("FORMAT", "S1", "S2") %in% names(direct$data)))
  empty <- retrieve_vcf(lazy, chrom = "chr1", start = 50L, end = 60L,
                         as = "VariantTrack", verbose = FALSE)
  expect_equal(nrow(empty$data), 0L)
  expect_identical(empty$meta$sample_names, c("S2", "S1"))
  expect_true(all(c("FORMAT", "S1", "S2") %in% names(empty$data)))
  expect_equal(nrow(validate_vcf(empty)$invalid_records), 0L)
  observed <- summary(lazy, chrom = "chr1", start = 1L, end = 9L)
  expect_identical(observed$n_variants, 1L)
})
