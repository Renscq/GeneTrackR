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


test_that("WIG export preserves interval spans", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(10L, 20L, 30L, 40L),
    end = c(12L, 22L, 31L, 42L),
    value = c(1.5, 2.5, 3.5, 4.5),
    strand = "*"
  )
  bg <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA"),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )

  outdir <- tempfile(pattern = "gtr_wig_")
  dir.create(outdir)
  written <- write_bwg(bg, outdir = outdir, format = "wig", overwrite = TRUE)
  reread <- read_bwg(written$file, format = "wig", mode = "memory", verbose = FALSE)

  observed <- reread$data[, c("chrom", "start", "end", "value"), with = FALSE]
  expected <- signal[, c("chrom", "start", "end", "value"), with = FALSE]
  expect_equal(observed, expected)
})

test_that("lazy signal copies preserve compression suffixes and validate source files", {
  src <- tempfile(fileext = ".bedgraph.gz")
  con <- gzfile(src, open = "wt")
  writeLines(c("chr1\t0\t10\t1", "chr1\t10\t20\t2"), con)
  close(con)

  bg <- read_bwg(src, format = "auto", mode = "lazy", use_tabix = "no", verbose = FALSE)
  outdir <- tempfile(pattern = "gtr_lazy_copy_")
  dir.create(outdir)
  written <- write_bwg(bg, outdir = outdir, format = "bedgraph", overwrite = TRUE)

  expect_true(grepl("\\.bedgraph\\.gz$", written$file))
  expect_true(file.exists(written$file))

  unlink(src)
  expect_error(
    write_bwg(bg, outdir = tempfile(pattern = "gtr_missing_source_"), format = "bedgraph"),
    "Source signal file does not exist"
  )
})
