test_that("feature and VCF writers create standard output files", {
  gp <- read_genepred(gtr_extdata("gtr_demo.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)

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


test_that("bundled libBigWig backend writes and rereads bigWig", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(10L, 20L, 100L),
    end = c(12L, 22L, 101L),
    value = c(1.5, 2.5, 3.5),
    strand = "+"
  )
  bg <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", strand = "+"),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )
  chrom_sizes <- data.table::data.table(
    chrom = c("chr1", "chr2"),
    size = c(1000L, 2000L)
  )

  outdir <- tempfile(pattern = "gtr_bigwig_")
  dir.create(outdir)
  written <- write_bwg(
    bg,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes,
    overwrite = TRUE
  )

  expect_true(file.exists(written$file))
  expect_equal(written$format, "bigwig")

  expect_warning(
    reread <- read_bwg(
      written$file,
      format = "bigwig",
      sample_names = "sampleA",
      strand = "+",
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )

  observed <- reread$data[, .(chrom, start, end, value)]
  expected <- signal[, .(chrom, start, end, value)]
  data.table::setorderv(observed, c("chrom", "start", "end"))
  data.table::setorderv(expected, c("chrom", "start", "end"))
  expect_equal(observed[, .(chrom, start, end)], expected[, .(chrom, start, end)])
  expect_equal(observed$value, expected$value, tolerance = 1e-6)
})


test_that("bigWig export rejects overlaps before calling libBigWig", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(10L, 12L),
    end = c(15L, 20L),
    value = c(1, 2),
    strand = "+"
  )
  bg <- BwgTrack(
    samples = data.table::data.table(sample_id = "sampleA", strand = "+"),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )

  expect_error(
    write_bwg(
      bg,
      outdir = tempfile(pattern = "gtr_bigwig_overlap_"),
      format = "bigwig",
      chrom_sizes = data.frame(chrom = "chr1", size = 1000),
      overwrite = TRUE
    ),
    "overlapping intervals"
  )
})


test_that("bundled libBigWig backend writes multiple strand-specific samples", {
  signal <- data.table::data.table(
    sample_id = c("RNA_seq_plus", "RNA_seq_plus", "RNA_seq_minus", "RNA_seq_minus"),
    chrom = c("chr1", "chr1", "chr1", "chr1"),
    start = c(10L, 20L, 40L, 50L),
    end = c(12L, 22L, 42L, 52L),
    value = c(10, 20, 15, 25),
    strand = c("+", "+", "-", "-")
  )
  bg <- BwgTrack(
    samples = data.table::data.table(
      sample_id = c("RNA_seq_plus", "RNA_seq_minus"),
      strand = c("+", "-")
    ),
    data = signal,
    meta = list(mode = "memory", coordinate = "1-based closed")
  )

  outdir <- tempfile(pattern = "gtr_bigwig_multi_")
  dir.create(outdir)
  written <- write_bwg(
    bg,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = "chr1", size = 1000),
    overwrite = TRUE
  )

  expect_equal(written$sample_id, c("RNA_seq_plus", "RNA_seq_minus"))
  expect_true(all(file.exists(written$file)))
  expect_true(all(written$format == "bigwig"))

  expect_warning(
    plus <- read_bwg(
      written$file[written$sample_id == "RNA_seq_plus"],
      format = "bigwig",
      sample_names = "RNA_seq_plus",
      strand = "+",
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )
  expect_warning(
    minus <- read_bwg(
      written$file[written$sample_id == "RNA_seq_minus"],
      format = "bigwig",
      sample_names = "RNA_seq_minus",
      strand = "-",
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )

  expect_equal(
    plus$data[, .(chrom, start, end, value)],
    signal[sample_id == "RNA_seq_plus", .(chrom, start, end, value)],
    tolerance = 1e-6
  )
  expect_equal(
    minus$data[, .(chrom, start, end, value)],
    signal[sample_id == "RNA_seq_minus", .(chrom, start, end, value)],
    tolerance = 1e-6
  )
})


test_that("demo strand-specific RNA-seq tracks round-trip through bundled libBigWig", {
  rnaseq <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph"
    )),
    sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
    strand = c("+", "-"),
    mode = "memory",
    verbose = FALSE
  )

  chrom_sizes <- gtr_extdata("gtr_demo.chrom.sizes")
  outdir <- tempfile(pattern = "gtr_demo_bigwig_")
  dir.create(outdir)

  written <- write_bwg(
    rnaseq,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes,
    overwrite = TRUE
  )

  expect_equal(written$sample_id, c("RNA_seq_plus", "RNA_seq_minus"))
  expect_true(all(file.exists(written$file)))
  expect_true(all(file.info(written$file)$size > 0))

  expect_warning(
    reread <- read_bwg(
      written$file,
      format = "bigwig",
      sample_names = written$sample_id,
      strand = c("+", "-"),
      mode = "memory",
      verbose = FALSE
    ),
    "Full-memory bigWig loading"
  )

  observed <- data.table::copy(reread$data)
  expected <- data.table::copy(rnaseq$data)
  data.table::setorderv(observed, c("sample_id", "chrom", "start", "end"))
  data.table::setorderv(expected, c("sample_id", "chrom", "start", "end"))

  expect_equal(
    observed[, .(sample_id, chrom, start, end)],
    expected[, .(sample_id, chrom, start, end)]
  )
  expect_equal(observed$value, expected$value, tolerance = 1e-6)
})


test_that("write_bwg exposes a single bundled bigWig backend", {
  args <- names(formals(write_bwg))
  expect_false("bigwig_backend" %in% args)
  expect_false("bedGraphToBigWig" %in% args)
})
