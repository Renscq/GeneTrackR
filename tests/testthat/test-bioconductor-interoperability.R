test_that("core genomic track objects convert to GRanges", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE,
    progress = FALSE
  )
  signal <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph"
    )),
    format = "bedgraph",
    sample_names = c("RNA_plus", "RNA_minus"),
    strand = c("+", "-"),
    mode = "memory",
    verbose = FALSE
  )

  gene_gr <- as_granges(gp, level = "gene")
  signal_gr <- as_granges(signal)
  variant_gr <- as_granges(vcf)

  expect_true(inherits(gene_gr, "GRanges"))
  expect_true(inherits(signal_gr, "GRanges"))
  expect_true(inherits(variant_gr, "GRanges"))

  signal_meta <- as.data.frame(S4Vectors::mcols(signal_gr))
  expect_true(all(c("sample_id", "score") %in% names(signal_meta)))
  expect_true(all(as.character(GenomicRanges::strand(signal_gr)) %in% c("+", "-", "*")))

  variant_meta <- as.data.frame(S4Vectors::mcols(variant_gr))
  expect_true(all(c(
    "pos", "variant_id", "ref", "alt", "variant_type"
  ) %in% names(variant_meta)))
  expect_true(all(as.character(GenomicRanges::strand(variant_gr)) == "*"))
})

test_that("lazy track conversion requires regional retrieval", {
  lazy_signal <- BwgTrack(
    samples = data.table::data.table(
      sample_id = "sampleA",
      file = "sampleA.bigwig",
      format = "bigwig",
      strand = "*"
    ),
    data = NULL,
    meta = list(mode = "lazy")
  )
  lazy_variant <- VariantTrack(
    data = NULL,
    meta = list(
      lazy = TRUE,
      source_file = "variants.vcf.gz",
      format = "VCF"
    )
  )

  expect_error(
    as_granges(lazy_signal),
    "retrieve_bwg"
  )
  expect_error(
    as_granges(lazy_variant),
    "retrieve_vcf"
  )
})

test_that("retrieve_vcf can return GRanges", {
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE,
    progress = FALSE
  )

  gr <- retrieve_vcf(
    vcf,
    chrom = "chr1",
    start = 12339700,
    end = 12343200,
    as = "GRanges",
    verbose = FALSE
  )
  dt <- retrieve_vcf(
    vcf,
    chrom = "chr1",
    start = 12339700,
    end = 12343200,
    as = "data.table",
    verbose = FALSE
  )

  expect_true(inherits(gr, "GRanges"))
  expect_identical(as.character(GenomicRanges::seqnames(gr)), as.character(dt$chrom))
  expect_identical(as.integer(IRanges::start(gr)), as.integer(dt$pos))
  expect_identical(nrow(as.data.frame(S4Vectors::mcols(gr))), nrow(dt))
})

test_that("dependency placement matches runtime usage", {
  desc <- utils::packageDescription("GeneTrackR")
  imports <- trimws(strsplit(desc[["Imports"]], ",", fixed = TRUE)[[1L]])
  suggests <- trimws(strsplit(desc[["Suggests"]], ",", fixed = TRUE)[[1L]])

  expect_true(all(c(
    "GenomicRanges", "IRanges", "S4Vectors", "grDevices", "tools"
  ) %in% imports))
  expect_true("Rsamtools" %in% suggests)
  expect_false("Rsamtools" %in% imports)
  expect_identical(intersect(imports, suggests), character())
})
