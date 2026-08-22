test_that("new deterministic demo data has expected core dimensions", {
  gp <- read_genepred(gtr_extdata("gtr_demo.genePredExt"), format = "genePredExt", verbose = FALSE)
  gtf <- read_gtf(gtr_extdata("gtr_demo.gtf"), verbose = FALSE, progress = FALSE)
  gff <- read_gff(gtr_extdata("gtr_demo.gff3"), verbose = FALSE, progress = FALSE)
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)
  pheno <- read_pheno(gtr_extdata("gtr_demo_pheno.tsv"), verbose = FALSE)

  expect_equal(nrow(gp$genes), 20L)
  expect_equal(nrow(gp$transcripts), 24L)
  expect_equal(nrow(gtf$genes), 20L)
  expect_equal(nrow(gtf$transcripts), 24L)
  expect_equal(nrow(gff$genes), 20L)
  expect_equal(nrow(gff$transcripts), 24L)
  expect_equal(nrow(vcf$data), 50L)
  expect_equal(length(vcf$meta$sample_names), 36L)
  expect_equal(nrow(pheno), 36L)
  expect_setequal(pheno$sample_id, vcf$meta$sample_names)
  expect_false(identical(pheno$sample_id, vcf$meta$sample_names))
})


test_that("GeneA demo locus contains four balanced designed haplotypes", {
  gp <- read_genepred(gtr_extdata("gtr_demo.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)

  hap <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "code"
  )

  expect_equal(nrow(hap$variants), 11L)
  expect_equal(nrow(hap$haplotypes), 4L)
  expect_equal(sort(hap$haplotypes$sample_n), rep(9L, 4L))
})


test_that("demo LD truth includes perfect high-LD block and two-variant case", {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)

  ld_high <- compute_ld_block(
    vcf,
    chrom = "chr1",
    start = 12342620L,
    end = 12343180L,
    method = "r2",
    verbose = FALSE
  )
  expect_equal(nrow(ld_high$variants), 6L)
  expect_true(all(abs(ld_high$data$r2 - 1) < 1e-12, na.rm = TRUE))
  expect_false(anyNA(ld_high$data$r2))

  ld_pair <- compute_ld_block(
    vcf,
    chrom = "chr2",
    start = 16995001L,
    end = 17006000L,
    method = "r2",
    verbose = FALSE
  )
  expect_equal(nrow(ld_pair$variants), 2L)
  expect_equal(nrow(ld_pair$data), 1L)
  expect_equal(ld_pair$variants$variant_id, c("varPair01", "varPair02"))
})


test_that("demo phenotype encodes designed effects and a negative control", {
  model_file <- system.file("scripts", "demo_model", "samples.tsv", package = "GeneTrackR", mustWork = TRUE)
  samples <- data.table::fread(model_file, data.table = TRUE)
  pheno <- read_pheno(gtr_extdata("gtr_demo_pheno.tsv"), verbose = FALSE)
  x <- merge(samples, pheno, by = "sample_id", sort = FALSE)

  seed_means <- x[, mean(seed_weight), by = hap_group][order(hap_group), V1]
  protein_means <- x[, mean(protein_content), by = hap_group][order(hap_group), V1]
  flowering_means <- x[, mean(flowering_time), by = hap_group][order(hap_group), V1]

  expect_equal(round(seed_means, 6), c(20, 25, 31, 23))
  expect_equal(round(protein_means, 6), c(38, 44, 44, 38))
  expect_equal(round(flowering_means, 6), rep(45, 4L))
})


test_that("demo RNA-seq signal is exon-enriched and intron-depleted", {
  rnaseq <- data.table::fread(
    gtr_extdata("gtr_demo_rnaseq.bedgraph"),
    col.names = c("chrom", "start", "end", "value")
  )

  expect_true(all(rnaseq$start < rnaseq$end))
  expect_true(all(rnaseq$value > 0))

  genea_exon <- rnaseq[chrom == "chr1" & start < 12341000L & end > 12340000L]
  genea_intron <- rnaseq[chrom == "chr1" & start < 12342000L & end > 12341000L]
  genei_lncRNA <- rnaseq[chrom == "chr1" & start < 15006000L & end > 15000000L]

  expect_gt(nrow(genea_exon), 0L)
  expect_equal(nrow(genea_intron), 0L)
  expect_gt(nrow(genei_lncRNA), 0L)
})


test_that("demo Ribo-seq signal is single-base CDS density with 3-nt periodicity", {
  riboseq <- data.table::fread(
    gtr_extdata("gtr_demo_riboseq.bedgraph"),
    col.names = c("chrom", "start", "end", "value")
  )

  expect_true(all(riboseq$end - riboseq$start == 1L))
  expect_true(all(riboseq$value > 0))

  # GeneI is a designed lncRNA and therefore has no Ribo-seq CDS signal.
  expect_equal(
    nrow(riboseq[chrom == "chr1" & start < 15006000L & end > 15000000L]),
    0L
  )

  genea <- riboseq[chrom == "chr1" & start >= 12340500L & end <= 12351500L]
  genea[, pos := start + 1L]
  start_peak <- genea[pos == 12340501L, value]
  stop_peak <- genea[pos >= 12351498L & pos <= 12351500L, max(value)]
  periodic <- genea[pos >= 12340531L & pos <= 12340990L]
  periodic[, phase := (pos - 12340501L) %% 3L]
  phase_means <- periodic[, mean(value), by = phase][order(phase)]$V1

  expect_equal(length(start_peak), 1L)
  expect_gt(start_peak, max(periodic$value))
  expect_gt(stop_peak, max(periodic$value))
  expect_true(phase_means[1L] > phase_means[2L])
  expect_true(phase_means[2L] > phase_means[3L])
})
