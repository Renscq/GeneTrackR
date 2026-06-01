test_that("hap_gene_variant and hap_region_variant build HapVariant objects", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)

  hap_gene <- hap_gene_variant(vcf, annotation = gp, gene_id = "GeneA", genotype_mode = "string", min_variant_number = 1)
  hap_region <- hap_region_variant(vcf, chrom = "chr1", start = 1, end = 1200, genotype_mode = "code", min_variant_number = 1)

  expect_s3_class(hap_gene, "HapVariant")
  expect_s3_class(hap_region, "HapVariant")
  expect_gt(nrow(hap_gene$variants), 0L)
  expect_gt(nrow(hap_gene$haplotypes), 0L)
  expect_gt(nrow(hap_region$sample_haplotypes), 0L)
})

test_that("plot_hap_variant returns a browser-like haplotype figure", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  hap <- hap_gene_variant(vcf, annotation = gp, gene_id = "GeneA", genotype_mode = "string", min_variant_number = 1)

  p <- plot_hap_variant(
    hap,
    annotation = gp,
    min_hap_samples = 1,
    show_reference_row = TRUE,
    show_gene_pos_axis = TRUE,
    gene_pos_x_angle = 90,
    gene_track_legend_position = "top"
  )
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("phenotype readers and haplotype phenotype plots are stable", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  pheno <- read_pheno(gtr_extdata("example_pheno.tsv"), verbose = FALSE)
  hap <- hap_gene_variant(vcf, annotation = gp, gene_id = "GeneA", genotype_mode = "code", min_variant_number = 1)

  pheno_sum <- summary_pheno(pheno)
  expect_true(all(c("trait", "type", "missing_n") %in% names(pheno_sum)))
  expect_s3_class(plot_pheno(pheno, traits = c("plant_height", "seed_weight")), "ggplot")

  res <- plot_hap_pheno(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
  expect_s3_class(res, "GeneTrackRPhenoPlot")
  expect_true(inherits(res$figure, "ggplot"))
  expect_true(is.data.frame(res$pvalue))
  expect_true(is.data.frame(res$summary))
})

test_that("single variant phenotype plots return figure and p-value tables", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  pheno <- read_pheno(gtr_extdata("example_pheno.tsv"), verbose = FALSE)

  res <- plot_variant_pheno(
    variant = vcf,
    phenotype = pheno,
    variant_id = "rsA1",
    traits = "plant_height",
    min_group_samples = 1
  )

  expect_s3_class(res, "GeneTrackRPhenoPlot")
  expect_true(inherits(res$figure, "ggplot"))
  expect_true(is.data.frame(res$pvalue))
  expect_true(is.data.frame(res$variant_data))
})
