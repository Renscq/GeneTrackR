test_that("documentation workflow uses current public object contracts", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  gtf <- read_gtf(
    gtr_extdata("gtr_demo.gtf"),
    verbose = FALSE,
    progress = FALSE
  )
  gff <- read_gff(
    gtr_extdata("gtr_demo.gff3"),
    verbose = FALSE,
    progress = FALSE
  )
  features <- read_bed(
    gtr_extdata("gtr_demo_features.bed"),
    verbose = FALSE,
    progress = FALSE
  )
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE,
    progress = FALSE
  )
  pheno <- read_pheno(
    gtr_extdata("gtr_demo_pheno.tsv"),
    verbose = FALSE,
    progress = FALSE
  )

  expect_s3_class(gp, "GenePred")
  expect_s3_class(gtf, "FeatureTrack")
  expect_s3_class(gff, "FeatureTrack")
  expect_s3_class(features, "FeatureTrack")
  expect_s3_class(vcf, "VariantTrack")
  expect_s3_class(pheno, "data.table")

  genea <- retrieve_feature(gp, gene_id = "GeneA")
  genea_variant_table <- retrieve_vcf(
    vcf,
    chrom = "chr1",
    start = 12339700,
    end = 12352000,
    verbose = FALSE,
    progress = FALSE
  )
  genea_variants <- retrieve_vcf(
    vcf,
    chrom = "chr1",
    start = 12339700,
    end = 12352000,
    as = "VariantTrack",
    verbose = FALSE,
    progress = FALSE
  )
  expect_s3_class(genea, "Feature")
  expect_s3_class(genea, "GenePred")
  expect_s3_class(genea_variant_table, "data.table")
  expect_s3_class(genea_variants, "VariantTrack")

  gene_figure <- plot_gene(
    gp,
    gene_id = "GeneA",
    gene_palette = "Set2",
    gene_border_color = "black"
  )
  expect_s3_class(gene_figure, "ggplot")

  hap <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "string",
    min_variant_number = 1
  )
  expect_s3_class(hap, "HapVariant")
  expect_true(all(c(
    "region", "variants", "genotype_long", "genotype_wide",
    "haplotypes", "sample_haplotypes", "meta"
  ) %in% names(hap)))

  hap_figure <- plot_hap_variant(
    hap,
    annotation = gp,
    min_hap_samples = 3
  )
  expect_true(inherits(hap_figure, "patchwork") || inherits(hap_figure, "ggplot"))

  hap_pheno <- plot_hap_pheno(
    hap,
    phenotype = pheno,
    traits = "seed_weight",
    min_hap_samples = 3
  )
  expect_s3_class(hap_pheno, "GeneTrackRPhenoPlot")
  expect_true(all(c("figure", "pvalue", "summary", "bracket", "plot_data") %in% names(hap_pheno)))

  variant_pheno <- plot_variant_pheno(
    vcf,
    phenotype = pheno,
    variant_id = "varA03",
    traits = "protein_content",
    genotype_mode = "code",
    min_group_samples = 3
  )
  expect_s3_class(variant_pheno, "GeneTrackRPhenoPlot")
  expect_true(all(c(
    "figure", "pvalue", "summary", "bracket", "plot_data", "variant_data"
  ) %in% names(variant_pheno)))

  ld <- compute_ld_block(
    vcf,
    chrom = "chr1",
    start = 12342620,
    end = 12343180,
    method = "r2",
    verbose = FALSE
  )
  expect_s3_class(ld, "LDTrack")
  expect_true(all(c("data", "matrix", "variants", "region", "genotype", "figure", "meta") %in% names(ld)))

  ld <- plot_ld_block(
    ld,
    annotation = gp,
    show_region = TRUE,
    show_variant_labels = FALSE,
    verbose = FALSE
  )
  expect_s3_class(ld, "LDTrack")
  expect_false(is.null(ld$figure))

  refined <- refine_haplotype(
    hap,
    phenotype = pheno,
    traits = "seed_weight",
    min_hap_samples = 3
  )
  expect_s3_class(refined, "HapRefined")
  expect_true(all(c(
    "original_hap", "refined_hap", "refined_haplotypes",
    "sample_refined_haplotypes", "haplotype_map", "trait_group_map",
    "phenotype_summary", "pairwise_test", "plot_data", "parameters"
  ) %in% names(refined)))
  expect_s3_class(refined$refined_hap, "HapVariant")

  refined_variant_figure <- plot_refined_hap_variant(
    refined,
    annotation = gp,
    min_hap_samples = 3
  )
  expect_true(inherits(refined_variant_figure, "patchwork") || inherits(refined_variant_figure, "ggplot"))

  refined_pheno <- plot_refined_hap_pheno(
    refined,
    phenotype = pheno,
    traits = "seed_weight",
    min_hap_samples = 3
  )
  expect_s3_class(refined_pheno, "GeneTrackRPhenoPlot")
  expect_true(all(c("figure", "pvalue", "summary", "bracket", "plot_data") %in% names(refined_pheno)))

  effect <- plot_variant_effect(
    hap,
    phenotype = pheno,
    traits = "protein_content",
    min_group_samples = 3,
    x_axis = "position"
  )
  expect_s3_class(effect, "GeneTrackRVariantEffectPlot")
  expect_true(all(c("figure", "effect", "plot_data") %in% names(effect)))
})
