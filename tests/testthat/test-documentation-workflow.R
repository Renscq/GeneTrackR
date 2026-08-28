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
    end = 12355500,
    variant_type = "snp",
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

test_that("browser-track documentation workflow uses current public APIs", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  signal_all <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph",
      "gtr_demo_riboseq_plus.bedgraph",
      "gtr_demo_riboseq_minus.bedgraph"
    )),
    format = "bedgraph",
    sample_names = c(
      "RNA_seq_plus",
      "RNA_seq_minus",
      "Ribo_seq_plus",
      "Ribo_seq_minus"
    ),
    strand = c("+", "-", "+", "-"),
    mode = "memory",
    verbose = FALSE
  )
  features <- read_bed(
    gtr_extdata("gtr_demo_features.bed"),
    verbose = FALSE,
    progress = FALSE
  )
  variants <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  p_gene <- plot_tracks(
    annotation = gp,
    signal = signal_all,
    gene_id = "GeneA",
    samples = c("RNA_seq_plus", "Ribo_seq_plus"),
    strand = "+",
    signal_type = "bar"
  )
  expect_true(inherits(p_gene, "patchwork") || inherits(p_gene, "ggplot"))

  p_transcript <- plot_tracks(
    annotation = gp,
    signal = signal_all,
    transcript_id = "TxA1",
    samples = c("RNA_seq_plus", "Ribo_seq_plus"),
    strand = "+",
    signal_type = "bar",
    ribo_signal_type = "auto"
  )
  expect_true(inherits(p_transcript, "patchwork") || inherits(p_transcript, "ggplot"))

  p_complete <- plot_tracks(
    annotation = gp,
    signal = signal_all,
    features = features,
    variants = variants,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    samples = c("RNA_seq_plus", "Ribo_seq_plus"),
    strand = "+",
    signal_type = "bar",
    variant_palette = "Paired",
    variant_colors = c(SNP = "#1F78B4", INS = "#33A02C", DEL = "#E31A1C"),
    highlight = data.frame(start = 12342620, end = 12343180),
    layout = "gene_top",
    heights = c(signal = 4, gene = 1.5, feature = 0.9, variant = 0.8)
  )
  expect_true(inherits(p_complete, "patchwork") || inherits(p_complete, "ggplot"))
})

test_that("haplotype documentation workflow preserves the designed GeneA truth", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  hap_gene <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "string"
  )

  expect_s3_class(hap_gene, "HapVariant")
  expect_equal(hap_gene$meta$variant_n, 11L)
  expect_equal(hap_gene$meta$sample_n, 36L)
  expect_equal(hap_gene$meta$haplotype_n, 4L)
  expect_equal(sort(hap_gene$haplotypes$sample_n), rep(9L, 4L))
  expect_true(all(c(
    "region", "variants", "genotype_long", "genotype_wide",
    "haplotypes", "sample_haplotypes", "meta"
  ) %in% names(hap_gene)))

  hap_tx <- hap_gene_variant(
    vcf,
    annotation = gp,
    transcript_id = "TxA1",
    genotype_mode = "string"
  )
  expect_s3_class(hap_tx, "HapVariant")
  expect_identical(hap_tx$region$locator, "transcript")

  hap_region <- hap_region_variant(
    vcf,
    chrom = "chr1",
    start = 12340001,
    end = 12352000,
    genotype_mode = "code"
  )
  expect_s3_class(hap_region, "HapVariant")
  expect_equal(hap_region$meta$variant_n, 11L)

  hap_relaxed <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    upstream = 500,
    genotype_mode = "string",
    min_variant_number = 11
  )
  expect_true("varAup01" %in% hap_relaxed$variants$variant_id)
  expect_equal(hap_relaxed$meta$sample_n, 36L)

  p <- plot_hap_variant(
    hap_gene,
    annotation = gp,
    min_hap_samples = 3,
    variant_label = "pos",
    gene_track_legend_position = "top"
  )
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("phenotype documentation workflow preserves sample matching and result contracts", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE,
    progress = FALSE
  )
  pheno <- read_pheno(
    gtr_extdata("gtr_demo_pheno.tsv"),
    verbose = FALSE,
    progress = FALSE
  )
  hap_gene <- hap_gene_variant(
    vcf = vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "string"
  )

  expect_equal(nrow(pheno), 36L)
  expect_true(setequal(
    as.character(pheno$sample_id),
    as.character(hap_gene$sample_haplotypes$sample_id)
  ))

  pheno_summary <- summary_pheno(pheno)
  expect_identical(
    pheno_summary[trait == "flower_color", type],
    "categorical"
  )
  expect_true(all(
    pheno_summary[trait %in% c(
      "seed_weight",
      "protein_content",
      "plant_height",
      "flowering_time"
    ), type] == "numeric"
  ))

  hap_seed <- plot_hap_pheno(
    hap = hap_gene,
    phenotype = pheno,
    traits = "seed_weight",
    min_hap_samples = 3,
    p_adjust = "BH",
    p_value_type = "adjusted"
  )
  expect_s3_class(hap_seed, "GeneTrackRPhenoPlot")
  expect_true(all(c(
    "figure", "pvalue", "summary", "bracket", "plot_data"
  ) %in% names(hap_seed)))
  expect_true(all(c("p_value", "p_adj") %in% names(hap_seed$pvalue)))
  expect_equal(length(unique(as.character(hap_seed$summary$hap_id))), 4L)

  variant_protein <- plot_variant_pheno(
    variant = vcf,
    phenotype = pheno,
    variant_id = "varA03",
    traits = "protein_content",
    genotype_mode = "code",
    min_group_samples = 3,
    p_adjust = "BH",
    p_value_type = "adjusted"
  )
  expect_s3_class(variant_protein, "GeneTrackRPhenoPlot")
  expect_true("variant_data" %in% names(variant_protein))
  expect_identical(
    as.character(variant_protein$variant_data$variant_id[1L]),
    "varA03"
  )
})

test_that("haplotype refinement documentation workflow has deterministic demo outcomes", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  vcf <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    keep_genotype = TRUE,
    verbose = FALSE,
    progress = FALSE
  )
  pheno <- read_pheno(
    gtr_extdata("gtr_demo_pheno.tsv"),
    verbose = FALSE,
    progress = FALSE
  )
  hap_gene <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "string"
  )

  refined_protein <- refine_haplotype(
    hap_gene,
    phenotype = pheno,
    traits = "protein_content",
    min_hap_samples = 3,
    test_method = "t.test",
    p_adjust = "BH",
    alpha = 0.05,
    effect_threshold = 0.5
  )
  expect_s3_class(refined_protein, "HapRefined")
  expect_equal(nrow(refined_protein$refined_haplotypes), 2L)
  expect_equal(sort(refined_protein$refined_haplotypes$sample_n), c(18L, 18L))
  expect_equal(sum(refined_protein$pairwise_test$can_merge, na.rm = TRUE), 2L)

  grouped_members <- refined_protein$haplotype_map[, .(
    original_haps = paste(sort(hap_id), collapse = ";")
  ), by = refined_hap_id]
  expect_setequal(
    grouped_members$original_haps,
    c("Hap1;Hap3", "Hap2;Hap4")
  )
  expect_setequal(
    as.character(refined_protein$refined_haplotypes$varA03),
    c("C", "G")
  )

  grouped_plot <- plot_refined_hap_variant(
    refined_protein,
    annotation = gp,
    min_hap_samples = 3,
    collapse_refined = FALSE
  )
  collapsed_plot <- plot_refined_hap_variant(
    refined_protein,
    annotation = gp,
    min_hap_samples = 3,
    collapse_refined = TRUE
  )
  expect_true(inherits(grouped_plot, "patchwork") || inherits(grouped_plot, "ggplot"))
  expect_true(inherits(collapsed_plot, "patchwork") || inherits(collapsed_plot, "ggplot"))

  refined_pheno <- plot_refined_hap_pheno(
    refined_protein,
    phenotype = pheno,
    traits = "protein_content",
    min_hap_samples = 3,
    p_value_type = "adjusted"
  )
  expect_s3_class(refined_pheno, "GeneTrackRPhenoPlot")

  refined_seed_weight <- refine_haplotype(
    hap_gene,
    phenotype = pheno,
    traits = "seed_weight",
    min_hap_samples = 3,
    effect_threshold = 0.5
  )
  expect_equal(nrow(refined_seed_weight$refined_haplotypes), 4L)

  refined_flowering <- refine_haplotype(
    hap_gene,
    phenotype = pheno,
    traits = "flowering_time",
    min_hap_samples = 3,
    effect_threshold = 0.5
  )
  expect_equal(nrow(refined_flowering$refined_haplotypes), 1L)

  refined_multi_trait <- refine_haplotype(
    hap_gene,
    phenotype = pheno,
    traits = c("protein_content", "seed_weight"),
    min_hap_samples = 3,
    effect_threshold = 0.5
  )
  expect_equal(nrow(refined_multi_trait$refined_haplotypes), 4L)
})
