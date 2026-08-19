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


test_that("custom missing genotype labels remain missing for haplotype filtering", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)

  hap <- hap_region_variant(
    vcf,
    chrom = "chr1",
    start = 1,
    end = 1200,
    genotype_mode = "code",
    missing_genotype = "-",
    min_variant_number = 1
  )

  s10 <- hap$sample_haplotypes[hap$sample_haplotypes$sample_id == "S10", ]
  missing_gt <- hap$genotype_long[
    hap$genotype_long$sample_id == "S10" & hap$genotype_long$variant_id == "rs001_3",
  ]
  expect_equal(s10$non_missing_variant_n, 4L)
  expect_identical(missing_gt$genotype, "-")
  expect_true(missing_gt$genotype_missing)

  hap_complete <- hap_region_variant(
    vcf,
    chrom = "chr1",
    start = 1,
    end = 1200,
    genotype_mode = "code",
    missing_genotype = "-"
  )
  expect_false("S10" %in% hap_complete$sample_haplotypes$sample_id)
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

test_that("compact indel labels remain aligned with REF and ALT rows", {
  expect_equal(
    format_hap_allele(c("A", "AT", "GTTACA", "i2", "i6", "i10")),
    c("A", "i2", "i6", "i2", "i6", "i10")
  )

  labels <- format_hap_table_genotype_label(
    genotype = c("i2", "i6", "i10"),
    ref = c("AT", "A", "A"),
    alt = c("A", "GTTACA", "AAAAAAAAAA"),
    genotype_mode = "string"
  )
  expect_equal(labels, c("i2", "i6", "i10"))
})

test_that("haplotype allele fill classes use a stable five-color mapping", {
  expect_identical(
    normalize_hap_table_fill_class(
      c("A", "T", "C", "G", "i2", "i5", "i10", "NA"),
      genotype_mode = "string"
    ),
    c("A", "T", "C", "G", "indel", "indel", "indel", NA_character_)
  )

  colors <- make_hap_table_fill_colors(
    table_palette = "Paired",
    table_alpha = 1
  )
  expect_identical(names(colors), c("A", "T", "C", "G", "indel"))
  expect_length(unique(unname(colors)), 5L)

  custom_colors <- make_hap_table_fill_colors(
    table_colors = c(
      A = "#111111",
      T = "#222222",
      C = "#333333",
      G = "#444444",
      indel = "#555555"
    ),
    table_alpha = 1
  )
  expect_identical(
    unname(custom_colors),
    c("#111111FF", "#222222FF", "#333333FF", "#444444FF", "#555555FF")
  )
})

test_that("haplotype gene-track variant legend uses an independent color guide", {
  gene_data <- list(
    transcripts = data.table::data.table(
      transcript_id = "Tx1",
      track_y = 1,
      tx_xstart = 0.5,
      tx_xend = 2.5
    ),
    segments = data.table::data.table(
      transcript_id = "Tx1",
      track_y = 1,
      feature = "CDS",
      plot_start = 0.5,
      plot_end = 2.5
    ),
    region = list(chrom = "chr1", start = 100, end = 200),
    axis_breaks = data.table::data.table(
      x = c(0.5, 2.5),
      pos = c(100, 200),
      label = c("100", "200")
    )
  )
  vars <- data.table::data.table(
    pos = c(125, 175),
    gene_x = c(1, 2),
    variant_type = c("snp", "indel"),
    variant_type_label = c("SNP", "Ind")
  )
  p <- draw_hap_gene_track(
    gene_data = gene_data,
    vars = vars,
    x_limits = c(0.5, 2.5),
    x_breaks = c(1, 2),
    x_labels = c("125", "175"),
    gene_track_legend_position = "top"
  )

  expect_identical(p$scales$get_scales("fill")$guide, "none")
  color_guide <- p$scales$get_scales("colour")$guide
  expect_false(is.null(color_guide))
  expect_false(identical(color_guide, "none"))
  expect_identical(
    p$scales$get_scales("colour")$limits,
    c("SNP", "Ind")
  )
  point_layers <- vapply(
    p$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1L)
  )
  expect_true(any(point_layers))
  expect_true(any(vapply(
    p$layers[point_layers],
    function(layer) isTRUE(layer$show.legend),
    logical(1L)
  )))
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

  trait_order <- c("seed_weight", "plant_height")
  multi <- plot_hap_pheno(
    hap,
    phenotype = pheno,
    traits = trait_order,
    min_hap_samples = 1,
    facet_ncol = 1
  )
  expect_null(multi$figure$labels$title)
  expect_identical(levels(multi$figure$data$trait), trait_order)
  expect_identical(multi$figure$facet$params$ncol, 1L)
})

test_that("custom missing genotype labels are excluded from variant phenotype groups", {
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)

  direct <- extract_single_variant_genotype(
    variant = vcf,
    variant_id = "rs001_2",
    missing_genotype = "-"
  )
  expect_true(is.na(direct$genotype$genotype_group[direct$genotype$sample_id == "S12"]))

  hap <- hap_region_variant(
    vcf,
    chrom = "chr1",
    start = 1,
    end = 1200,
    genotype_mode = "code",
    missing_genotype = "-",
    min_variant_number = 1
  )
  from_hap <- extract_single_variant_genotype(
    variant = hap,
    variant_id = "rs001_2"
  )
  expect_true(is.na(from_hap$genotype$genotype_group[from_hap$genotype$sample_id == "S12"]))
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

test_that("plot_variant_effect uses signed-effect colors", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("example_haplotype.vcf"), mode = "memory", verbose = FALSE)
  pheno <- read_pheno(gtr_extdata("example_pheno.tsv"), verbose = FALSE)
  hap <- hap_gene_variant(
    vcf,
    annotation = gp,
    gene_id = "GeneA",
    genotype_mode = "code",
    min_variant_number = 1
  )

  res <- plot_variant_effect(
    hap,
    phenotype = pheno,
    traits = "plant_height",
    min_group_samples = 1
  )

  expect_true(inherits(res$figure, "ggplot"))
  color_scale <- res$figure$scales$get_scales("colour")
  expect_false(is.null(color_scale))
  expect_identical(color_scale$name, "Signed effect")

  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    res$figure$layers
  )
  expect_gt(length(point_layers), 0L)
  expect_false(is.null(point_layers[[1L]]$mapping$colour))
  expect_null(point_layers[[1L]]$aes_params$colour)

  axis_text_x <- res$figure$theme$axis.text.x
  expect_false(is.null(axis_text_x))
  expect_equal(axis_text_x$angle, 90)
  expect_equal(axis_text_x$hjust, 1)
  expect_equal(axis_text_x$vjust, 0.5)
})

