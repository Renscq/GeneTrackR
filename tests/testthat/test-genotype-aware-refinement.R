refinement_demo <- function() {
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory",
                  verbose = FALSE, progress = FALSE)
  annotation <- read_genepred(gtr_extdata("gtr_demo.genePredExt"),
                              format = "genePredExt", verbose = FALSE, progress = FALSE)
  hap <- hap_variant(vcf, annotation = annotation, gene_id = "GeneA")
  pheno <- read_pheno(gtr_extdata("gtr_demo_pheno.tsv"), verbose = FALSE)
  list(hap = hap, pheno = pheno, annotation = annotation)
}

test_that("genotype-aware refinement protects the GeneA signal and retains mappings", {
  demo <- refinement_demo()
  ids <- demo$hap$sample_haplotypes[, c("sample_id", "hap_id")]
  pheno <- as.data.frame(ids)
  pheno$value <- ave(seq_len(nrow(pheno)), pheno$hap_id, FUN = seq_along) / 100 +
    ifelse(pheno$hap_id %in% c("Hap3", "Hap4"), 10, 0)
  result <- refine_haplotype(
    demo$hap, pheno, traits = "value", strategy = "genotype_aware",
    control = list(direction = "higher", equivalence_margin = 0.5,
                   protect_variants = "varA03", max_diff_blocks = 2L)
  )
  expect_s3_class(result, "HapRefined")
  expect_equal(sort(as.integer(table(result$variant_map$block_id))), c(1L, 1L, 9L))
  main_block <- result$variant_map$block_id[result$variant_map$variant_id == "varA03"]
  expect_true(result$variant_rank$protected[match(main_block, result$variant_rank$block_id)])
  expect_equal(length(unique(result$sample_map$sample_id)), 36L)
  expect_equal(nrow(result$refined_haplotypes), 2L)
  expect_gt(result$haplotype_rank$rank_score[1L], result$haplotype_rank$rank_score[2L])
  expect_true(any(result$merge_audit$reason == "protected_genotype"))
  expect_equal(sum(result$merge_audit$can_merge), 2L)
  expect_gt(nrow(result$variant_contrasts), 0L)
  expect_s3_class(plot_refined_hap_variant(result, variant_view = "blocks"), "ggplot")
  expect_s3_class(plot_refined_hap_variant(result, annotation = NULL), "ggplot")
})

test_that("holdout phenotype values cannot change discovery rankings or group membership", {
  demo <- refinement_demo()
  sample_haps <- as.data.frame(demo$hap$sample_haplotypes)
  held <- sample_haps$sample_id[!duplicated(sample_haps$hap_id)]
  pheno <- as.data.frame(demo$pheno)
  control <- list(direction = "higher", equivalence_margin = 0.5,
                  validation_samples = held, max_diff_blocks = 2L)
  fit <- function(x) refine_haplotype(demo$hap, x, traits = "protein_content",
                                     strategy = "genotype_aware", control = control)
  first <- fit(pheno)
  pheno$protein_content[pheno$sample_id %in% held] <- 999
  second <- fit(pheno)
  expect_identical(first$haplotype_map, second$haplotype_map)
  expect_identical(first$haplotype_rank, second$haplotype_rank)
  expect_identical(first$variant_rank, second$variant_rank)
  expect_setequal(first$validation$predictions$sample_id, held)
  expect_gt(second$validation$refined_rmse, first$validation$refined_rmse)
})

test_that("refined variant blocks share original coordinates and allele labels", {
  demo <- refinement_demo()
  result <- refine_haplotype(
    demo$hap, demo$pheno, traits = "protein_content", strategy = "genotype_aware",
    control = list(direction = "higher", equivalence_margin = 0.5,
                   protect_variants = "varA03", max_diff_blocks = 2L)
  )
  original <- plot_hap_variant(demo$hap, annotation = demo$annotation)
  refined <- plot_refined_hap_variant(
    result, annotation = demo$annotation, variant_view = "blocks", collapse_refined = TRUE
  )
  expect_s3_class(refined, "patchwork")
  raw_vars <- attr(original, "variant_data")
  refined_vars <- attr(refined, "variant_data")
  expect_identical(refined_vars$variant_id, raw_vars$variant_id)
  expect_identical(refined_vars$pos, raw_vars$pos)
  expect_equal(refined_vars$gene_x, raw_vars$gene_x)
  singleton <- refined_vars$variant_id %in% c("varA02", "varA05")
  connector <- ggplot2::ggplot_build(refined[[2L]])$data
  branches <- connector[[1L]]
  trunk <- connector[[2L]]
  crossing_gap <- connector[[3L]]
  direct <- connector[[4L]]
  junction <- connector[[5L]]
  expect_equal(branches$x, refined_vars$gene_x[!singleton])
  expect_equal(branches$xend, refined_vars$table_x[!singleton])
  expect_equal(trunk$x, unique(branches$xend))
  expect_true(all(branches$yend == trunk$y) && all(trunk$yend == 0))
  expect_equal(direct$x, refined_vars$gene_x[singleton])
  expect_equal(direct$xend, refined_vars$table_x[singleton])
  expect_true(all(direct$yend == 0) && all(direct$colour != branches$colour[1L]))
  expect_equal(crossing_gap$xend, direct$xend)
  expect_true(all(crossing_gap$linewidth > direct$linewidth))
  expect_equal(junction$x, trunk$x)
  expect_equal(junction$y, trunk$y)
  expect_equal(refined_vars$table_x[refined_vars$variant_id %in% c("varA02", "varA05")],
               raw_vars$table_x[raw_vars$variant_id %in% c("varA02", "varA05")])
  expect_length(unique(refined_vars$table_x[!refined_vars$variant_id %in% c("varA02", "varA05")]), 1L)
  cells <- attr(refined, "plot_data")
  expect_setequal(unique(cells$variant_label), c("rV001_9", "varA02", "varA05"))
  expect_true(all(c("refined_hap_id", "row", "block_id", "state") %in% names(cells)))
  expect_true(any(cells$state == "1/1" & cells$genotype_label == "G-G-A-T-T-C-G-C-i3", na.rm = TRUE))
  expect_true(any(grepl("Hap3: T\nHap4: C", cells$genotype_label, fixed = TRUE)))
  expect_false(any(grepl("mixed", cells$genotype_label, fixed = TRUE)))
  expect_true(any(cells$genotype_label == "A-C-G-A-C-G-T-A-A" &
                    cells$genotype_fill_class == "A", na.rm = TRUE))
  expect_equal(cells$allele_label[cells$variant_label == "rV001_9" & cells$row_label == "REF"],
               "A-C-G-A-C-G-T-A-A")
  expect_equal(cells$allele_label[cells$variant_label == "rV001_9" & cells$row_label == "ALT"],
               "G-G-A-T-T-C-G-C-i3")
  table_layers <- ggplot2::ggplot_build(refined[[3L]])$data
  tiles <- unique(table_layers[[1L]][, c("x", "xmin", "xmax")])
  tiles <- tiles[order(tiles$xmin), ]
  expect_true(all(head(tiles$xmax, -1L) < tail(tiles$xmin, -1L)))
  expect_equal(tiles$xmax[tiles$x == trunk$x] - tiles$xmin[tiles$x == trunk$x], 2.96)
  expect_equal(tiles$x[tiles$x != trunk$x], as.numeric(raw_vars$table_x[singleton]))
  expect_true(any(table_layers[[3L]]$label == "G-G-A-T-T-C-G-C-i3"))
  expanded <- plot_refined_hap_variant(
    result, variant_view = "blocks", variant_label = "pos", show_reference_row = FALSE
  )
  expanded_cells <- attr(expanded, "plot_data")
  expect_setequal(unique(expanded_cells$variant_label),
                  c("rV001_9", as.character(raw_vars$pos[raw_vars$variant_id %in% c("varA02", "varA05")])))
  expect_true(all(expanded_cells$hap_id %in% result$sample_map$pattern_id))
})

test_that("block sequences decode each member in its own allele orientation", {
  demo <- refinement_demo()
  controls <- list(direction = "higher", protect_variants = "varA03")
  baseline <- refine_haplotype(demo$hap, demo$pheno, traits = "protein_content",
                               strategy = "genotype_aware", control = controls)
  controls$genotypes <- baseline$genotype_dosage
  controls$genotypes[, "varA03"] <- 2 - controls$genotypes[, "varA03"]
  flipped <- refine_haplotype(demo$hap, demo$pheno, traits = "protein_content",
                              strategy = "genotype_aware", control = controls)
  expect_true(flipped$variant_map$flipped[flipped$variant_map$variant_id == "varA03"])
  plot <- plot_refined_hap_variant(flipped, variant_view = "blocks", collapse_refined = TRUE)
  cells <- attr(plot, "plot_data")
  expect_true(any(cells$state == "1/1" & cells$genotype_label == "G-C-A-T-T-C-G-C-i3", na.rm = TRUE))
})

test_that("chained pairwise compatibility cannot merge an incompatible third pattern", {
  edges <- data.table::data.table(
    i = c(1L, 2L), j = c(2L, 3L), n_differing = c(1L, 1L),
    can_merge = c(TRUE, TRUE), p_value = c(0.001, 0.002)
  )
  expect_identical(refine_merge_clusters(3L, edges), c(1L, 1L, 3L))
})

test_that("explicit combination candidates are retained without locking their genotypes", {
  demo <- refinement_demo()
  result <- refine_haplotype(
    demo$hap, demo$pheno, traits = "flowering_time", strategy = "genotype_aware",
    control = list(direction = "lower", candidate_variants = "varA02")
  )
  map <- result$variant_map
  block <- map$block_id[map$variant_id == "varA02"]
  expect_true(result$variant_rank$in_panel[result$variant_rank$block_id == block])
  expect_false(any(map$protected[map$block_id == block]))
  expect_equal(nrow(result$refined_haplotypes), 4L)
})

test_that("genotype-aware refinement separates missing phenotypes from carrier counts", {
  demo <- refinement_demo()
  keep <- demo$hap$sample_haplotypes$sample_id[demo$hap$sample_haplotypes$hap_id == "Hap1"]
  pheno <- as.data.frame(demo$pheno)
  pheno$protein_content[!pheno$sample_id %in% keep] <- NA_real_
  result <- refine_haplotype(
    demo$hap, pheno, traits = "protein_content", strategy = "genotype_aware",
    control = list(direction = "lower", equivalence_margin = 0.5)
  )
  expect_equal(sum(result$haplotype_rank$n_phenotype), length(keep))
  expect_equal(sum(result$haplotype_rank$n_carrier), 36L)
  expect_equal(nrow(result$haplotype_map), 4L)
  expect_equal(sum(!is.na(result$haplotype_rank$rank)), 1L)
  expect_error(refine_haplotype(
    demo$hap, rbind(pheno, pheno[1L, ]), traits = "protein_content",
    strategy = "genotype_aware", control = list(direction = "higher")
  ), "unique")
})

test_that("legacy refinement returns a single group without pairwise comparisons", {
  demo <- refinement_demo()
  keep <- demo$hap$sample_haplotypes$sample_id[demo$hap$sample_haplotypes$hap_id == "Hap1"]
  result <- refine_haplotype(
    demo$hap, demo$pheno[demo$pheno$sample_id %in% keep, ], traits = "protein_content"
  )
  expect_equal(nrow(result$refined_haplotypes), 1L)
  expect_equal(nrow(result$pairwise_test), 0L)
})
