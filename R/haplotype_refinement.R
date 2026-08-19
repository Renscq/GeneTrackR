# Author: Rensc
# Date: 2026-07-30
# Version: dev003
# Function: Refine haplotypes and prioritize phenotype-associated variant effects
# Input: HapVariant objects and phenotype tables
# Output: Refined haplotype objects, variant-effect tables, and ggplot figures

#' Refine haplotypes by phenotype similarity
#'
#' @description
#' Collapse original haplotypes into refined haplotype groups according to
#' phenotype similarity. Haplotypes are merged only when their pairwise
#' phenotype difference is not significant and, optionally, the absolute mean
#' difference is smaller than `effect_threshold`.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_hap_samples Minimum sample number required for an original haplotype.
#' @param test_method Pairwise test method. One of `t.test`, `wilcox.test`, or `ks.test`.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param alpha Adjusted p-value cutoff. Pairs with adjusted p-value larger than
#' `alpha` are considered statistically indistinguishable.
#' @param effect_threshold Optional maximum absolute mean difference for merging
#' haplotypes. If NULL, only the significance criterion is used.
#' @param group_prefix Prefix for refined haplotype IDs.
#' @param mixed_label Label used for variant states that are heterogeneous within
#' a refined haplotype group.
#' @return A HapRefined object. The object contains `refined_hap`, which is a
#' HapVariant-compatible object and can be passed to `plot_hap_pheno()` and
#' `plot_hap_variant()`.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' pheno <- read_pheno(pheno_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string")
#' refined <- refine_haplotype(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' refined$refined_haplotypes
#' plot_hap_pheno(refined$refined_hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' @export
refine_haplotype <- function(hap,
                             phenotype,
                             traits = NULL,
                             sample_col = "sample_id",
                             min_hap_samples = 2L,
                             test_method = c("t.test", "wilcox.test", "ks.test"),
                             p_adjust = "BH",
                             alpha = 0.05,
                             effect_threshold = NULL,
                             group_prefix = "RHap",
                             mixed_label = "mixed") {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  test_method <- match.arg(test_method)

  min_hap_samples <- as.integer(min_hap_samples)[1L]
  stop_if_not(!is.na(min_hap_samples) && min_hap_samples >= 1L, "`min_hap_samples` must be a positive integer.")
  alpha <- as.numeric(alpha)[1L]
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  if (!is.null(effect_threshold)) {
    effect_threshold <- as.numeric(effect_threshold)[1L]
    stop_if_not(!is.na(effect_threshold) && effect_threshold >= 0, "`effect_threshold` must be NULL or a non-negative number.")
  }
  group_prefix <- as.character(group_prefix)[1L]
  if (is.na(group_prefix) || group_prefix == "") group_prefix <- "RHap"
  mixed_label <- as.character(mixed_label)[1L]
  if (is.na(mixed_label) || mixed_label == "") mixed_label <- "mixed"

  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  if (!identical(sample_col, "sample_id")) {
    data.table::setnames(pheno, sample_col, "sample_id")
  }
  pheno[, "sample_id" := as.character(sample_id)]

  hap_sample <- data.table::as.data.table(hap$sample_haplotypes)
  stop_if_not(all(c("sample_id", "hap_id") %in% names(hap_sample)), "`hap$sample_haplotypes` must contain sample_id and hap_id columns.")
  hap_sample[, "sample_id" := as.character(sample_id)]
  hap_sample[, "hap_id" := as.character(hap_id)]

  dt <- merge(hap_sample[, .(sample_id, hap_id)], pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between haplotypes and phenotype table.")

  hap_count <- dt[, .(sample_n = .N), by = hap_id]
  hap_keep <- hap_count[sample_n >= min_hap_samples, hap_id]
  dt <- dt[hap_id %in% hap_keep]
  hap_count <- hap_count[hap_id %in% hap_keep]
  stop_if_not(length(unique(dt$hap_id)) >= 1L, "No haplotypes with enough samples remained.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & !trait %in% c("hap_id", "sample_n"), trait]
  }
  traits <- as.character(traits)
  stop_if_not(length(traits) > 0L, "No numeric phenotype traits were selected.")
  stop_if_not(all(traits %in% names(dt)), "Some traits were not found in phenotype table.")

  long <- data.table::melt(
    dt,
    id.vars = c("sample_id", "hap_id"),
    measure.vars = traits,
    variable.name = "trait",
    value.name = "value",
    variable.factor = FALSE
  )
  long[, "value" := suppressWarnings(as.numeric(value))]
  long <- long[!is.na(value)]
  stop_if_not(nrow(long) > 0L, "No non-missing phenotype values were available.")

  summary_dt <- long[, .(
    sample_n = .N,
    mean_value = mean(value, na.rm = TRUE),
    median_value = stats::median(value, na.rm = TRUE),
    sd_value = stats::sd(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ), by = .(trait, hap_id)]

  test_dt <- pairwise_hap_test(long, method = test_method, p_adjust = p_adjust)
  if (nrow(test_dt) > 0L) {
    mean_tab <- summary_dt[, .(trait, hap_id, mean_value)]
    test_dt <- merge(test_dt, mean_tab, by.x = c("trait", "group1"), by.y = c("trait", "hap_id"), all.x = TRUE, sort = FALSE)
    data.table::setnames(test_dt, "mean_value", "mean_group1")
    test_dt <- merge(test_dt, mean_tab, by.x = c("trait", "group2"), by.y = c("trait", "hap_id"), all.x = TRUE, sort = FALSE)
    data.table::setnames(test_dt, "mean_value", "mean_group2")
    test_dt[, "mean_diff" := mean_group2 - mean_group1]
    test_dt[, "abs_mean_diff" := abs(mean_diff)]
    test_dt[, "can_merge" := !is.na(p_adj) & p_adj > alpha]
    if (!is.null(effect_threshold)) {
      test_dt[, "can_merge" := can_merge & !is.na(abs_mean_diff) & abs_mean_diff <= effect_threshold]
    }
  } else {
    test_dt[, c("mean_group1", "mean_group2", "mean_diff", "abs_mean_diff", "can_merge") := list(numeric(), numeric(), numeric(), numeric(), logical())]
  }

  trait_maps <- lapply(traits, function(trait_i) {
    make_refined_groups_for_trait(
      summary_dt = summary_dt[trait == trait_i],
      test_dt = test_dt[trait == trait_i],
      alpha = alpha,
      effect_threshold = effect_threshold
    )
  })
  trait_map <- data.table::rbindlist(trait_maps, use.names = TRUE, fill = TRUE)

  wide_map <- data.table::dcast(
    trait_map,
    hap_id ~ trait,
    value.var = "trait_group",
    fill = NA_character_
  )
  trait_cols <- setdiff(names(wide_map), "hap_id")
  wide_map[, "refined_key" := do.call(paste, c(lapply(.SD, as.character), sep = "|")), .SDcols = trait_cols]
  key_dt <- unique(wide_map[, .(refined_key)])
  key_dt[, "refined_hap_id" := paste0(group_prefix, seq_len(.N))]
  wide_map <- merge(wide_map, key_dt, by = "refined_key", all.x = TRUE, sort = FALSE)
  refined_map <- wide_map[, c("hap_id", "refined_hap_id", "refined_key", trait_cols), with = FALSE]
  data.table::setorder(refined_map, refined_hap_id, hap_id)

  sample_refined <- merge(hap_sample, refined_map[, .(hap_id, refined_hap_id)], by = "hap_id", all.x = FALSE, sort = FALSE)
  data.table::setnames(sample_refined, "hap_id", "original_hap_id")
  data.table::setnames(sample_refined, "refined_hap_id", "hap_id")
  data.table::setcolorder(sample_refined, c("sample_id", "hap_id", "original_hap_id", setdiff(names(sample_refined), c("sample_id", "hap_id", "original_hap_id"))))

  haps <- data.table::as.data.table(hap$haplotypes)
  haps[, "hap_id" := as.character(hap_id)]
  variant_ids <- get_hap_variant_columns(hap)
  haps <- merge(haps, refined_map[, .(hap_id, refined_hap_id)], by = "hap_id", all.x = FALSE, sort = FALSE)

  refined_haps <- haps[, .(
    original_hap_n = data.table::uniqueN(hap_id),
    original_hap_ids = paste(unique(hap_id), collapse = ";"),
    sample_n = sum(as.integer(sample_n), na.rm = TRUE),
    samples = paste(unique(unlist(strsplit(paste(samples, collapse = ";"), ";", fixed = TRUE))), collapse = ";")
  ), by = refined_hap_id]

  if (length(variant_ids) > 0L) {
    variant_state <- haps[, lapply(.SD, collapse_refined_variant_state, mixed_label = mixed_label), by = refined_hap_id, .SDcols = variant_ids]
    refined_haps <- merge(refined_haps, variant_state, by = "refined_hap_id", all.x = TRUE, sort = FALSE)
  }
  data.table::setnames(refined_haps, "refined_hap_id", "hap_id")
  data.table::setorder(refined_haps, -sample_n, hap_id)

  refined_sample_haplotypes <- sample_refined[, .(sample_id, hap_id, original_hap_id)]
  refined_sample_haplotypes[, "hap_pattern" := hap_id]
  refined_sample_haplotypes[, "non_missing_variant_n" := hap$meta$min_variant_number %||% length(variant_ids)]

  refined_hap <- hap
  refined_hap$haplotypes <- refined_haps[]
  refined_hap$sample_haplotypes <- refined_sample_haplotypes[]
  refined_hap$meta$original_haplotype_n <- hap$meta$haplotype_n %||% nrow(hap$haplotypes)
  refined_hap$meta$haplotype_n <- nrow(refined_haps)
  refined_hap$meta$refined <- TRUE
  refined_hap$meta$refine_traits <- traits
  refined_hap$meta$refine_alpha <- alpha
  refined_hap$meta$refine_effect_threshold <- effect_threshold
  class(refined_hap) <- "HapVariant"

  out <- list(
    original_hap = hap,
    refined_hap = refined_hap,
    refined_haplotypes = refined_haps[],
    sample_refined_haplotypes = sample_refined[],
    haplotype_map = refined_map[],
    trait_group_map = trait_map[],
    phenotype_summary = summary_dt[],
    pairwise_test = test_dt[],
    plot_data = long[],
    parameters = list(
      traits = traits,
      min_hap_samples = min_hap_samples,
      test_method = test_method,
      p_adjust = p_adjust,
      alpha = alpha,
      effect_threshold = effect_threshold,
      mixed_label = mixed_label
    )
  )
  class(out) <- c("HapRefined", "list")
  out
}

#' Plot phenotype distributions for refined haplotypes
#'
#' @description
#' Draw phenotype distributions for refined haplotype groups generated by
#' `refine_haplotype()`. This function is a refined-haplotype wrapper around
#' `plot_hap_pheno()`, but all commonly used phenotype-plot parameters are
#' exposed explicitly so that the help page shows the complete usage.
#'
#' @param refined_hap A `HapRefined` object returned by `refine_haplotype()` or
#' a refined `HapVariant` object.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible
#' data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_hap_samples Minimum sample number required for a refined haplotype group.
#' @param plot_type Plot type. One of `violin`, `boxplot`, or `violin_boxplot`.
#' @param test_method Pairwise test method. One of `t.test`, `wilcox.test`, or `ks.test`.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param p_label P-value label style. One of `stars`, `number`, or `both`.
#' @param p_cutoff Significance cutoff used for displaying pairwise comparisons.
#' @param p_value_type Which p-value is used for filtering and labeling. One of `raw` or `adjusted`.
#' @param show_signif_only Logical. Whether to only display significant comparisons.
#' @param show_points Logical. Whether to show sample points.
#' @param show_outliers Logical. Whether to show boxplot outliers.
#' @param fill_palette RColorBrewer palette name used for median-based refined haplotype fills.
#' @param fill_colors Optional custom fill colors for refined haplotypes.
#' @param fill_alpha Alpha value for violin/boxplot fill colors.
#' @param violin_width Violin plot width.
#' @param box_width Boxplot width.
#' @param bracket_step Fraction of y-range used to separate significance brackets.
#' @param bracket_tip_fraction Fraction of bracket vertical spacing used for the short downward bracket tips.
#' @param x_text_angle Rotation angle for refined haplotype labels on the x-axis.
#' @param strip_label_width Maximum character width for wrapping long facet strip labels.
#' @param strip_text_lineheight Line height for wrapped facet strip labels.
#' @param strip_fill Facet strip background fill color.
#' @param strip_border_color Facet strip border color. Use NULL to remove the border.
#' @param text_size Text size.
#' @return A `GeneTrackRPhenoPlot` object returned by `plot_hap_pheno()`.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' pheno <- read_pheno(pheno_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string")
#' refined <- refine_haplotype(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' plot_refined_hap_pheno(refined, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' @export
plot_refined_hap_pheno <- function(refined_hap,
                                   phenotype,
                                   traits = NULL,
                                   sample_col = "sample_id",
                                   min_hap_samples = 2L,
                                   plot_type = c("violin_boxplot", "violin", "boxplot"),
                                   test_method = c("t.test", "wilcox.test", "ks.test"),
                                   p_adjust = "BH",
                                   p_label = c("stars", "number", "both"),
                                   p_cutoff = 0.05,
                                   p_value_type = c("raw", "adjusted"),
                                   show_signif_only = TRUE,
                                   show_points = FALSE,
                                   show_outliers = FALSE,
                                   fill_palette = "RdBu",
                                   fill_colors = NULL,
                                   fill_alpha = 0.75,
                                   violin_width = 0.9,
                                   box_width = 0.18,
                                   bracket_step = 0.08,
                                   bracket_tip_fraction = 0.12,
                                   x_text_angle = 90,
                                   strip_label_width = 24,
                                   strip_text_lineheight = 0.9,
                                   strip_fill = "white",
                                   strip_border_color = NULL,
                                   text_size = 14) {
  hap_obj <- get_refined_hap_variant(refined_hap)
  plot_hap_pheno(
    hap = hap_obj,
    phenotype = phenotype,
    traits = traits,
    sample_col = sample_col,
    min_hap_samples = min_hap_samples,
    plot_type = plot_type,
    test_method = test_method,
    p_adjust = p_adjust,
    p_label = p_label,
    p_cutoff = p_cutoff,
    p_value_type = p_value_type,
    show_signif_only = show_signif_only,
    show_points = show_points,
    show_outliers = show_outliers,
    fill_palette = fill_palette,
    fill_colors = fill_colors,
    fill_alpha = fill_alpha,
    violin_width = violin_width,
    box_width = box_width,
    bracket_step = bracket_step,
    bracket_tip_fraction = bracket_tip_fraction,
    x_text_angle = x_text_angle,
    strip_label_width = strip_label_width,
    strip_text_lineheight = strip_text_lineheight,
    strip_fill = strip_fill,
    strip_border_color = strip_border_color,
    text_size = text_size
  )
}

#' Plot grouped variant table for refined haplotypes
#'
#' @description
#' Draw a refined haplotype-variant table. By default, when the input is a
#' `HapRefined` object, the original haplotype rows are retained and arranged
#' into refined haplotype blocks. For example, if ten original haplotypes are
#' refined into three groups, the table shows the ten original haplotypes split
#' into three vertically separated blocks. This makes it easier to inspect which
#' variants are shared or variable within each refined group.
#'
#' If `collapse_refined = TRUE`, the function falls back to a collapsed refined
#' table in which each refined haplotype is shown as one row, using the same
#' visual grammar as `plot_hap_variant()`.
#'
#' @param refined_hap A `HapRefined` object returned by `refine_haplotype()` or
#' a refined `HapVariant` object. The grouped display requires a `HapRefined`
#' object generated by the current version of `refine_haplotype()`.
#' @param annotation Optional Feature/GenePred annotation object. When supplied,
#' a compact gene track is drawn above the refined haplotype table.
#' @param show_gene_model Logical. Whether to draw the gene track when `annotation` is supplied.
#' @param min_hap_samples Minimum sample number required for an original haplotype row.
#' @param show_reference_row Logical. Whether to add two reference rows showing REF and ALT alleles for each variant.
#' @param collapse_refined Logical. Whether to plot one collapsed row per refined haplotype group.
#' @param group_gap Numeric vertical gap inserted between refined haplotype blocks.
#' @param show_group_label Logical. Whether to show refined group labels on the left side of the table.
#' @param group_label_prefix Prefix used for left-side group labels.
#' @param group_label_hjust Horizontal adjustment for refined group labels. Larger
#' values move the group labels further left and reduce overlap with original
#' haplotype row labels.
#' @param group_label_size Text size for refined group labels. If NULL, it is
#' derived from `text_size`.
#' @param group_label_left_margin Left plot margin used when refined group labels
#' are shown. Increase this value when group labels are clipped.
#' @param variant_label Column used for variant labels. One of `variant_id`, `pos`, or an existing column in `hap$variants`.
#' @param show_gene_pos_axis Logical. Whether to show genomic coordinate labels above the gene track.
#' @param gene_pos_axis_n Approximate number of genomic coordinate ticks above the gene track.
#' @param gene_pos_axis_label Optional x-axis title for genomic coordinate labels.
#' @param gene_pos_x_angle Angle of gene-position x-axis labels.
#' @param gene_track_legend_position Legend position for variant-type markers in the gene track. One of `right`, `top`, or `none`.
#' @param text_size Base text size for gene-track labels, axes, legends, and table axes.
#' @param table_x_angle Angle of haplotype table x-axis labels.
#' @param genotype_text_size Genotype cell text size only.
#' @param gene_track_height Relative height of the gene track panel.
#' @param connector_height Relative height of the connector panel.
#' @param table_height Relative height of the haplotype table panel. If NULL, it is inferred from the number of displayed rows and group gaps.
#' @param exon_height Height of exon boxes in the compact gene track.
#' @param cds_height Height of CDS boxes in the compact gene track.
#' @param gene_palette RColorBrewer palette name used for gene model feature fills.
#' @param gene_colors Optional custom fill colors for gene model features.
#' @param gene_border_color Border color for gene model feature boxes.
#' @param table_palette RColorBrewer palette name used for haplotype table fills.
#' For allele-string haplotypes, the first five colors are always assigned in
#' the fixed order A, T, C, G, and indel.
#' @param table_colors Optional custom fill colors for haplotype table genotypes.
#' For allele-string haplotypes, use a vector named `A`, `T`, `C`, `G`, and
#' `indel`, or supply five unnamed colors in that order.
#' @param table_alpha Alpha value for haplotype table fill colors.
#' @param reference_fill Background fill color for the REF and ALT reference rows.
#' @param variant_palette RColorBrewer palette name used for solid variant-type triangle marker colors.
#' @param variant_colors Optional custom colors for solid variant-type triangle markers.
#' @return A patchwork object with attributes `plot_data`, `variant_data`,
#' `gene_data`, `refined_haplotype_map`, and `refined_haplotypes` when available.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' pheno <- read_pheno(pheno_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string")
#' refined <- refine_haplotype(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' plot_refined_hap_variant(refined, annotation = anno, min_hap_samples = 1)
#' @export
plot_refined_hap_variant <- function(refined_hap,
                                     annotation = NULL,
                                     show_gene_model = TRUE,
                                     min_hap_samples = 5L,
                                     show_reference_row = TRUE,
                                     collapse_refined = FALSE,
                                     group_gap = 0.8,
                                     show_group_label = TRUE,
                                     group_label_prefix = "",
                                     group_label_hjust = 3.4,
                                     group_label_size = NULL,
                                     group_label_left_margin = NULL,
                                     variant_label = c("variant_id", "pos"),
                                     show_gene_pos_axis = TRUE,
                                     gene_pos_axis_n = 5L,
                                     gene_pos_axis_label = NULL,
                                     gene_pos_x_angle = 0,
                                     gene_track_legend_position = c("right", "top", "none"),
                                     text_size = 14,
                                     table_x_angle = 90,
                                     genotype_text_size = 3.2,
                                     gene_track_height = 1.25,
                                     connector_height = 0.35,
                                     table_height = NULL,
                                     exon_height = 0.22,
                                     cds_height = 0.44,
                                     gene_palette = "Paired",
                                     gene_colors = NULL,
                                     gene_border_color = NA,
                                     table_palette = "RdBu",
                                     table_colors = NULL,
                                     table_alpha = 0.6,
                                     reference_fill = "white",
                                     variant_palette = "Set1",
                                     variant_colors = NULL) {
  if (!inherits(refined_hap, "HapRefined") || isTRUE(collapse_refined)) {
    hap_obj <- get_refined_hap_variant(refined_hap)
    p <- plot_hap_variant(
      hap = hap_obj,
      annotation = annotation,
      show_gene_model = show_gene_model,
      min_hap_samples = min_hap_samples,
      show_reference_row = show_reference_row,
      variant_label = variant_label,
      show_gene_pos_axis = show_gene_pos_axis,
      gene_pos_axis_n = gene_pos_axis_n,
      gene_pos_axis_label = gene_pos_axis_label,
      gene_pos_x_angle = gene_pos_x_angle,
      gene_track_legend_position = gene_track_legend_position,
      text_size = text_size,
      table_x_angle = table_x_angle,
      genotype_text_size = genotype_text_size,
      gene_track_height = gene_track_height,
      connector_height = connector_height,
      table_height = table_height,
      exon_height = exon_height,
      cds_height = cds_height,
      gene_palette = gene_palette,
      gene_colors = gene_colors,
      gene_border_color = gene_border_color,
      table_palette = table_palette,
      table_colors = table_colors,
      table_alpha = table_alpha,
      reference_fill = reference_fill,
      variant_palette = variant_palette,
      variant_colors = variant_colors
    )
    if (inherits(refined_hap, "HapRefined")) {
      attr(p, "refined_haplotype_map") <- refined_hap$haplotype_map
      attr(p, "refined_haplotypes") <- refined_hap$refined_haplotypes
    }
    return(p)
  }

  stop_if_not(inherits(refined_hap$original_hap, "HapVariant"),
              "`refined_hap` does not contain the original HapVariant object. Please rerun refine_haplotype() with the updated GeneTrackR version.")

  original_hap <- refined_hap$original_hap
  vars <- data.table::as.data.table(original_hap$variants)
  haps <- data.table::as.data.table(original_hap$haplotypes)
  map <- data.table::as.data.table(refined_hap$haplotype_map)

  stop_if_not(nrow(vars) > 0L, "No variants are available in the original haplotype object.")
  stop_if_not(nrow(haps) > 0L, "No haplotypes are available in the original haplotype object.")
  stop_if_not(all(c("hap_id", "refined_hap_id") %in% names(map)), "`refined_hap$haplotype_map` must contain hap_id and refined_hap_id columns.")

  min_hap_samples <- as.integer(min_hap_samples)[1L]
  stop_if_not(!is.na(min_hap_samples) && min_hap_samples >= 1L, "`min_hap_samples` must be a positive integer.")
  show_reference_row <- isTRUE(show_reference_row)
  show_group_label <- isTRUE(show_group_label)
  group_gap <- as.numeric(group_gap)[1L]
  if (is.na(group_gap) || group_gap < 0) group_gap <- 0.8
  group_label_prefix <- as.character(group_label_prefix)[1L]
  if (is.na(group_label_prefix)) group_label_prefix <- ""
  group_label_hjust <- as.numeric(group_label_hjust)[1L]
  if (is.na(group_label_hjust) || group_label_hjust <= 0) group_label_hjust <- 3.4
  if (is.null(group_label_size)) {
    group_label_size <- text_size / 3.5
  }
  group_label_size <- as.numeric(group_label_size)[1L]
  if (is.na(group_label_size) || group_label_size <= 0) group_label_size <- text_size / 3.5
  if (is.null(group_label_left_margin)) {
    group_label_left_margin <- 88 + max(0, nchar(group_label_prefix) - 4) * 3
  }
  group_label_left_margin <- as.numeric(group_label_left_margin)[1L]
  if (is.na(group_label_left_margin) || group_label_left_margin < 8) group_label_left_margin <- 88

  gene_border_color <- normalize_border_color(gene_border_color)
  table_alpha <- as.numeric(table_alpha)[1L]
  if (is.na(table_alpha)) table_alpha <- 0.6
  table_alpha <- max(0, min(1, table_alpha))
  genotype_text_size <- as.numeric(genotype_text_size)[1L]
  if (is.na(genotype_text_size) || genotype_text_size <= 0) genotype_text_size <- 3.2
  table_x_angle <- as.numeric(table_x_angle)[1L]
  if (is.na(table_x_angle)) table_x_angle <- 90
  table_x_angle <- max(0, min(180, table_x_angle))
  gene_pos_x_angle <- as.numeric(gene_pos_x_angle)[1L]
  if (is.na(gene_pos_x_angle)) gene_pos_x_angle <- 0
  gene_pos_x_angle <- max(0, min(180, gene_pos_x_angle))
  reference_fill <- as.character(reference_fill)[1L]
  if (is.na(reference_fill) || reference_fill == "") reference_fill <- "white"

  haps[, "hap_id" := as.character(hap_id)]
  map[, "hap_id" := as.character(hap_id)]
  map[, "refined_hap_id" := as.character(refined_hap_id)]
  haps <- merge(haps, map[, .(hap_id, refined_hap_id)], by = "hap_id", all.x = FALSE, sort = FALSE)
  haps <- haps[sample_n >= min_hap_samples]
  stop_if_not(nrow(haps) > 0L, "No original haplotypes remained after `min_hap_samples` filtering.")

  variant_ids <- vars[["variant_id"]]
  variant_ids <- variant_ids[variant_ids %in% names(haps)]
  stop_if_not(length(variant_ids) > 0L, "No variant genotype columns were found in the original haplotype table.")

  vars <- vars[match(variant_ids, variant_id)]
  vars[, "variant_index" := seq_len(.N)]

  variant_label <- match.arg(variant_label)
  gene_track_legend_position <- match.arg(gene_track_legend_position)
  show_gene_pos_axis <- isTRUE(show_gene_pos_axis)
  gene_pos_axis_n <- as.integer(gene_pos_axis_n)[1L]
  if (is.na(gene_pos_axis_n) || gene_pos_axis_n < 2L) gene_pos_axis_n <- 5L
  if (is.null(gene_pos_axis_label)) gene_pos_axis_label <- NA_character_
  gene_pos_axis_label <- as.character(gene_pos_axis_label)[1L]

  if (!variant_label %in% names(vars)) variant_label <- "variant_id"
  vars[, "variant_label" := as.character(.SD[[1L]]), .SDcols = variant_label]
  vars[is.na(variant_label) | variant_label == "", "variant_label" := as.character(variant_id)]

  refined_order <- NULL
  if (!is.null(refined_hap$refined_haplotypes) && "hap_id" %in% names(refined_hap$refined_haplotypes)) {
    refined_order <- as.character(refined_hap$refined_haplotypes$hap_id)
  }
  refined_order <- unique(c(refined_order, haps$refined_hap_id))
  refined_order <- refined_order[refined_order %in% haps$refined_hap_id]
  haps[, "refined_order" := match(refined_hap_id, refined_order)]
  haps[, "hap_sort_number" := parse_hap_numeric_id(hap_id)]
  data.table::setorderv(haps, c("refined_order", "sample_n", "hap_sort_number", "hap_id"), order = c(1L, -1L, 1L, 1L))

  display_rows <- haps[, .(hap_id, refined_hap_id, sample_n)]
  display_rows[, "display_index" := seq_len(.N)]
  display_rows[, "group_index" := match(refined_hap_id, refined_order)]
  group_n <- length(refined_order)
  row_n <- nrow(display_rows)
  display_rows[, "hap_y" := row_n - display_index + 1 + (group_n - group_index) * group_gap]
  display_rows[, "row_label" := paste0(hap_id, " (n=", sample_n, ")")]

  id_vars <- intersect(c("hap_id", "sample_n", "samples", "refined_hap_id"), names(haps))
  table_long <- data.table::melt(
    haps,
    id.vars = id_vars,
    measure.vars = variant_ids,
    variable.name = "variant_id",
    value.name = "genotype",
    variable.factor = FALSE
  )
  table_long <- merge(
    table_long,
    vars[, .(variant_id, variant_index, variant_label, ref, alt)],
    by = "variant_id",
    all.x = TRUE,
    sort = FALSE
  )
  table_long <- merge(
    table_long,
    display_rows[, .(hap_id, refined_hap_id, hap_y, row_label, group_index)],
    by = c("hap_id", "refined_hap_id"),
    all.x = TRUE,
    sort = FALSE
  )
  table_long[, "genotype_label" := format_hap_table_genotype_label(
    genotype = genotype,
    ref = ref,
    alt = alt,
    genotype_mode = original_hap$meta$genotype_mode %||% NA_character_
  )]
  table_long[is.na(genotype_label), "genotype_label" := "NA"]
  table_long[, "genotype_fill_class" := normalize_hap_table_fill_class(
    genotype_label,
    genotype_mode = original_hap$meta$genotype_mode %||% NA_character_
  )]

  y_breaks <- display_rows$hap_y
  y_labels <- display_rows$row_label
  y_min <- min(display_rows$hap_y) - 0.5
  y_max <- max(display_rows$hap_y) + 0.5
  table_genotype <- table_long
  table_reference <- table_long[0]

  if (show_reference_row) {
    alt_y <- y_max + 0.65
    ref_y <- y_max + 1.55
    ref_row_label <- "REF"
    alt_row_label <- "ALT"
    alt_label_values <- vapply(strsplit(as.character(vars$alt), ",", fixed = TRUE), function(x) {
      x <- x[!is.na(x) & nzchar(x)]
      x <- format_hap_allele(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0L) return("NA")
      paste(x, collapse = ",")
    }, character(1L))
    ref_row <- vars[, .(
      variant_id,
      variant_index,
      variant_label,
      row_label = ref_row_label,
      refined_hap_id = ref_row_label,
      hap_y = ref_y,
      genotype = format_hap_allele(ref),
      genotype_label = format_hap_allele(ref),
      allele_label = format_hap_allele(ref)
    )]
    alt_row <- vars[, .(
      variant_id,
      variant_index,
      variant_label,
      row_label = alt_row_label,
      refined_hap_id = alt_row_label,
      hap_y = alt_y,
      genotype = alt_label_values,
      genotype_label = alt_label_values,
      allele_label = alt_label_values
    )]
    ref_alt <- data.table::rbindlist(list(alt_row, ref_row), fill = TRUE)
    for (nm in setdiff(names(table_long), names(ref_alt))) ref_alt[, (nm) := NA]
    data.table::setcolorder(ref_alt, names(table_long))
    table_reference <- ref_alt
    y_breaks <- c(y_breaks, alt_y, ref_y)
    y_labels <- c(y_labels, alt_row_label, ref_row_label)
    y_max <- ref_y + 0.5
  }
  if (!"allele_label" %in% names(table_reference)) {
    table_reference[, "allele_label" := character()]
  }

  group_label_dt <- display_rows[, .(
    y = mean(range(hap_y)),
    first_y = max(hap_y),
    last_y = min(hap_y),
    row_n = .N
  ), by = refined_hap_id]
  group_label_dt[, "group_label" := paste0(group_label_prefix, refined_hap_id)]
  data.table::setorder(group_label_dt, -y)
  sep_dt <- group_label_dt[order(-y)]
  if (nrow(sep_dt) > 1L) {
    sep_dt <- sep_dt[-nrow(sep_dt)]
    sep_dt[, "sep_y" := last_y - group_gap / 2]
  } else {
    sep_dt <- sep_dt[0]
    sep_dt[, "sep_y" := numeric()]
  }

  x_limits <- c(0.5, length(variant_ids) + 0.5)
  x_breaks <- vars$variant_index
  x_labels <- vars$variant_label
  table_mapper <- make_hap_table_mapper(vars)
  gene_mapper <- make_hap_gene_mapper(vars, original_hap$region)
  vars[, "table_x" := table_mapper(pos)]
  vars[, "gene_x" := gene_mapper(pos)]

  gene_data <- NULL
  p_gene <- NULL
  if (isTRUE(show_gene_model) && !is.null(annotation)) {
    gene_data <- prepare_hap_gene_track(annotation, original_hap, vars, gene_mapper)
    if (!is.null(gene_data) && nrow(gene_data$transcripts) > 0L) {
      p_gene <- draw_hap_gene_track(
        gene_data = gene_data,
        vars = vars,
        x_limits = x_limits,
        x_breaks = x_breaks,
        x_labels = x_labels,
        text_size = text_size,
        gene_text_size = text_size,
        exon_height = exon_height,
        cds_height = cds_height,
        gene_palette = gene_palette,
        gene_colors = gene_colors,
        gene_border_color = gene_border_color,
        variant_palette = variant_palette,
        variant_colors = variant_colors,
        gene_track_legend_position = gene_track_legend_position,
        show_gene_pos_axis = show_gene_pos_axis,
        gene_pos_axis_n = gene_pos_axis_n,
        gene_pos_axis_label = gene_pos_axis_label,
        gene_pos_x_angle = gene_pos_x_angle
      )
    }
  }

  p_connector <- ggplot2::ggplot(vars) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$gene_x, xend = .data$table_x, y = 1, yend = 0),
      linewidth = 0.35,
      color = "grey35",
      lineend = "butt"
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 12, 0, 8))

  p_table <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = table_genotype,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, fill = .data$genotype_fill_class),
      color = "grey80",
      linewidth = 0.28,
      width = 0.96,
      height = 0.86
    ) +
    ggplot2::geom_tile(
      data = table_reference,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y),
      fill = reference_fill,
      color = "grey70",
      linewidth = 0.28,
      width = 0.96,
      height = 0.86
    ) +
    ggplot2::geom_text(
      data = table_genotype,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, label = .data$genotype_label),
      size = genotype_text_size,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = table_reference,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, label = .data$allele_label),
      size = genotype_text_size * 0.95,
      color = "black",
      vjust = 0.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = x_limits,
      expand = c(0, 0),
      position = "bottom"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(y_min, y_max),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = NULL, y = if (show_group_label) NULL else "Original haplotype") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = text_size, angle = table_x_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 12, 8, if (show_group_label) group_label_left_margin else 8)
    )

  if (nrow(sep_dt) > 0L) {
    p_table <- p_table + ggplot2::geom_hline(
      data = sep_dt,
      ggplot2::aes(yintercept = .data$sep_y),
      linewidth = 0.45,
      color = "grey45",
      linetype = "dashed"
    )
  }
  if (isTRUE(show_group_label) && nrow(group_label_dt) > 0L) {
    p_table <- p_table + ggplot2::geom_text(
      data = group_label_dt,
      ggplot2::aes(x = 0.5, y = .data$y, label = .data$group_label),
      inherit.aes = FALSE,
      hjust = group_label_hjust,
      size = group_label_size,
      fontface = "bold",
      color = "black"
    )
  }

  p_table <- apply_hap_table_fill_scale(
    p_table,
    values = table_genotype$genotype_label,
    table_palette = table_palette,
    table_colors = table_colors,
    table_alpha = table_alpha,
    fixed_allele_classes = identical(
      as.character(original_hap$meta$genotype_mode %||% NA_character_)[1L],
      "string"
    )
  )

  if (is.null(table_height)) {
    table_height <- max(1.2, row_n * 0.38 + max(0, group_n - 1) * group_gap * 0.28 + if (show_reference_row) 0.85 else 0)
  }

  if (!is.null(p_gene)) {
    p <- patchwork::wrap_plots(
      p_gene,
      p_connector,
      p_table,
      ncol = 1,
      heights = c(gene_track_height, connector_height, table_height)
    )
  } else {
    p_variant <- ggplot2::ggplot(vars, ggplot2::aes(x = .data$gene_x, y = 1)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$variant_label),
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        size = text_size / 3.5,
        color = "black"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.7, 1.45), expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        text = ggplot2::element_text(color = "black"),
        plot.margin = ggplot2::margin(8, 12, 0, 8)
      )
    p <- patchwork::wrap_plots(
      p_variant,
      p_connector,
      p_table,
      ncol = 1,
      heights = c(0.75, connector_height, table_height)
    )
  }

  attr(p, "plot_data") <- data.table::rbindlist(list(table_genotype, table_reference), fill = TRUE)[]
  attr(p, "variant_data") <- vars[]
  attr(p, "gene_data") <- gene_data
  attr(p, "refined_haplotype_map") <- refined_hap$haplotype_map
  attr(p, "refined_haplotypes") <- refined_hap$refined_haplotypes
  attr(p, "group_data") <- group_label_dt[]
  p
}

#' Plot phenotype distributions for refined haplotypes
#'
#' @description
#' Deprecated compatibility alias of `plot_refined_hap_pheno()`.
#'
#' @param refined_hap A `HapRefined` object returned by `refine_haplotype()` or
#' a refined `HapVariant` object.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible
#' data.frame.
#' @param ... Additional parameters passed to `plot_refined_hap_pheno()`.
#' @return A `GeneTrackRPhenoPlot` object returned by `plot_hap_pheno()`.
#' @export
plot_hap_refined <- function(refined_hap, phenotype, ...) {
  .Deprecated("plot_refined_hap_pheno")
  plot_refined_hap_pheno(refined_hap = refined_hap, phenotype = phenotype, ...)
}


get_refined_hap_variant <- function(refined_hap) {
  if (inherits(refined_hap, "HapRefined")) {
    stop_if_not(
      inherits(refined_hap$refined_hap, "HapVariant"),
      "`refined_hap$refined_hap` must be a HapVariant object."
    )
    return(refined_hap$refined_hap)
  }
  if (inherits(refined_hap, "HapVariant")) {
    is_refined <- isTRUE(refined_hap$meta$refined)
    if (!is_refined) {
      warning("`refined_hap` is a HapVariant object without `meta$refined = TRUE`. It will be plotted as supplied.", call. = FALSE)
    }
    return(refined_hap)
  }
  stop("`refined_hap` must be a HapRefined object returned by refine_haplotype() or a refined HapVariant object.", call. = FALSE)
}

#' Plot phenotype effect size of each natural variant
#'
#' @description
#' Calculate per-variant phenotype effect sizes from a HapVariant object and draw
#' an ordered effect-size plot. For a two-genotype variant, the effect is the
#' mean difference between the second and first genotype groups. For variants
#' with more than two genotype groups, the effect is the maximum group mean minus
#' the minimum group mean. Point color encodes the signed effect direction and
#' magnitude, with blue for negative effects and red for positive effects.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_group_samples Minimum sample number required for a genotype group.
#' @param test_method Significance test method. One of `t.test`, `wilcox.test`,
#' `anova`, or `kruskal.test`.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param effect_type Which effect value is plotted. One of `absolute` or `signed`.
#' @param top_n Optional number of top variants to label by absolute effect size.
#' @param variant_label Column used for variant labels. One of `variant_id`, `pos`,
#' or an existing column in `hap$variants`.
#' @param x_axis X-axis type. One of `index` or `position`.
#' @param point_size Point size.
#' @param label_size Text size for top-variant labels.
#' @param text_size Base text size.
#' @return A list with `figure`, `effect`, and `plot_data` elements.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' pheno <- read_pheno(pheno_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "code")
#' plot_variant_effect(hap, phenotype = pheno, traits = "plant_height", min_group_samples = 1)
#' @export
plot_variant_effect <- function(hap,
                                phenotype,
                                traits = NULL,
                                sample_col = "sample_id",
                                min_group_samples = 2L,
                                test_method = c("t.test", "wilcox.test", "anova", "kruskal.test"),
                                p_adjust = "BH",
                                effect_type = c("absolute", "signed"),
                                top_n = 10L,
                                variant_label = c("variant_id", "pos"),
                                x_axis = c("index", "position"),
                                point_size = 2.2,
                                label_size = 3.2,
                                text_size = 14) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  test_method <- match.arg(test_method)
  effect_type <- match.arg(effect_type)
  x_axis <- match.arg(x_axis)

  min_group_samples <- as.integer(min_group_samples)[1L]
  stop_if_not(!is.na(min_group_samples) && min_group_samples >= 1L, "`min_group_samples` must be a positive integer.")
  top_n <- as.integer(top_n)[1L]
  if (is.na(top_n) || top_n < 0L) top_n <- 0L

  variant_ids <- get_hap_variant_columns(hap)
  stop_if_not(length(variant_ids) > 0L, "No variant genotype columns were found in `hap$genotype_wide` or `hap$haplotypes`.")

  geno <- data.table::as.data.table(hap$genotype_wide)
  stop_if_not("sample_id" %in% names(geno), "`hap$genotype_wide` must contain a sample_id column.")
  geno <- geno[, c("sample_id", variant_ids), with = FALSE]
  geno[, "sample_id" := as.character(sample_id)]

  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  if (!identical(sample_col, "sample_id")) {
    data.table::setnames(pheno, sample_col, "sample_id")
  }
  pheno[, "sample_id" := as.character(sample_id)]

  dt <- merge(geno, pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between genotypes and phenotype table.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & !trait %in% variant_ids, trait]
  }
  traits <- as.character(traits)
  stop_if_not(length(traits) > 0L, "No numeric phenotype traits were selected.")
  stop_if_not(all(traits %in% names(dt)), "Some traits were not found in phenotype table.")

  geno_long <- data.table::melt(
    dt,
    id.vars = "sample_id",
    measure.vars = variant_ids,
    variable.name = "variant_id",
    value.name = "genotype_group",
    variable.factor = FALSE
  )
  geno_long[, "genotype_group" := as.character(genotype_group)]
  geno_long <- geno_long[!is.na(genotype_group) & genotype_group != ""]

  pheno_long <- data.table::melt(
    dt[, c("sample_id", traits), with = FALSE],
    id.vars = "sample_id",
    measure.vars = traits,
    variable.name = "trait",
    value.name = "value",
    variable.factor = FALSE
  )
  pheno_long[, "value" := suppressWarnings(as.numeric(value))]
  pheno_long <- pheno_long[!is.na(value)]

  long <- merge(geno_long, pheno_long, by = "sample_id", allow.cartesian = TRUE)
  stop_if_not(nrow(long) > 0L, "No non-missing genotype-phenotype records were available.")

  group_count <- long[, .(sample_n = .N), by = .(trait, variant_id, genotype_group)]
  keep_groups <- group_count[sample_n >= min_group_samples, .(trait, variant_id, genotype_group)]
  long <- merge(long, keep_groups, by = c("trait", "variant_id", "genotype_group"), all.x = FALSE, sort = FALSE)
  stop_if_not(nrow(long) > 0L, "No genotype groups remained after `min_group_samples` filtering.")

  effect_dt <- long[, calculate_one_variant_effect(.SD, test_method = test_method), by = .(trait, variant_id)]
  effect_dt <- effect_dt[group_n >= 2L]
  stop_if_not(nrow(effect_dt) > 0L, "At least two genotype groups with enough samples are required for variant-effect calculation.")
  effect_dt[, "p_adj" := stats::p.adjust(p_value, method = p_adjust), by = trait]

  vars <- data.table::as.data.table(hap$variants)
  if (!"variant_id" %in% names(vars)) {
    vars[, "variant_id" := variant_ids]
  }
  variant_label <- match.arg(variant_label)
  if (!variant_label %in% names(vars)) variant_label <- "variant_id"
  vars[, "variant_label" := as.character(.SD[[1L]]), .SDcols = variant_label]
  vars[is.na(variant_label) | variant_label == "", "variant_label" := as.character(variant_id)]
  if (!"pos" %in% names(vars)) vars[, "pos" := seq_len(.N)]
  vars[, "variant_index" := match(variant_id, variant_ids)]

  effect_dt <- merge(
    effect_dt,
    vars[, intersect(c("variant_id", "variant_label", "chrom", "pos", "type", "variant_index"), names(vars)), with = FALSE],
    by = "variant_id",
    all.x = TRUE,
    sort = FALSE
  )
  effect_dt[, "plot_effect" := if (identical(effect_type, "absolute")) abs_effect else effect]
  effect_dt[, "plot_log10_padj" := -log10(pmax(p_adj, .Machine$double.xmin))]
  finite_alpha <- effect_dt[is.finite(plot_log10_padj), max(plot_log10_padj, na.rm = TRUE)]
  if (!is.finite(finite_alpha)) finite_alpha <- 1
  effect_dt[!is.finite(plot_log10_padj), "plot_log10_padj" := finite_alpha]
  effect_dt[, "x_value" := if (identical(x_axis, "position")) as.numeric(pos) else as.numeric(variant_index)]
  effect_dt[is.na(x_value), "x_value" := seq_len(.N)]

  label_dt <- effect_dt[order(trait, -abs_effect)]
  if (top_n > 0L) {
    label_dt <- label_dt[, head(.SD, top_n), by = trait]
  } else {
    label_dt <- label_dt[0]
  }

  p <- ggplot2::ggplot(effect_dt, ggplot2::aes(x = .data$x_value, y = .data$plot_effect)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(
      size = .data$abs_effect,
      alpha = .data$plot_log10_padj,
      color = .data$effect
    )) +
    ggplot2::scale_color_gradient2(
      name = "Signed effect",
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      na.value = "grey70"
    ) +
    ggplot2::scale_size_continuous(name = "Absolute effect", range = c(point_size * 0.5, point_size * 1.8)) +
    ggplot2::scale_alpha_continuous(name = "-log10(adj. P)", range = c(0.45, 1), na.value = 0.45) +
    ggplot2::labs(
      x = if (identical(x_axis, "position")) "Genomic position" else "Variant index",
      y = if (identical(effect_type, "absolute")) "Absolute phenotype effect" else "Signed phenotype effect"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(size = text_size, color = "black"),
      axis.text.x = ggplot2::element_text(
        size = text_size, angle = 90, hjust = 1, vjust = 0.5, color = "black"
      ),
      axis.title = ggplot2::element_text(size = text_size, color = "black"),
      legend.title = ggplot2::element_text(size = text_size * 0.9, color = "black"),
      legend.text = ggplot2::element_text(size = text_size * 0.85, color = "black"),
      strip.text = ggplot2::element_text(size = text_size, color = "black"),
      strip.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (nrow(label_dt) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_dt,
      ggplot2::aes(label = .data$variant_label),
      inherit.aes = TRUE,
      size = label_size,
      vjust = -0.6,
      check_overlap = TRUE
    )
  }

  if (length(unique(effect_dt$trait)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free_y")
  }

  attr(p, "effect_table") <- effect_dt[]
  attr(p, "plot_data") <- long[]

  out <- list(
    figure = p,
    effect = effect_dt[],
    plot_data = long[]
  )
  class(out) <- c("GeneTrackRVariantEffectPlot", "list")
  out
}

#' @export
print.HapRefined <- function(x, ...) {
  cat("<HapRefined>\n")
  cat("  original haplotypes: ", format(length(unique(x$haplotype_map$hap_id)), big.mark = ","), "\n", sep = "")
  cat("  refined haplotypes : ", format(nrow(x$refined_haplotypes), big.mark = ","), "\n", sep = "")
  cat("  traits             : ", paste(x$parameters$traits, collapse = ", "), "\n", sep = "")
  invisible(x)
}

make_refined_groups_for_trait <- function(summary_dt,
                                          test_dt,
                                          alpha = 0.05,
                                          effect_threshold = NULL) {
  dt <- data.table::copy(summary_dt)
  if (nrow(dt) == 0L) return(data.table::data.table())
  data.table::setorder(dt, mean_value, hap_id)
  haps <- as.character(dt$hap_id)
  group_id <- character(length(haps))
  group_members <- list()
  group_n <- 0L

  for (i in seq_along(haps)) {
    hap_i <- haps[i]
    assigned <- FALSE
    if (group_n > 0L) {
      for (g in seq_len(group_n)) {
        members <- group_members[[g]]
        ok <- all(vapply(members, function(member) {
          can_merge_pair(hap_i, member, test_dt = test_dt, alpha = alpha, effect_threshold = effect_threshold)
        }, logical(1L)))
        if (isTRUE(ok)) {
          group_members[[g]] <- c(group_members[[g]], hap_i)
          group_id[i] <- paste0("G", g)
          assigned <- TRUE
          break
        }
      }
    }
    if (!assigned) {
      group_n <- group_n + 1L
      group_members[[group_n]] <- hap_i
      group_id[i] <- paste0("G", group_n)
    }
  }

  dt[, "trait_group" := group_id]
  dt[, .(trait, hap_id, trait_group, mean_value, median_value, sample_n)]
}

can_merge_pair <- function(hap1,
                           hap2,
                           test_dt,
                           alpha = 0.05,
                           effect_threshold = NULL) {
  if (identical(hap1, hap2)) return(TRUE)
  if (nrow(test_dt) == 0L) return(FALSE)
  row <- test_dt[(group1 == hap1 & group2 == hap2) | (group1 == hap2 & group2 == hap1)]
  if (nrow(row) == 0L) return(FALSE)
  row <- row[1L]
  ok <- !is.na(row$p_adj) && row$p_adj > alpha
  if (!isTRUE(ok)) return(FALSE)
  if (!is.null(effect_threshold)) {
    ok <- !is.na(row$abs_mean_diff) && row$abs_mean_diff <= effect_threshold
  }
  isTRUE(ok)
}

collapse_refined_variant_state <- function(x, mixed_label = "mixed") {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  if (length(x) == 1L) return(x[1L])
  mixed_label
}

get_hap_variant_columns <- function(hap) {
  vars <- data.table::as.data.table(hap$variants)
  variant_ids <- character()
  if ("variant_id" %in% names(vars)) {
    variant_ids <- as.character(vars$variant_id)
  }
  candidate_tables <- list(hap$genotype_wide, hap$haplotypes)
  for (tab in candidate_tables) {
    dt <- data.table::as.data.table(tab)
    cols <- intersect(variant_ids, names(dt))
    if (length(cols) > 0L) return(cols)
  }
  meta_cols <- c("sample_id", "hap_id", "hap_pattern", "non_missing_variant_n", "sample_n", "samples", "original_hap_n", "original_hap_ids")
  dt <- data.table::as.data.table(hap$genotype_wide)
  setdiff(names(dt), meta_cols)
}

calculate_one_variant_effect <- function(dt, test_method = "t.test") {
  dt <- data.table::as.data.table(dt)
  group_summary <- dt[, .(
    sample_n = .N,
    mean_value = mean(value, na.rm = TRUE),
    median_value = stats::median(value, na.rm = TRUE),
    sd_value = stats::sd(value, na.rm = TRUE)
  ), by = genotype_group]
  data.table::setorder(group_summary, genotype_group)
  group_n <- nrow(group_summary)
  if (group_n < 2L) {
    return(data.table::data.table())
  }

  min_idx <- which.min(group_summary$mean_value)
  max_idx <- which.max(group_summary$mean_value)
  low_group <- group_summary$genotype_group[min_idx]
  high_group <- group_summary$genotype_group[max_idx]
  effect <- group_summary$mean_value[max_idx] - group_summary$mean_value[min_idx]

  if (group_n == 2L) {
    effect <- group_summary$mean_value[2L] - group_summary$mean_value[1L]
    low_group <- group_summary$genotype_group[1L]
    high_group <- group_summary$genotype_group[2L]
  }

  p_value <- run_variant_effect_test(dt, method = test_method)

  data.table::data.table(
    group_n = group_n,
    sample_n = nrow(dt),
    effect = effect,
    abs_effect = abs(effect),
    low_group = as.character(low_group),
    high_group = as.character(high_group),
    low_group_mean = group_summary$mean_value[match(low_group, group_summary$genotype_group)],
    high_group_mean = group_summary$mean_value[match(high_group, group_summary$genotype_group)],
    p_value = p_value,
    test_method = test_method,
    group_summary = paste(
      paste0(group_summary$genotype_group, ":n=", group_summary$sample_n, ",mean=", signif(group_summary$mean_value, 4)),
      collapse = ";"
    )
  )
}

run_variant_effect_test <- function(dt, method = "t.test") {
  dt <- data.table::as.data.table(dt)
  groups <- unique(as.character(dt$genotype_group))
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2L) return(NA_real_)
  tryCatch({
    if (length(groups) == 2L && identical(method, "t.test")) {
      stats::t.test(value ~ genotype_group, data = dt)$p.value
    } else if (length(groups) == 2L && identical(method, "wilcox.test")) {
      stats::wilcox.test(value ~ genotype_group, data = dt)$p.value
    } else if (identical(method, "kruskal.test")) {
      stats::kruskal.test(value ~ genotype_group, data = dt)$p.value
    } else {
      stats::anova(stats::aov(value ~ genotype_group, data = dt))[["Pr(>F)"]][1L]
    }
  }, error = function(e) NA_real_)
}
