# Author: Rensc
# Date: 2026-05-31
# Version: 0.1.0
# Function: Plot phenotype distributions grouped by a single variant genotype
# Input: VariantTrack/VCF/HapVariant objects and phenotype tables
# Output: ggplot figures with pairwise test annotations

#' Plot phenotype values grouped by a single variant genotype
#'
#' @description
#' Draws phenotype distributions for genotype groups at one variant site. The
#' function is designed as the single-variant counterpart of `plot_hap_pheno()`.
#' Groups are ordered by sample number from left to right, and fill colors are
#' mapped to the median phenotype value of each genotype group.
#'
#' @param variant A VariantTrack object, a VCF file path, a HapVariant object, or
#' a VCF-like data.frame/data.table containing one or more variants with genotype
#' sample columns.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param variant_id Optional variant ID to select.
#' @param chrom Optional chromosome name used to select a variant.
#' @param pos Optional 1-based variant position. This is equivalent to setting
#' both `start` and `end` to the same value.
#' @param start Optional 1-based start position for region-based selection.
#' @param end Optional 1-based end position for region-based selection.
#' @param samples Optional sample names to keep.
#' @param genotype_mode Genotype representation. `code` converts genotypes to 0/1
#' states, where 0 means reference genotype and 1 means any alternate allele is
#' present. `string` converts genotypes to compact allele labels; long InDel
#' alleles are compressed as `iN`, where `N` is allele length.
#' @param missing_genotype Missing genotype label. Default is `NA_character_`.
#' @param min_group_samples Minimum sample number required for a genotype group.
#' @param plot_type Plot type. One of `violin`, `boxplot`, or `violin_boxplot`.
#' @param test_method Pairwise test method. One of `t.test`, `wilcox.test`, or `ks.test`.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param p_label P-value label style. One of `stars`, `number`, or `both`.
#' @param p_cutoff Significance cutoff used for displaying pairwise comparisons.
#' @param p_value_type Which p-value is used for filtering and labeling. One of `raw` or `adjusted`.
#' @param show_signif_only Logical. Whether to only display significant comparisons.
#' @param show_points Logical. Whether to show sample points.
#' @param show_outliers Logical. Whether to show boxplot outliers.
#' @param fill_palette RColorBrewer palette name used for median-based genotype fills.
#' @param fill_colors Optional custom fill colors for genotype groups.
#' @param fill_alpha Alpha value for violin/boxplot fill colors.
#' @param violin_width Violin plot width.
#' @param box_width Boxplot width.
#' @param bracket_step Fraction of y-range used to separate significance brackets.
#' @param bracket_tip_fraction Fraction of bracket vertical spacing used for the short downward bracket tips.
#' @param x_text_angle Rotation angle for genotype labels on the x-axis.
#' @param strip_label_width Maximum character width for wrapping long facet strip labels.
#' @param strip_fill Strip background fill color. Default is white.
#' @param strip_border_color Strip border color. Default NULL removes the strip border.
#' @param text_size Text size.
#' @return A list with `figure` and `pvalue` elements. Additional elements include
#' `summary`, `bracket`, `plot_data`, and `variant_data`.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' pheno <- read_pheno(pheno_file)
#' plot_variant_pheno(vcf, phenotype = pheno, variant_id = "var1",
#'                    traits = "plant_height", min_group_samples = 1)
#' plot_variant_pheno(vcf, phenotype = pheno, chrom = "chr1", pos = 120,
#'                    traits = "seed_weight", genotype_mode = "string",
#'                    min_group_samples = 1, test_method = "wilcox.test")
#' @export
plot_variant_pheno <- function(variant,
                               phenotype,
                               traits = NULL,
                               sample_col = "sample_id",
                               variant_id = NULL,
                               chrom = NULL,
                               pos = NULL,
                               start = NULL,
                               end = NULL,
                               samples = NULL,
                               genotype_mode = c("code", "string"),
                               missing_genotype = NA_character_,
                               min_group_samples = 2L,
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
                               strip_fill = "white",
                               strip_border_color = NULL,
                               strip_text_lineheight = 0.9,
                               text_size = 14) {
  genotype_mode <- match.arg(genotype_mode)
  plot_type <- match.arg(plot_type)
  test_method <- match.arg(test_method)
  p_label <- match.arg(p_label)
  p_value_type <- match.arg(p_value_type)

  min_group_samples <- as.integer(min_group_samples)[1L]
  stop_if_not(!is.na(min_group_samples) && min_group_samples >= 1L, "`min_group_samples` must be a positive integer.")
  p_cutoff <- as.numeric(p_cutoff)[1L]
  if (is.na(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) {
    p_cutoff <- 0.05
  }
  fill_alpha <- as.numeric(fill_alpha)[1L]
  if (is.na(fill_alpha)) fill_alpha <- 0.75
  fill_alpha <- max(0, min(1, fill_alpha))
  x_text_angle <- as.numeric(x_text_angle)[1L]
  if (is.na(x_text_angle)) x_text_angle <- 90
  strip_label_width <- as.integer(strip_label_width)[1L]
  if (is.na(strip_label_width) || strip_label_width < 1L) strip_label_width <- 24L
  strip_text_lineheight <- as.numeric(strip_text_lineheight)[1L]
  if (is.na(strip_text_lineheight) || strip_text_lineheight <= 0) strip_text_lineheight <- 0.9
  strip_text_lineheight <- 0.9
  strip_fill <- as.character(strip_fill)[1L]
  if (is.na(strip_fill) || strip_fill == "") strip_fill <- "white"
  if (is.null(strip_border_color) || length(strip_border_color) == 0L || is.na(strip_border_color[1L]) || strip_border_color[1L] == "") {
    strip_border_color <- NA_character_
  } else {
    strip_border_color <- as.character(strip_border_color)[1L]
  }

  variant_info <- extract_single_variant_genotype(
    variant = variant,
    variant_id = variant_id,
    chrom = chrom,
    pos = pos,
    start = start,
    end = end,
    samples = samples,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype
  )

  geno <- variant_info$genotype
  stop_if_not(nrow(geno) > 0L, "No genotype records were available for the selected variant.")
  geno <- geno[!is.na(genotype_group)]
  stop_if_not(nrow(geno) > 0L, "No non-missing genotypes were available for the selected variant.")

  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  if (!identical(sample_col, "sample_id")) {
    data.table::setnames(pheno, sample_col, "sample_id")
  }
  pheno[, "sample_id" := as.character(sample_id)]

  dt <- merge(geno[, .(sample_id, genotype_group)], pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between variant genotypes and phenotype table.")

  group_count <- dt[, .(sample_n = .N), by = genotype_group]
  group_keep <- group_count[sample_n >= min_group_samples, genotype_group]
  dt <- dt[genotype_group %in% group_keep]
  group_count <- group_count[genotype_group %in% group_keep]
  stop_if_not(length(unique(dt$genotype_group)) >= 2L, "At least two genotype groups with enough samples are required.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & !trait %in% c("genotype_group", "sample_n"), trait]
  }
  traits <- as.character(traits)
  stop_if_not(length(traits) > 0L, "No numeric phenotype traits were selected.")
  stop_if_not(all(traits %in% names(dt)), "Some traits were not found in phenotype table.")

  long <- data.table::melt(
    dt,
    id.vars = c("sample_id", "genotype_group"),
    measure.vars = traits,
    variable.name = "trait",
    value.name = "value",
    variable.factor = FALSE
  )
  long[, "value" := suppressWarnings(as.numeric(value))]
  long <- long[!is.na(value)]
  stop_if_not(nrow(long) > 0L, "No non-missing phenotype values were available.")

  group_count <- group_count[order(-sample_n, genotype_group)]
  group_order <- group_count$genotype_group
  group_axis_labels <- paste0(as.character(group_count$genotype_group), " (", group_count$sample_n, ")")
  names(group_axis_labels) <- as.character(group_count$genotype_group)

  long[, "hap_id" := factor(genotype_group, levels = group_order)]

  summary_dt <- long[, .(
    sample_n = .N,
    median_value = stats::median(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ), by = .(trait, hap_id)]

  fill_dt <- make_hap_pheno_fill_table(summary_dt, fill_palette = fill_palette, fill_colors = fill_colors, alpha = fill_alpha)
  long <- merge(long, fill_dt[, .(trait, hap_id, hap_fill)], by = c("trait", "hap_id"), all.x = TRUE, sort = FALSE)

  test_dt <- pairwise_hap_test(long, method = test_method, p_adjust = p_adjust)
  bracket_dt <- make_hap_brackets(
    long = long,
    test_dt = test_dt,
    p_cutoff = p_cutoff,
    p_label = p_label,
    p_value_type = p_value_type,
    show_signif_only = show_signif_only,
    bracket_step = bracket_step,
    bracket_tip_fraction = bracket_tip_fraction
  )


  variant_title <- make_variant_pheno_title(variant_info$variant)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$hap_id, y = .data$value, fill = .data$hap_fill))

  if (plot_type %in% c("violin", "violin_boxplot")) {
    p <- p + ggplot2::geom_violin(
      width = violin_width,
      trim = FALSE,
      color = "black",
      linewidth = 0.35,
      alpha = fill_alpha
    )
  }

  if (plot_type %in% c("boxplot", "violin_boxplot")) {
    outlier_shape <- if (isTRUE(show_outliers)) 19 else NA
    p <- p + ggplot2::geom_boxplot(
      width = box_width,
      outlier.shape = outlier_shape,
      color = "black",
      linewidth = 0.35,
      alpha = min(1, fill_alpha + 0.15)
    )
  }

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      size = 1.4,
      alpha = 0.65,
      color = "black"
    )
  }

  if (length(traits) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.data$trait),
      scales = "free_y",
      labeller = ggplot2::as_labeller(function(x) wrap_strip_labels(x, width = strip_label_width))
    )
  }

  if (nrow(bracket_dt) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        data = bracket_dt,
        ggplot2::aes(x = .data$x_mid, y = .data$label_y, label = .data$label),
        inherit.aes = FALSE,
        color = "black",
        size = text_size / 3.2,
        vjust = 0
      )
  }

  p <- p +
    ggplot2::scale_x_discrete(labels = group_axis_labels) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = "Genotype", y = "Phenotype value", title = variant_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(size = text_size, color = "black", hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = text_size, angle = x_text_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title = ggplot2::element_text(size = text_size, color = "black"),
      strip.text = ggplot2::element_text(size = text_size, color = "black", lineheight = strip_text_lineheight, margin = ggplot2::margin(3, 3, 3, 3)),
      strip.background = ggplot2::element_rect(fill = strip_fill, color = strip_border_color),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none"
    )

  attr(p, "plot_data") <- long[]
  attr(p, "summary_table") <- summary_dt[]
  attr(p, "test_table") <- test_dt[]
  attr(p, "bracket_table") <- bracket_dt[]
  attr(p, "variant_data") <- variant_info$variant[]

  out <- list(
    figure = p,
    pvalue = test_dt[],
    summary = summary_dt[],
    bracket = bracket_dt[],
    plot_data = long[],
    variant_data = variant_info$variant[]
  )
  class(out) <- c("GeneTrackRPhenoPlot", "list")
  out
}

extract_single_variant_genotype <- function(variant,
                                            variant_id = NULL,
                                            chrom = NULL,
                                            pos = NULL,
                                            start = NULL,
                                            end = NULL,
                                            samples = NULL,
                                            genotype_mode = "code",
                                            missing_genotype = NA_character_) {
  if (inherits(variant, "HapVariant")) {
    return(extract_single_variant_from_hap(
      hap = variant,
      variant_id = variant_id,
      chrom = chrom,
      pos = pos,
      start = start,
      end = end,
      samples = samples
    ))
  }

  vt <- coerce_single_variant_input(
    variant = variant,
    variant_id = variant_id,
    chrom = chrom,
    pos = pos,
    start = start,
    end = end
  )

  sample_cols <- get_vcf_sample_columns(vt)
  stop_if_not(length(sample_cols) > 0L, "No VCF sample genotype columns were found.")
  if (!is.null(samples)) {
    samples <- as.character(samples)
    missing_samples <- setdiff(samples, sample_cols)
    stop_if_not(length(missing_samples) == 0L, paste0("Samples not found in VCF: ", paste(missing_samples, collapse = ", ")))
    sample_cols <- samples
  }

  stop_if_not(nrow(vt$data) == 1L, "Exactly one variant must be selected for `plot_variant_pheno()`.")
  geno <- extract_vcf_genotype_long(
    vt$data,
    sample_cols = sample_cols,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype
  )
  geno[, "genotype_group" := as.character(genotype)]
  geno[is.na(genotype_group), "genotype_group" := NA_character_]
  list(variant = vt$data[], genotype = geno[, .(sample_id, genotype_group)])
}

coerce_single_variant_input <- function(variant,
                                        variant_id = NULL,
                                        chrom = NULL,
                                        pos = NULL,
                                        start = NULL,
                                        end = NULL) {
  if (is.character(variant) && length(variant) == 1L && file.exists(variant)) {
    if (!is.null(pos)) {
      start <- as.integer(pos)[1L]
      end <- as.integer(pos)[1L]
    }
    if (!is.null(chrom) && !is.null(start) && !is.null(end)) {
      vt <- retrieve_vcf(variant, chrom = chrom, start = start, end = end, as = "VariantTrack")
    } else {
      vt <- read_vcf(variant, keep_genotype = TRUE, verbose = FALSE)
    }
  } else if (inherits(variant, "VariantTrack")) {
    vt <- variant
  } else if (is.data.frame(variant) || data.table::is.data.table(variant)) {
    vt <- VariantTrack(variant)
  } else {
    stop("`variant` must be a VariantTrack object, a VCF file path, a HapVariant object, or a VCF-like data.frame.", call. = FALSE)
  }

  dt <- data.table::as.data.table(vt$data)
  if (!is.null(variant_id)) {
    query_variant_id <- as.character(variant_id)
    dt <- dt[as.character(dt[["variant_id"]]) %in% query_variant_id]
  }
  if (!is.null(pos)) {
    start <- as.integer(pos)[1L]
    end <- as.integer(pos)[1L]
  }
  if (!is.null(chrom)) {
    dt <- dt[as.character(dt[["chrom"]]) == as.character(chrom)[1L]]
  }
  if (!is.null(start) && !is.null(end)) {
    dt <- dt[as.integer(dt[["pos"]]) >= as.integer(start)[1L] & as.integer(dt[["pos"]]) <= as.integer(end)[1L]]
  }
  stop_if_not(nrow(dt) > 0L, "No variant matched the selected locator.")
  stop_if_not(nrow(dt) == 1L, "The selected locator matched multiple variants. Please provide `variant_id` or exact `chrom` + `pos`.")

  out <- vt
  out$data <- dt[]
  out
}

extract_single_variant_from_hap <- function(hap,
                                            variant_id = NULL,
                                            chrom = NULL,
                                            pos = NULL,
                                            start = NULL,
                                            end = NULL,
                                            samples = NULL) {
  vars <- data.table::as.data.table(hap$variants)
  geno <- data.table::as.data.table(hap$genotype_long)
  if (!is.null(variant_id)) {
    query_variant_id <- as.character(variant_id)
    vars <- vars[as.character(vars[["variant_id"]]) %in% query_variant_id]
  }
  if (!is.null(pos)) {
    start <- as.integer(pos)[1L]
    end <- as.integer(pos)[1L]
  }
  if (!is.null(chrom)) {
    query_chrom <- as.character(chrom)[1L]
    vars <- vars[as.character(vars[["chrom"]]) == query_chrom]
  }
  if (!is.null(start) && !is.null(end)) {
    vars <- vars[as.integer(vars[["pos"]]) >= as.integer(start)[1L] & as.integer(vars[["pos"]]) <= as.integer(end)[1L]]
  }
  stop_if_not(nrow(vars) > 0L, "No variant matched the selected locator in the HapVariant object.")
  stop_if_not(nrow(vars) == 1L, "The selected locator matched multiple variants. Please provide `variant_id` or exact `chrom` + `pos`.")

  selected_id <- as.character(vars$variant_id[1L])
  geno <- geno[as.character(geno[["variant_id"]]) == selected_id]
  if (!is.null(samples)) {
    geno <- geno[as.character(sample_id) %in% as.character(samples)]
  }
  stop_if_not(nrow(geno) > 0L, "No genotype records were found for the selected variant.")
  geno[, "genotype_group" := as.character(genotype)]
  list(variant = vars[], genotype = geno[, .(sample_id, genotype_group)])
}

make_variant_pheno_title <- function(variant_dt) {
  if (is.null(variant_dt) || nrow(variant_dt) == 0L) return(NULL)
  x <- variant_dt[1L]
  paste0(
    as.character(x$variant_id), " (",
    as.character(x$chrom), ":", as.integer(x$pos), ", ",
    format_hap_ref_alt(x$ref, x$alt), ")"
  )
}
