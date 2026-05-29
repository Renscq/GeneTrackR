# Author: Rensc
# Date: 2026-05-29
# Version: 0.1.0
# Function: Plot haplotype-variant and haplotype-phenotype figures
# Input: HapVariant and phenotype objects
# Output: ggplot or patchwork figures

#' Plot variants and haplotype table
#'
#' @description
#' Draws a gene/region model, variant positions, and a haplotype genotype table.
#' Variant columns are aligned by genomic position.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param annotation Optional annotation object used to draw gene models.
#' @param show_gene_model Logical. Whether to draw a gene model panel.
#' @param text_size Text size.
#' @param table_text_size Haplotype table text size.
#' @return A patchwork object.
#' @export
plot_hap_variant <- function(hap,
                             annotation = NULL,
                             show_gene_model = !is.null(annotation),
                             text_size = 14,
                             table_text_size = 3.2) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")

  region <- hap$region
  vars <- data.table::as.data.table(hap$variants)
  haps <- data.table::as.data.table(hap$haplotypes)
  variant_ids <- vars[["variant_id"]]

  p_list <- list()
  heights <- numeric()

  if (isTRUE(show_gene_model) && !is.null(annotation)) {
    p_gene <- plot_region(
      annotation,
      chrom = region$chrom,
      start = region$start,
      end = region$end,
      mode = "overlap",
      collapse = "none",
      label_position = "axis",
      text_size = text_size
    )
    p_list$gene <- p_gene
    heights <- c(heights, 1.0)
  }

  p_var <- plot_variant(
    VariantTrack(vars, meta = list(format = "hap_variant")),
    chrom = region$chrom,
    start = region$start,
    end = region$end,
    color_by = "variant_type",
    label_by = "none",
    text_size = text_size,
    track_name = "Variant"
  )
  p_list$variant <- p_var
  heights <- c(heights, 0.8)

  table_long <- data.table::melt(
    haps,
    id.vars = intersect(c("hap_id", "sample_n", "samples"), names(haps)),
    measure.vars = variant_ids,
    variable.name = "variant_id",
    value.name = "genotype",
    variable.factor = FALSE
  )
  table_long <- merge(
    table_long,
    vars[, .(variant_id, pos)],
    by = "variant_id",
    all.x = TRUE
  )
  table_long[, "hap_label" := paste0(hap_id, " (n=", sample_n, ")")]
  table_long[, "hap_label" := factor(hap_label, levels = rev(unique(hap_label)))]

  p_table <- ggplot2::ggplot(table_long, ggplot2::aes(x = .data$pos, y = .data$hap_label, fill = .data$genotype)) +
    ggplot2::geom_tile(color = "grey85", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$genotype), size = table_text_size) +
    ggplot2::coord_cartesian(xlim = c(region$start, region$end)) +
    ggplot2::labs(x = paste0("Chromosome ", region$chrom, " position (bp)"), y = "Haplotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = text_size),
      axis.title = ggplot2::element_text(size = text_size),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = text_size),
      panel.grid = ggplot2::element_blank()
    )

  p_list$haplotype <- p_table
  heights <- c(heights, max(1.2, 0.25 * nrow(haps)))

  patchwork::wrap_plots(p_list, ncol = 1) + patchwork::plot_layout(heights = heights)
}

#' Plot phenotype values grouped by haplotype
#'
#' @description
#' Merges haplotype assignments and phenotype values, then draws boxplots with
#' pairwise t-test summaries.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_hap_samples Minimum sample number required for a haplotype group.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param show_points Logical. Whether to show sample points.
#' @param text_size Text size.
#' @return A ggplot object with attributes `plot_data` and `test_table`.
#' @export
plot_hap_pheno <- function(hap,
                           phenotype,
                           traits = NULL,
                           sample_col = "sample_id",
                           min_hap_samples = 2L,
                           p_adjust = "BH",
                           show_points = TRUE,
                           text_size = 14) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  data.table::setnames(pheno, sample_col, "sample_id")

  hap_sample <- data.table::as.data.table(hap$sample_haplotypes)
  dt <- merge(hap_sample[, .(sample_id, hap_id)], pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between haplotypes and phenotype table.")

  hap_keep <- dt[, .N, by = hap_id][N >= as.integer(min_hap_samples), hap_id]
  dt <- dt[hap_id %in% hap_keep]
  stop_if_not(length(unique(dt$hap_id)) >= 2L, "At least two haplotypes with enough samples are required.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & trait != "hap_id", trait]
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
  long <- long[!is.na(value)]
  stop_if_not(nrow(long) > 0L, "No non-missing phenotype values were available.")

  test_dt <- pairwise_hap_ttest(long, p_adjust = p_adjust)
  label_dt <- make_hap_test_labels(long, test_dt)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$hap_id, y = .data$value, fill = .data$hap_id)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free_y") +
    ggplot2::labs(x = "Haplotype", y = "Phenotype value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = text_size, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = text_size),
      axis.title = ggplot2::element_text(size = text_size),
      strip.text = ggplot2::element_text(size = text_size),
      legend.position = "none"
    )

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.75)
  }

  if (nrow(label_dt) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_dt,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      size = text_size / 3.2,
      hjust = 0.5,
      vjust = 1
    )
  }

  attr(p, "plot_data") <- long[]
  attr(p, "test_table") <- test_dt[]
  p
}

pairwise_hap_ttest <- function(long, p_adjust = "BH") {
  out <- long[, {
    haps <- unique(hap_id)
    if (length(haps) < 2L) {
      return(data.table::data.table(group1 = character(), group2 = character(), p_value = numeric()))
    }
    pairs <- utils::combn(haps, 2L, simplify = FALSE)
    res <- lapply(pairs, function(pair) {
      x <- value[hap_id == pair[1L]]
      y <- value[hap_id == pair[2L]]
      if (length(x) < 2L || length(y) < 2L) {
        pval <- NA_real_
      } else {
        pval <- tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
      }
      data.table::data.table(group1 = pair[1L], group2 = pair[2L], p_value = pval)
    })
    data.table::rbindlist(res)
  }, by = trait]
  out[, "p_adj" := stats::p.adjust(p_value, method = p_adjust), by = trait]
  out[, "label" := paste0(group1, " vs ", group2, ": p=", format_p_value(p_adj))]
  out[]
}

make_hap_test_labels <- function(long, test_dt) {
  if (nrow(test_dt) == 0L) return(data.table::data.table())
  best <- test_dt[!is.na(p_adj)]
  if (nrow(best) == 0L) return(data.table::data.table())
  data.table::setorder(best, trait, p_adj)
  best <- best[, .SD[1L], by = trait]
  y_dt <- long[, .(y = max(value, na.rm = TRUE)), by = trait]
  y_dt[, "y" := y + abs(y) * 0.08 + 1e-8]
  x_dt <- long[, .(x = mean(seq_along(unique(hap_id)))), by = trait]
  merge(best[, .(trait, label)], merge(y_dt, x_dt, by = "trait"), by = "trait")
}

format_p_value <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3g", p)))
}
