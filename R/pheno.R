# Author: Rensc
# Date: 2026-05-29
# Version: 0.1.0
# Function: Read, summarize, and plot phenotype tables
# Input: Phenotype files
# Output: Phenotype tables and ggplot figures

#' Read a phenotype table
#'
#' @description
#' Reads a phenotype table where the first column contains taxa/sample IDs and
#' all remaining columns are phenotype traits.
#'
#' @param file Phenotype file path.
#' @param sep Field separator. Use `auto` for automatic detection by `fread()`.
#' @param sample_col Optional sample column name. If NULL, the first column is used.
#' @param na_strings Strings treated as missing values.
#' @return A data.table with the first column standardized as `sample_id`.
#' @examples
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' pheno <- read_pheno(pheno_file)
#' head(pheno)
#' summary_pheno(pheno)
#' plot_pheno(pheno, traits = c("plant_height", "seed_weight"))
#' @export
read_pheno <- function(file,
                       sep = "auto",
                       sample_col = NULL,
                       na_strings = c("NA", "NaN", "", ".", "null", "NULL")) {
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))
  dt <- data.table::fread(
    file,
    sep = sep,
    na.strings = na_strings,
    data.table = TRUE,
    showProgress = FALSE
  )
  stop_if_not(ncol(dt) >= 2L, "Phenotype file must contain a sample column and at least one trait column.")
  if (is.null(sample_col)) {
    sample_col <- names(dt)[1L]
  }
  stop_if_not(sample_col %in% names(dt), paste0("Sample column was not found: ", sample_col))
  data.table::setnames(dt, sample_col, "sample_id")
  dt[, "sample_id" := as.character(sample_id)]
  stop_if_not(!anyDuplicated(dt$sample_id), "Duplicated sample IDs were found in phenotype table.")
  dt[]
}

#' Summarize phenotype traits
#'
#' @param pheno A phenotype data.frame/data.table returned by `read_pheno()`.
#' @param sample_col Sample column name.
#' @return A data.table summarizing trait type, missing values, and basic statistics.
#' @export
summary_pheno <- function(pheno, sample_col = "sample_id") {
  dt <- data.table::as.data.table(pheno)
  stop_if_not(sample_col %in% names(dt), paste0("Sample column was not found: ", sample_col))
  traits <- setdiff(names(dt), sample_col)
  out <- lapply(traits, function(trait) {
    x <- dt[[trait]]
    is_num <- is.numeric(x) || is.integer(x)
    missing_n <- sum(is.na(x))
    non_missing <- x[!is.na(x)]
    data.table::data.table(
      trait = trait,
      type = if (is_num) "numeric" else "categorical",
      sample_n = nrow(dt),
      missing_n = missing_n,
      missing_rate = missing_n / nrow(dt),
      unique_n = length(unique(non_missing)),
      min = if (is_num && length(non_missing) > 0L) suppressWarnings(min(non_missing)) else NA_real_,
      mean = if (is_num && length(non_missing) > 0L) suppressWarnings(mean(non_missing)) else NA_real_,
      median = if (is_num && length(non_missing) > 0L) suppressWarnings(stats::median(non_missing)) else NA_real_,
      max = if (is_num && length(non_missing) > 0L) suppressWarnings(max(non_missing)) else NA_real_
    )
  })
  data.table::rbindlist(out, fill = TRUE)
}

#' Plot phenotype distributions
#'
#' @param pheno A phenotype data.frame/data.table.
#' @param traits Optional trait names. If NULL, all traits are plotted.
#' @param sample_col Sample column name.
#' @param bins Histogram bins for numeric traits.
#' @param text_size Text size.
#' @return A ggplot object.
#' @export
plot_pheno <- function(pheno,
                       traits = NULL,
                       sample_col = "sample_id",
                       bins = 30L,
                       text_size = 14) {
  dt <- data.table::as.data.table(pheno)
  stop_if_not(sample_col %in% names(dt), paste0("Sample column was not found: ", sample_col))
  if (is.null(traits)) traits <- setdiff(names(dt), sample_col)
  traits <- as.character(traits)
  stop_if_not(all(traits %in% names(dt)), "Some traits were not found in phenotype table.")

  info <- summary_pheno(dt, sample_col = sample_col)
  numeric_traits <- info[trait %in% traits & type == "numeric", trait]
  categorical_traits <- info[trait %in% traits & type != "numeric", trait]

  plots <- list()

  if (length(numeric_traits) > 0L) {
    x <- data.table::melt(
      dt,
      id.vars = sample_col,
      measure.vars = numeric_traits,
      variable.name = "trait",
      value.name = "value",
      variable.factor = FALSE
    )
    plots$numeric <- ggplot2::ggplot(x, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(bins = as.integer(bins), na.rm = TRUE) +
      ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free") +
      ggplot2::labs(x = "Phenotype value", y = "Sample count") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = text_size),
        axis.title = ggplot2::element_text(size = text_size),
        strip.text = ggplot2::element_text(size = text_size)
      )
  }

  if (length(categorical_traits) > 0L) {
    x <- data.table::melt(
      dt,
      id.vars = sample_col,
      measure.vars = categorical_traits,
      variable.name = "trait",
      value.name = "value",
      variable.factor = FALSE
    )
    plots$categorical <- ggplot2::ggplot(x, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_bar(na.rm = TRUE) +
      ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free") +
      ggplot2::labs(x = "Phenotype class", y = "Sample count") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = text_size, angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(size = text_size),
        axis.title = ggplot2::element_text(size = text_size),
        strip.text = ggplot2::element_text(size = text_size)
      )
  }

  if (length(plots) == 0L) stop("No traits were available for plotting.", call. = FALSE)
  if (length(plots) == 1L) return(plots[[1L]])
  patchwork::wrap_plots(plots, ncol = 1)
}
