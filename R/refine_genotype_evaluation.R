# Author: Rensc
# Date: 2026-11-09
# Version: dev001
# Function: Rank observed genotype combinations and evaluate frozen refinements
# Input: Aligned genotypes, phenotypes, covariates, and optional kinship
# Output: Combination ranks and paired holdout prediction summaries

# The cohort is fixed across subsets whenever it supports estimable candidates.
# Carrier IDs cover all supplied genotypes; model carrier IDs cover the cohort.
refine_combination_rank <- function(dosage, y, covariates, kinship, direction,
                                    min_n, max_order, sample_ids) {
  variants <- as.character(colnames(dosage))
  direction <- match.arg(direction, c("higher", "lower"))
  if (is.null(covariates)) {
    covariates <- matrix(numeric(), nrow = length(y), ncol = 0L)
  }
  stop_if_not(nrow(dosage) == length(y) && nrow(covariates) == length(y) &&
                length(sample_ids) == length(y), "Combination inputs must have aligned rows.")
  stop_if_not(!anyNA(sample_ids) && !anyDuplicated(sample_ids),
              "Combination sample IDs must be nonmissing and unique.")
  stop_if_not(ncol(dosage) <= 10L && max_order %in% 1:3,
              "Combination search supports up to ten variants and orders one to three.")
  stop_if_not(length(variants) == ncol(dosage) && !anyNA(variants) && !anyDuplicated(variants),
              "Combination variants must have nonmissing, unique column IDs.")
  stop_if_not(all(is.na(dosage) | dosage %in% 0:2),
              "Combination dosages must be 0, 1, 2, or missing.")

  # Stable row and variant order makes labels, ties, and alias choices reproducible.
  sample_order <- order(sample_ids, method = "radix")
  dosage <- dosage[sample_order, order(variants, method = "radix"), drop = FALSE]
  variants <- as.character(colnames(dosage))
  y <- y[sample_order]
  sample_ids <- sample_ids[sample_order]
  covariates <- covariates[sample_order, , drop = FALSE]
  if (!is.null(kinship)) kinship <- kinship[sample_order, sample_order, drop = FALSE]
  eligible <- is.finite(y) & rowSums(!is.finite(covariates)) == 0L
  common_rows <- which(eligible & rowSums(is.na(dosage)) == 0L)
  genotype <- matrix(c("0/0", "0/1", "1/1")[as.integer(dosage) + 1L], nrow(dosage), ncol(dosage))
  subsets <- unlist(lapply(seq_len(min(max_order, ncol(dosage))), function(k) {
    utils::combn(seq_along(variants), k, simplify = FALSE)
  }), recursive = FALSE)

  fit <- function(rows, state) {
    k <- if (is.null(kinship)) NULL else kinship[rows, rows, drop = FALSE]
    refine_model(y[rows], state, covariates[rows, , drop = FALSE], k)
  }

  # A common-cohort single-locus fit is reused for every multi-locus subset.
  single_cohorts <- list()
  single_results <- list()
  best_single <- function(rows) {
    if (!is.null(kinship)) return(list(variant = NA_character_, r2 = NA_real_))
    cached <- which(vapply(single_cohorts, identical, logical(1L), rows))
    if (length(cached)) return(single_results[[cached[1L]]])

    r2 <- vapply(seq_along(variants), function(j) {
      state <- genotype[rows, j]
      counts <- table(state)
      if (anyNA(state) || length(counts) < 2L ||
          !any(counts >= min_n & length(rows) - counts >= min_n)) return(NA_real_)
      fit(rows, state)$partial_r2
    }, numeric(1L))
    result <- list(variant = NA_character_, r2 = NA_real_)
    if (any(is.finite(r2))) {
      top <- max(r2, na.rm = TRUE)
      best <- which(abs(r2 - top) <= 1e-12)[1L]
      result <- list(variant = variants[best], r2 = r2[best])
    }
    single_cohorts[[length(single_cohorts) + 1L]] <<- rows
    single_results[[length(single_results) + 1L]] <<- result
    result
  }

  rank_subsets <- function(common) {
    result <- lapply(subsets, function(cols) {
      states <- do.call(paste, c(lapply(cols, function(col) genotype[, col]), sep = "|"))
      states[rowSums(is.na(dosage[, cols, drop = FALSE])) > 0L] <- NA_character_
      rows <- if (common) common_rows else which(eligible & !is.na(states))
      counts <- table(states[rows])
      states_to_rank <- names(counts)[counts >= min_n & length(rows) - counts >= min_n]
      if (!length(states_to_rank)) return(NULL)

      model <- fit(rows, states[rows])
      weights <- stats::setNames(as.numeric(counts) / length(rows), names(counts))
      best <- if (length(cols) > 1L && is.finite(model$partial_r2)) best_single(rows) else
        list(variant = NA_character_, r2 = NA_real_)
      subset <- paste(variants[cols], collapse = ";")
      comparison_key <- paste(rows, collapse = ",")
      lapply(states_to_rank, function(state) {
        contrast <- refine_difference(model, state, weights)
        if (!all(is.finite(contrast[c("estimate", "lower", "upper")]))) return(NULL)

        carrier <- which(states == state)
        model_carrier <- rows[states[rows] == state]
        alleles <- strsplit(state, "|", fixed = TRUE)[[1L]]
        combination <- paste(paste(variants[cols], alleles, sep = "="), collapse = ";")
        score <- if (direction == "higher") contrast[["lower"]] else -contrast[["upper"]]
        data.frame(
          combination = combination, subset = subset, state = state, order = length(cols),
          comparable = common, exploratory = TRUE, comparison_n = length(rows), sample_n = length(model_carrier),
          carrier_n = length(carrier), phenotype_n = sum(is.finite(y[carrier])),
          adjusted_mean = unname(model$means[state]),
          estimate = contrast[["estimate"]], lower = contrast[["lower"]], upper = contrast[["upper"]],
          score = score, partial_r2 = model$partial_r2, best_single = best$variant,
          single_delta_r2 = model$partial_r2 - best$r2,
          carrier_ids = I(list(sample_ids[carrier])), model_carrier_ids = I(list(sample_ids[model_carrier])),
          comparison_ids = I(list(sample_ids[rows])),
          carrier_key = paste(comparison_key, paste(carrier, collapse = ","), sep = ":")
        )
      })
    })
    do.call(rbind, unlist(result, recursive = FALSE))
  }

  common <- length(common_rows) >= 2L * min_n
  result <- if (length(subsets)) rank_subsets(common) else NULL
  if (is.null(result) && common && length(common_rows) < sum(eligible)) {
    common <- FALSE
    result <- rank_subsets(common)
  }
  if (is.null(result)) {
    result <- data.frame(
      combination = character(), subset = character(), state = character(), order = integer(),
      comparable = logical(), exploratory = logical(), comparison_n = integer(), sample_n = integer(),
      carrier_n = integer(), phenotype_n = integer(),
      adjusted_mean = numeric(), estimate = numeric(), lower = numeric(), upper = numeric(),
      score = numeric(), partial_r2 = numeric(), best_single = character(), single_delta_r2 = numeric(),
      carrier_ids = I(list()), model_carrier_ids = I(list()), comparison_ids = I(list()),
      aliases = I(list()), rank = integer(), subset_rank = integer()
    )
    status <- if (length(subsets)) "no_estimable_candidates" else "no_candidates"
  } else {
    result <- result[order(result$order, result$combination, method = "radix"), ]
    aliases <- split(result$combination, match(result$carrier_key, unique(result$carrier_key)))
    result <- result[!duplicated(result$carrier_key), ]
    result$aliases <- I(lapply(aliases, function(x) x[-1L]))
    result$carrier_key <- NULL
    sort_subset <- if (common) rep("", nrow(result)) else result$subset
    result <- result[order(sort_subset, -result$score, -result$sample_n, result$combination, method = "radix"), ]
    result$rank <- if (common) seq_len(nrow(result)) else rep(NA_integer_, nrow(result))
    result$subset_rank <- as.integer(stats::ave(seq_len(nrow(result)), result$subset, FUN = seq_along))
    rownames(result) <- NULL
    status <- if (common) "common_cohort" else "subset_cohorts"
  }
  attr(result, "comparison_status") <- status
  attr(result, "common_n") <- length(common_rows)
  result
}

# Fit both models on the same training rows. Holdout predictions use their
# frozen group assignments and fixed effects, including when a GRM is supplied.
refine_holdout <- function(y, pattern, refined_group, covariates, holdout,
                           kinship = NULL, min_coverage = 0.8, max_rmse_increase = NULL) {
  n <- length(y)
  if (is.null(covariates)) covariates <- matrix(numeric(), nrow = n, ncol = 0L)
  stop_if_not(is.logical(holdout) && !anyNA(holdout) && length(holdout) == n &&
                length(pattern) == n && length(refined_group) == n && nrow(covariates) == n,
              "Holdout inputs must have aligned rows and a complete logical holdout mask.")
  stop_if_not(length(min_coverage) == 1L && is.finite(min_coverage) && min_coverage >= 0 && min_coverage <= 1,
              "Minimum holdout coverage must be between zero and one.")
  stop_if_not(is.null(max_rmse_increase) || (length(max_rmse_increase) == 1L &&
                is.finite(max_rmse_increase) && max_rmse_increase >= 0),
              "Maximum holdout RMSE increase must be nonnegative or NULL.")
  pattern <- as.character(pattern)
  refined_group <- as.character(refined_group)
  finite_covariates <- rowSums(!is.finite(covariates)) == 0L
  train <- which(!holdout & is.finite(y) & finite_covariates & !is.na(pattern) & !is.na(refined_group))
  test <- which(holdout)
  test_groups <- refined_group[test]
  test_groups[is.na(pattern[test]) | !pattern[test] %in% pattern[!holdout]] <- NA_character_
  k_train <- if (is.null(kinship)) NULL else kinship[train, train, drop = FALSE]

  predict <- function(group, test_group) {
    prediction <- rep(NA_real_, length(test))
    known <- which(!is.na(test_group) & finite_covariates[test])
    if (!length(train) || !length(known)) return(prediction)
    model <- refine_model(y[train], group[train], covariates[train, , drop = FALSE], k_train)
    prediction[known] <- unname(model$means[test_group[known]])
    if (ncol(covariates) && any(is.finite(prediction))) {
      x <- sweep(covariates[test[known], , drop = FALSE], 2L, model$covariate_center)
      prediction[known] <- prediction[known] + drop(x %*% model$covariate_coefficients)
    }
    prediction
  }

  original_prediction <- predict(pattern, pattern[test])
  refined_prediction <- predict(refined_group, test_groups)
  evaluable <- is.finite(y[test])
  original_known <- evaluable & is.finite(original_prediction)
  refined_known <- evaluable & is.finite(refined_prediction)
  compared <- original_known & refined_known
  n_evaluable <- sum(evaluable)
  n_common <- sum(compared)
  coverage <- if (n_evaluable) c(sum(original_known), sum(refined_known), n_common) / n_evaluable else
    rep(NA_real_, 3L)
  original_loss <- (y[test][compared] - original_prediction[compared])^2
  refined_loss <- (y[test][compared] - refined_prediction[compared])^2
  original_rmse <- if (n_common) sqrt(mean(original_loss)) else NA_real_
  refined_rmse <- if (n_common) sqrt(mean(refined_loss)) else NA_real_
  increase <- refined_rmse - original_rmse
  interval <- c(lower = NA_real_, upper = NA_real_)
  uncertainty_method <- "unavailable"
  if (is.null(kinship) && n_common >= 10L) {
    # Paired delta-method interval conditional on training, for independent test lines.
    influence <- function(loss, rmse) if (rmse > 0) (loss - mean(loss)) / (2 * rmse) else rep(0, n_common)
    se <- stats::sd(influence(refined_loss, refined_rmse) - influence(original_loss, original_rmse)) /
      sqrt(n_common)
    interval[] <- increase + c(-1, 1) * stats::qt(0.975, n_common - 1L) * se
    uncertainty_method <- "paired_delta_independent_lines"
  }

  status <- if (!length(test)) {
    "no_holdout"
  } else if (!n_evaluable) {
    "no_evaluable_holdout"
  } else if (coverage[3L] < min_coverage) {
    "failed_coverage"
  } else if (is.null(max_rmse_increase)) {
    "no_rmse_criterion"
  } else if (anyNA(interval)) {
    "insufficient_uncertainty"
  } else if (interval[["upper"]] <= max_rmse_increase) {
    "passed"
  } else if (interval[["lower"]] > max_rmse_increase) {
    "failed_rmse"
  } else {
    "inconclusive_rmse"
  }

  list(
    status = status,
    exploratory = is.null(max_rmse_increase) || anyNA(interval) || status == "inconclusive_rmse",
    method = if (is.null(kinship)) "fixed_effects_lm" else "fixed_effects_gls",
    uncertainty_method = uncertainty_method,
    n_holdout = length(test), n_evaluable = n_evaluable, n_common = n_common, n_train = length(train),
    original_coverage = coverage[1L], refined_coverage = coverage[2L], common_coverage = coverage[3L],
    original_rmse = original_rmse, refined_rmse = refined_rmse, rmse_increase = increase,
    rmse_increase_ci = interval, min_coverage = min_coverage, max_rmse_increase = max_rmse_increase,
    predictions = data.frame(
      row_index = test, observed = y[test], pattern = pattern[test], refined_group = test_groups,
      original_prediction = original_prediction, refined_prediction = refined_prediction, compared = compared
    )
  )
}
