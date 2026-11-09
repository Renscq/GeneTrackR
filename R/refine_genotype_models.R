# Author: Rensc
# Date: 2026-11-09
# Version: dev001
# Function: Fit adjusted genotype means and evaluate variant evidence
# Input: Aligned phenotype, genotype groups, covariates, and optional GRM
# Output: Model summaries, pairwise contrasts, and a bounded variant panel

# Fit a group model and standardize its means over the fitted cohort.
refine_model <- function(y, group, covariates, kinship = NULL) {
  refine_fit_model(y, group, covariates, kinship)
}

refine_fit_model <- function(y, group, covariates, kinship, kinship_cache = NULL) {
  n_total <- length(y)
  if (length(group) != n_total || !is.matrix(covariates) || nrow(covariates) != n_total) {
    stop("Group, phenotype, and covariate inputs must have aligned rows.", call. = FALSE)
  }
  groups <- sort(unique(group[!is.na(group)]), method = "radix")
  rows <- is.finite(y) & !is.na(group) & rowSums(!is.finite(covariates)) == 0L
  group <- match(group[rows], groups)
  counts <- stats::setNames(tabulate(group, nbins = length(groups)), groups)
  observed <- which(counts > 0L)
  n <- sum(rows)
  n_cov <- ncol(covariates)
  result <- list(
    groups = groups,
    counts = counts,
    means = stats::setNames(rep(NA_real_, length(groups)), groups),
    covariance = matrix(NA_real_, length(groups), length(groups), dimnames = list(groups, groups)),
    covariate_center = rep(NA_real_, n_cov),
    covariate_coefficients = rep(NA_real_, n_cov),
    rows = rows,
    n = n,
    df = NA_integer_,
    p_value = NA_real_,
    statistic = NA_real_,
    partial_r2 = NA_real_,
    method = if (is.null(kinship)) "linear" else "kinship_reml",
    kinship_fraction = NA_real_,
    status = "insufficient_df"
  )
  df <- n - length(observed) - n_cov
  if (df <= 0L || !length(observed)) {
    return(result)
  }

  covariates <- covariates[rows, , drop = FALSE]
  center <- colMeans(covariates)
  if (n_cov) {
    covariates <- sweep(covariates, 2L, center)
    covariate_scale <- apply(abs(covariates), 2L, max)
    covariate_scale[covariate_scale == 0] <- 1
    covariates <- sweep(covariates, 2L, covariate_scale, "/")
  } else {
    covariate_scale <- numeric()
  }
  y <- y[rows]
  y_center <- mean(y)
  y <- y - y_center
  y_scale <- max(abs(y))
  if (y_scale == 0) {
    y_scale <- 1
  }
  y <- y / y_scale

  indicators <- outer(group, observed, "==") * 1
  design <- cbind(indicators, covariates)
  null_design <- cbind(1, covariates)
  fit <- stats::lm.fit(design, y)
  if (fit$rank < ncol(design)) {
    result$status <- "rank_deficient"
    return(result)
  }
  if (!is.null(kinship)) {
    spectrum <- refine_kinship_spectrum(kinship, rows, kinship_cache)
    if (is.null(spectrum)) {
      result$status <- "invalid_kinship"
      return(result)
    }
    if (sum(fit$residuals ^ 2) > n * .Machine$double.eps ^ 2) {
      profile <- refine_kinship_weights(y, design, spectrum, df)
      if (!identical(profile$status, "ok")) {
        result$status <- profile$status
        return(result)
      }
      if (!is.null(profile$weights)) {
        design <- crossprod(profile$vectors, design) * profile$weights
        null_design <- crossprod(profile$vectors, null_design) * profile$weights
        y <- drop(crossprod(profile$vectors, y)) * profile$weights
        fit <- stats::lm.fit(design, y)
      }
      result$kinship_fraction <- profile$fraction
    }
  }

  if (fit$rank < ncol(design)) {
    result$status <- "rank_deficient"
    return(result)
  }
  group_coef <- seq_along(observed)
  result$means[observed] <- fit$coefficients[group_coef] * y_scale + y_center
  result$covariate_center <- center
  result$covariate_coefficients <- fit$coefficients[length(observed) + seq_len(n_cov)] *
    y_scale / covariate_scale
  rss <- sum(fit$residuals ^ 2)
  if (rss <= n * .Machine$double.eps ^ 2) {
    result$status <- "zero_residual_variance"
    return(result)
  }
  pivot_order <- order(fit$qr$pivot)
  sigma2 <- rss / df
  coef_cov <- chol2inv(qr.R(fit$qr))[pivot_order, pivot_order, drop = FALSE] * sigma2
  result$covariance[observed, observed] <- coef_cov[group_coef, group_coef, drop = FALSE] * y_scale ^ 2
  result$df <- df
  result$status <- "ok"

  if (length(observed) > 1L) {
    # Both models use identical rows and, for GRM, the same full-model REML covariance.
    # The plug-in GLS F test and t intervals use approximate residual degrees of freedom.
    null_rss <- sum(stats::lm.fit(null_design, y)$residuals ^ 2)
    rss_change <- max(0, null_rss - rss)
    group_df <- length(observed) - 1L
    result$statistic <- rss_change / group_df / sigma2
    result$p_value <- stats::pf(result$statistic, group_df, df, lower.tail = FALSE)
    if (is.null(kinship)) {
      result$partial_r2 <- rss_change / null_rss
    }
  }
  result
}

refine_kinship_spectrum <- function(kinship, rows, cache) {
  if (!is.null(cache) && identical(rows, cache$rows)) {
    return(cache$spectrum)
  }
  if (!is.matrix(kinship) || !is.numeric(kinship) ||
      !identical(dim(kinship), rep.int(length(rows), 2L)) || any(!is.finite(kinship)) ||
      !isSymmetric(kinship, tol = 1e-8, check.attributes = FALSE)) {
    return(NULL)
  }
  kinship <- kinship[rows, rows, drop = FALSE]
  spectrum <- eigen((kinship + t(kinship)) / 2, symmetric = TRUE)
  magnitude <- max(abs(spectrum$values))
  if (min(spectrum$values) < -magnitude * 1e-8) {
    return(NULL)
  }
  if (!is.null(cache)) {
    cache$rows <- rows
    cache$spectrum <- spectrum
  }
  spectrum
}

# Profile the restricted likelihood for V proportional to (1 - h) I + h K.
refine_kinship_weights <- function(y, design, spectrum, df) {
  magnitude <- max(abs(spectrum$values))
  # A scalar identity covariance only identifies total variance.
  if (diff(range(spectrum$values)) <= magnitude * 128 * .Machine$double.eps) {
    return(list(status = "ok", fraction = NA_real_))
  }
  # Normalize the GRM to unit average diagonal before expressing the variance fraction.
  values <- pmax(spectrum$values, 0)
  values <- values / mean(values)
  rotated_design <- crossprod(spectrum$vectors, design)
  rotated_y <- drop(crossprod(spectrum$vectors, y))
  objective <- function(fraction) {
    variances <- 1 - fraction + fraction * values
    weights <- 1 / sqrt(variances)
    fit <- stats::lm.fit(rotated_design * weights, rotated_y * weights)
    rss <- sum(fit$residuals ^ 2)
    if (fit$rank < ncol(design) || !is.finite(rss) || rss <= 0) {
      return(Inf)
    }
    df * log(rss / df) + sum(log(variances)) + 2 * sum(log(abs(diag(qr.R(fit$qr)))))
  }

  singular <- min(values) <= 128 * .Machine$double.eps
  upper <- if (singular) 1 - 1e-8 else 1
  grid <- seq(0, upper, length.out = 11L)
  losses <- vapply(grid, objective, numeric(1L))
  finite <- is.finite(losses)
  if (!any(finite) || diff(range(losses[finite])) < 1e-7) {
    return(list(status = "unidentifiable_kinship"))
  }
  # Refine each local grid minimum, including intervals touching a boundary.
  minima <- which(finite & losses <= c(Inf, utils::head(losses, -1L)) &
                    losses <= c(utils::tail(losses, -1L), Inf))
  best <- which.min(losses)
  fraction <- grid[best]
  best_loss <- losses[best]
  for (i in minima) {
    interval <- grid[c(max(1L, i - 1L), min(length(grid), i + 1L))]
    candidate <- stats::optimize(objective, interval, tol = 1e-9)
    if (candidate$objective < best_loss) {
      fraction <- candidate$minimum
      best_loss <- candidate$objective
    }
  }
  if (singular && fraction >= 1 - 2e-8) {
    return(list(status = "singular_kinship_boundary"))
  }
  list(status = "ok", fraction = fraction, vectors = spectrum$vectors,
       weights = 1 / sqrt(1 - fraction + fraction * values))
}

# Return an adjusted mean, a pairwise difference, or a difference from a weighted mean.
refine_difference <- function(model, first, second = NULL, level = 0.95) {
  missing <- c(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_, df = NA_real_)
  if (length(first) != 1L || !first %in% model$groups || !is.finite(model$means[[first]]) ||
      !is.finite(model$df)) {
    return(missing)
  }
  if (is.character(second) && length(second) == 1L) {
    second <- stats::setNames(1, second)
  }
  if (!is.null(second) &&
      (!is.numeric(second) || is.null(names(second)) || anyDuplicated(names(second)) ||
       any(!names(second) %in% model$groups) || any(!is.finite(second) | second < 0) ||
       abs(sum(second) - 1) > 1e-8 || any(!is.finite(model$means[names(second[second > 0])])))) {
    return(missing)
  }
  weights <- stats::setNames(numeric(length(model$groups)), model$groups)
  if (!is.null(second)) {
    weights[names(second)] <- -second
  }
  weights[first] <- weights[first] + 1
  weights <- weights[weights != 0]
  means <- model$means[names(weights)]
  covariance <- model$covariance[names(weights), names(weights), drop = FALSE]
  estimate <- sum(weights * means)
  se <- sqrt(max(0, drop(crossprod(weights, covariance %*% weights))))
  half_width <- stats::qt((1 + level) / 2, model$df) * se
  c(estimate = estimate, se = se, lower = estimate - half_width, upper = estimate + half_width,
    df = model$df)
}

refine_variant_models <- function(dosage, y, covariates, kinship, protected,
                                  effect_margin, alpha, min_n, max_markers,
                                  candidates = character()) {
  block_ids <- sort(as.character(colnames(dosage)), method = "radix")
  if (length(setdiff(protected, block_ids))) {
    stop("Manually protected variant blocks must occur in the dosage matrix.", call. = FALSE)
  }
  protected <- intersect(block_ids, protected)
  required <- union(protected, intersect(block_ids, candidates))
  if (length(required) > max_markers) {
    stop("The specified variant blocks exceed the conditional panel budget.", call. = FALSE)
  }
  # Cache the most recent cohort eigendecomposition without retaining multiple dense GRMs.
  kinship_cache <- if (!is.null(kinship)) new.env(parent = emptyenv()) else NULL
  marginal <- lapply(block_ids, function(id) {
    refine_fit_model(y, as.character(dosage[, id]), covariates, kinship, kinship_cache)
  })
  names(marginal) <- block_ids
  omnibus <- function(model) {
    counts <- model$counts
    if (length(counts) > 1L && all(counts >= min_n)) model$p_value else NA_real_
  }
  marginal_p <- vapply(marginal, omnibus, numeric(1L))
  marginal_p_adj <- stats::p.adjust(marginal_p, method = "holm", n = length(block_ids))
  screened <- which(!is.na(marginal_p_adj) & marginal_p_adj <= alpha & !block_ids %in% required)
  screened <- screened[order(marginal_p_adj[screened], marginal_p[screened], block_ids[screened], method = "radix")]
  panel <- c(required, block_ids[utils::head(screened, max_markers - length(required))])

  conditional <- lapply(block_ids, function(id) {
    background <- setdiff(panel, id)
    if (!length(background)) {
      return(marginal[[id]])
    }
    states <- dosage[, background, drop = FALSE]
    rows <- marginal[[id]]$rows & stats::complete.cases(states)
    # Drop the reference level of each nuisance marker on the target's fitted rows.
    dummy_columns <- lapply(seq_along(background), function(i) {
      values <- states[, i]
      levels <- sort(unique(values[rows]), method = "radix")
      outer(values, levels[-1L], "==") * 1
    })
    group <- as.character(dosage[, id])
    group[!stats::complete.cases(states)] <- NA_character_
    model <- refine_fit_model(y, group, do.call(cbind, c(list(covariates), dummy_columns)),
                              kinship, kinship_cache)
    # Preserve genotype levels observed before conditioning, including zero fitted counts.
    missing_groups <- setdiff(marginal[[id]]$groups, model$groups)
    if (length(missing_groups)) {
      groups <- marginal[[id]]$groups
      model$counts <- c(model$counts, stats::setNames(rep.int(0L, length(missing_groups)), missing_groups))[groups]
      model$means <- stats::setNames(model$means[groups], groups)
      covariance <- matrix(NA_real_, length(groups), length(groups), dimnames = list(groups, groups))
      covariance[model$groups, model$groups] <- model$covariance
      model$covariance <- covariance
      model$groups <- groups
    }
    model
  })
  names(conditional) <- block_ids
  conditional_p <- vapply(conditional, omnibus, numeric(1L))
  conditional_p_adj <- stats::p.adjust(conditional_p, method = "holm", n = length(block_ids))

  # Keep all contrasts so a later merge can inspect the exact genotype pair.
  contrast_table <- function(models, type) {
    rows <- lapply(names(models), function(id) {
      model <- models[[id]]
      if (length(model$groups) < 2L) {
        return(NULL)
      }
      pairs <- utils::combn(model$groups, 2L)
      intervals <- t(vapply(seq_len(ncol(pairs)), function(i) {
        refine_difference(model, pairs[1L, i], pairs[2L, i], level = 1 - alpha)
      }, c(estimate = 0, se = 0, lower = 0, upper = 0, df = 0)))
      data.frame(block_id = id, model = type, first = pairs[1L, ], second = pairs[2L, ],
                 n_first = unname(model$counts[pairs[1L, ]]), n_second = unname(model$counts[pairs[2L, ]]),
                 intervals)
    })
    empty <- data.frame(block_id = character(), model = character(), first = character(), second = character(),
                        n_first = integer(), n_second = integer(), estimate = numeric(), se = numeric(),
                        lower = numeric(), upper = numeric(), df = numeric())
    do.call(rbind, c(list(empty), rows))
  }
  marginal_contrasts <- contrast_table(marginal, "marginal")
  conditional_contrasts <- contrast_table(conditional, "conditional")
  automatic <- character()
  if (length(effect_margin) == 1L && is.finite(effect_margin)) {
    hits <- conditional_contrasts$n_first >= min_n & conditional_contrasts$n_second >= min_n &
      (conditional_contrasts$lower > effect_margin | conditional_contrasts$upper < -effect_margin)
    supported <- block_ids[which(conditional_p_adj <= alpha)]
    automatic <- intersect(conditional_contrasts$block_id[which(hits)], supported)
  }

  score_name <- if (is.null(kinship)) "partial_r2" else "statistic"
  scores <- function(models, p_value) {
    score <- unname(vapply(models, `[[`, numeric(1L), score_name))
    score[is.na(p_value)] <- NA_real_
    score
  }
  marginal_supported <- !is.na(marginal_p_adj) & marginal_p_adj <= alpha
  conditional_supported <- !is.na(conditional_p_adj) & conditional_p_adj <= alpha
  rank <- data.frame(
    block_id = block_ids,
    protected = block_ids %in% c(protected, automatic),
    in_panel = block_ids %in% panel,
    conditional_p = unname(conditional_p),
    conditional_p_adj = unname(conditional_p_adj),
    conditional_score = scores(conditional, conditional_p),
    marginal_p = unname(marginal_p),
    marginal_p_adj = unname(marginal_p_adj),
    marginal_score = scores(marginal, marginal_p),
    n = unname(vapply(conditional, `[[`, integer(1L), "n")),
    marginal_n = unname(vapply(marginal, `[[`, integer(1L), "n")),
    evidence = rep.int("insufficient", length(block_ids))
  )
  rank$evidence[marginal_supported] <- "marginal"
  rank$evidence[conditional_supported] <- "conditional"
  # Compare incremental linear scores only in a common model family and cohort.
  rank_supported <- function(index, models, scores, p_adj, common_model) {
    if (is.null(kinship) && common_model && length(index) > 1L && all(vapply(models[index], function(model) {
      identical(model$rows, models[[index[1L]]]$rows)
    }, logical(1L)))) {
      scores <- -scores[index]
    } else {
      scores <- rep.int(0, length(index))
    }
    index[order(scores, p_adj[index], block_ids[index], method = "radix")]
  }
  conditional_index <- which(conditional_supported)
  marginal_index <- which(marginal_supported & !conditional_supported)
  unsupported_index <- which(!marginal_supported & !conditional_supported)
  index <- c(rank_supported(conditional_index, conditional, rank$conditional_score, conditional_p_adj,
                            all(block_ids[conditional_index] %in% panel)),
             rank_supported(marginal_index, marginal, rank$marginal_score, marginal_p_adj, TRUE),
             unsupported_index[order(block_ids[unsupported_index], method = "radix")])
  rank <- rank[index, , drop = FALSE]
  rownames(rank) <- NULL
  list(rank = rank, models = conditional, marginal_models = marginal,
       panel = panel, contrasts = rbind(marginal_contrasts, conditional_contrasts))
}
