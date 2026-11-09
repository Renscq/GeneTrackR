# Author: Rensc
# Date: 2026-11-09
# Version: dev003
# Function: Assemble genotype-aware refinement, clustering, and rankings
# Input: HapVariant, genotype data, and one continuous phenotype
# Output: HapRefined with reversible mappings and audit tables

refine_genotype_aware <- function(hap, phenotype, trait, sample_col, min_n,
                                  alpha, group_prefix, mixed_label, control) {
  control <- refine_genotype_controls(control, min_n, alpha)
  stop_if_not(is.character(trait) && length(trait) == 1L && !is.na(trait),
              "Genotype-aware refinement requires one numeric trait.")
  stop_if_not(is.character(group_prefix) && length(group_prefix) == 1L &&
                !is.na(group_prefix) && nzchar(group_prefix),
              "`group_prefix` must be a nonempty string.")
  stop_if_not(is.character(mixed_label) && length(mixed_label) == 1L &&
                !is.na(mixed_label) && nzchar(mixed_label) &&
                !mixed_label %in% c("0/0", "0/1", "1/1", "Unknown"),
              "`mixed_label` must be distinct from the genotype labels.")
  raw <- refine_genotype_data(hap, control$genotypes, control$protect_variants)
  ids <- rownames(raw$dosage)
  pheno_ids <- refine_input_ids(phenotype, sample_col, "Phenotype")
  stop_if_not(trait %in% names(phenotype) && trait != sample_col &&
                is.numeric(phenotype[[trait]]),
              "The selected phenotype trait must be numeric.")
  y <- as.numeric(phenotype[[trait]][match(ids, pheno_ids)])
  y[!is.finite(y)] <- NA_real_
  stop_if_not(any(!is.na(y)), "No finite phenotype values match the genotype samples.")
  holdout <- ids %in% control$validation_samples
  stop_if_not(all(control$validation_samples %in% ids),
              "Validation sample IDs must occur in the genotype cohort.")
  design <- refine_covariates(control$covariates, ids, !holdout & !is.na(y))
  kinship <- refine_align_kinship(control$kinship, ids)
  y_train <- y
  y_train[holdout] <- NA_real_
  model_rows <- !is.na(y_train) & rowSums(!is.finite(design$matrix)) == 0L
  stop_if_not(any(model_rows), "No training samples have both phenotype and complete covariates.")

  keys <- sort(unique(raw$pattern), method = "radix")
  pattern_ids <- sprintf("PHap%03d", seq_along(keys))
  pattern_index <- match(raw$pattern, keys)
  pattern <- pattern_ids[pattern_index]
  representative <- match(keys, raw$pattern)
  states <- raw$blocks[representative, , drop = FALSE]
  supported <- raw$block_supported[colnames(raw$blocks)]
  names(supported) <- colnames(raw$blocks)
  locked <- unique(raw$variant_map$block_id[raw$variant_map$protected])
  candidate_index <- match(control$candidate_variants, raw$variant_map$variant_id)
  stop_if_not(!anyNA(candidate_index), "`candidate_variants` must contain known variant IDs.")
  candidates <- unique(raw$variant_map$block_id[candidate_index])
  stop_if_not(all(candidates %in% names(supported)[supported]),
              "Candidate variants require testable biallelic genotypes.")
  evidence <- refine_variant_models(
    raw$blocks[, supported, drop = FALSE], y_train, design$matrix, kinship,
    intersect(locked, names(supported)[supported]),
    control$variant_effect_margin, alpha, min_n, control$max_markers, candidates
  )
  variant_rank <- data.table::as.data.table(evidence$rank)
  if (any(!supported)) {
    variant_rank <- data.table::rbindlist(list(
      variant_rank,
      data.table::data.table(block_id = names(supported)[!supported],
                            protected = names(supported)[!supported] %in% locked,
                            evidence = "not_testable", in_panel = FALSE)
    ), use.names = TRUE, fill = TRUE)
  }
  variant_rank[["representative_id"]] <- raw$variant_map$representative_id[
    match(variant_rank$block_id, raw$variant_map$block_id)]
  variant_rank[["member_ids"]] <- lapply(variant_rank$block_id, function(block) {
    raw$variant_map$variant_id[raw$variant_map$block_id == block]
  })
  raw$variant_map[["protected"]] <- raw$variant_map$block_id %in%
    union(locked, variant_rank$block_id[variant_rank$protected %in% TRUE])

  n_model <- tabulate(pattern_index[model_rows], nbins = length(keys))
  eligible <- which(n_model >= min_n)
  model_pattern <- pattern
  model_pattern[!pattern_index %in% eligible] <- NA_character_
  rank_rows <- model_rows & !is.na(model_pattern)
  model <- refine_model(y_train, model_pattern, design$matrix, kinship)
  protected <- unique(raw$variant_map$block_id[raw$variant_map$protected])
  blocked <- rowSums(is.na(raw$dosage[representative, , drop = FALSE]) |
                     raw$unsupported[representative, , drop = FALSE]) > 0L
  audit <- refine_merge_audit(
    states, pattern_ids, eligible, blocked, protected, evidence$models,
    model, min_n, alpha, control
  )
  cluster <- refine_merge_clusters(length(keys), audit)
  audit[["merged"]] <- cluster[audit$i] == cluster[audit$j]
  group_levels <- sort(unique(cluster))
  groups <- sprintf("%s%03d_%d", group_prefix, seq_along(group_levels),
                    tabulate(match(cluster, group_levels)))
  group_for_pattern <- groups[match(cluster, group_levels)]
  refined <- group_for_pattern[pattern_index]
  original <- unname(raw$original_hap_id)
  map <- data.table::data.table(
    hap_id = pattern_ids, pattern_key = keys, refined_hap_id = group_for_pattern,
    original_hap_ids = vapply(seq_along(keys), function(i) {
      paste(sort(unique(original[pattern_index == i]), method = "radix"), collapse = ";")
    }, character(1L))
  )
  samples <- data.table::data.table(
    sample_id = ids, hap_id = refined, original_hap_id = original, pattern_id = pattern,
    phenotype_value = y, has_phenotype = !is.na(y), used_in_model = rank_rows,
    is_validation = holdout
  )

  reference <- refine_reference(control$reference_hap, pattern_ids, original, pattern,
                                n_model, min_n)
  group_reference <- group_for_pattern[match(reference, pattern_ids)]
  model_refined <- refined
  model_refined[is.na(model_pattern)] <- NA_character_
  refined_model <- refine_model(y_train, model_refined, design$matrix, kinship)
  hap_rank <- refine_hap_ranking(refined_model, refined, y, model_rows, group_reference,
                                 control$direction, min_n, holdout)
  original_rank <- refine_hap_ranking(model, pattern, y, model_rows, reference,
                                      control$direction, min_n, holdout)
  combination_rank <- NULL
  if (control$max_order > 0L) {
    panel <- intersect(evidence$panel, names(supported)[supported])
    combination_rank <- refine_combination_rank(
      raw$blocks[, panel, drop = FALSE], y_train, design$matrix, kinship,
      control$direction, min_n, control$max_order, ids
    )
  }
  validation <- refine_holdout(
    y, model_pattern, model_refined, design$matrix, holdout, kinship,
    control$validation_min_coverage, control$validation_max_rmse_increase
  )
  validation$predictions[["sample_id"]] <- ids[validation$predictions$row_index]
  validation$merge_equivalence <- refine_validate_merges(
    audit, cluster, y, pattern, design$matrix, holdout, kinship, min_n,
    alpha, control$equivalence_margin
  )
  tables <- refine_display_tables(hap, raw, samples, map, pattern_ids,
                                  group_for_pattern, mixed_label)

  long <- data.table::data.table(
    sample_id = ids[rank_rows], hap_id = refined[rank_rows],
    trait = rep(trait, sum(rank_rows)), value = y_train[rank_rows]
  )
  summary <- data.table::data.table(
    trait = trait, hap_id = pattern_ids, sample_n = n_model,
    mean_value = unname(model$means[pattern_ids])
  )
  out <- list(
    original_hap = hap, refined_hap = tables$refined_hap,
    refined_haplotypes = tables$refined_hap$haplotypes,
    sample_refined_haplotypes = samples, haplotype_map = map,
    trait_group_map = data.table::data.table(
      trait = trait, hap_id = pattern_ids, trait_group = group_for_pattern
    ),
    phenotype_summary = summary, pairwise_test = audit, plot_data = long,
    parameters = c(list(strategy = "genotype_aware", traits = trait,
                        min_hap_samples = min_n, alpha = alpha, mixed_label = mixed_label,
                        sample_col = sample_col, group_prefix = group_prefix),
                   control[setdiff(names(control), c("genotypes", "kinship", "covariates"))]),
    pattern_hap = tables$pattern_hap, variant_map = raw$variant_map,
    sample_map = samples, genotype_dosage = raw$dosage, block_dosage = raw$blocks,
    original_haplotype_rank = original_rank, haplotype_rank = hap_rank,
    variant_rank = variant_rank, variant_contrasts = data.table::as.data.table(evidence$contrasts),
    combination_rank = combination_rank,
    merge_audit = audit, validation = validation,
    qc = list(
      n_genotyped = length(ids), n_phenotype_input = length(pheno_ids),
      n_matched = sum(ids %in% pheno_ids), n_finite_phenotype = sum(!is.na(y)),
      n_training = sum(model_rows), n_rank_model = sum(rank_rows), n_validation = sum(holdout),
      genotype_without_phenotype = setdiff(ids, pheno_ids),
      phenotype_without_genotype = setdiff(pheno_ids, ids),
      covariates_without_genotype = design$extra, genotype_without_covariates = design$missing,
      incomplete_covariates = ids[rowSums(!is.finite(design$matrix)) > 0L],
      n_original_patterns = length(keys), n_original_variants = ncol(raw$dosage),
      n_blocks = ncol(raw$blocks), n_unsupported_calls = sum(raw$unsupported),
      genotype_source = raw$source, fitted_model_status = model$status,
      model_method = model$method, kinship_fraction = model$kinship_fraction,
      ranking_status = "exploratory",
      genotype_unit = if (any(raw$unsupported) || any(raw$dosage == 1, na.rm = TRUE)) {
        "genotype_pattern"
      } else {
        "haplotype"
      }
    )
  )
  structure(out, class = c("HapRefined", "list"))
}

refine_genotype_controls <- function(control, min_n, alpha) {
  defaults <- list(
    direction = NULL, equivalence_margin = NULL, variant_effect_margin = NULL,
    genotypes = NULL, covariates = NULL, kinship = NULL,
    protect_variants = character(), candidate_variants = character(),
    reference_hap = NULL, max_diff_blocks = 1L,
    max_markers = 10L, max_order = 3L, validation_samples = character(),
    validation_min_coverage = 0.8, validation_max_rmse_increase = NULL
  )
  stop_if_not(is.list(control) && (!length(control) ||
                (!is.null(names(control)) && !anyDuplicated(names(control)) &&
                   all(names(control) %in% names(defaults)))),
              "`control` must be a named list of recognized genotype-aware settings.")
  defaults[names(control)] <- control
  is_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  is_integer_between <- function(x, lower, upper) {
    is_scalar(x) && x == floor(x) && x >= lower && x <= upper
  }
  stop_if_not(is_integer_between(min_n, 2, .Machine$integer.max),
              "`min_hap_samples` must be an integer of at least two for genotype-aware refinement.")
  stop_if_not(is_scalar(alpha) && alpha > 0 && alpha < 0.5,
              "`alpha` must be between zero and 0.5 for genotype-aware refinement.")
  stop_if_not(length(defaults$direction) == 1L && defaults$direction %in% c("higher", "lower"),
              "Set `control$direction` to 'higher' or 'lower'.")
  for (name in c("equivalence_margin", "variant_effect_margin")) {
    x <- defaults[[name]]
    stop_if_not(is.null(x) || (is_scalar(x) && x > 0),
                sprintf("`control$%s` must be a positive finite number or NULL.", name))
  }
  if (is.null(defaults$variant_effect_margin)) {
    defaults["variant_effect_margin"] <- list(defaults$equivalence_margin)
  }
  stop_if_not(is_integer_between(defaults$max_diff_blocks, 0, .Machine$integer.max),
              "`control$max_diff_blocks` must be a nonnegative integer.")
  stop_if_not(is_integer_between(defaults$max_markers, 1, 10),
              "`control$max_markers` must be an integer from 1 to 10.")
  stop_if_not(is_integer_between(defaults$max_order, 0, 3),
              "`control$max_order` must be an integer from 0 to 3.")
  for (name in c("protect_variants", "candidate_variants", "validation_samples")) {
    x <- defaults[[name]]
    stop_if_not(is.character(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x),
                sprintf("`control$%s` must contain unique nonempty IDs.", name))
  }
  ref <- defaults$reference_hap
  stop_if_not(is.null(ref) || (is.character(ref) && length(ref) == 1L &&
                !is.na(ref) && nzchar(ref)), "`control$reference_hap` must be a single ID.")
  coverage <- defaults$validation_min_coverage
  rmse <- defaults$validation_max_rmse_increase
  stop_if_not(is_scalar(coverage) && coverage > 0 && coverage <= 1,
              "`control$validation_min_coverage` must be in (0, 1].")
  stop_if_not(is.null(rmse) || (is_scalar(rmse) && rmse >= 0),
              "`control$validation_max_rmse_increase` must be nonnegative or NULL.")
  defaults
}

refine_input_ids <- function(x, sample_col, label) {
  stop_if_not(is.data.frame(x) && length(sample_col) == 1L &&
                !anyDuplicated(names(x)) && sample_col %in% names(x),
              sprintf("%s must be a table with unique columns and a sample ID column.", label))
  ids <- as.character(x[[sample_col]])
  stop_if_not(!anyNA(ids) && all(nzchar(trimws(ids))) && !anyDuplicated(ids),
              sprintf("%s sample IDs must be unique and nonempty.", label))
  ids
}

refine_covariates <- function(covariates, ids, training) {
  if (is.null(covariates)) {
    return(list(matrix = matrix(numeric(), length(ids), 0L),
                missing = character(), extra = character()))
  }
  cov_ids <- refine_input_ids(covariates, "sample_id", "Covariate")
  idx <- match(ids, cov_ids)
  vars <- setdiff(names(covariates), "sample_id")
  cols <- lapply(vars, function(name) {
    x <- covariates[[name]][idx]
    if (is.numeric(x) && is.null(dim(x))) {
      x[!is.finite(x)] <- NA_real_
      known <- unique(x[training & !is.na(x)])
      if (length(known) > 1L) return(matrix(x, ncol = 1L))
      return(matrix(ifelse(!is.na(x) & x %in% known, 0, NA_real_), ncol = 1L))
    }
    stop_if_not(is.factor(x) || is.character(x) || is.logical(x),
                sprintf("Covariate '%s' must be numeric or categorical.", name))
    x <- as.character(x)
    levels <- sort(unique(x[training & !is.na(x) & nzchar(x)]), method = "radix")
    x[!x %in% levels] <- NA_character_
    if (length(levels) < 2L) return(matrix(ifelse(is.na(x), NA_real_, 0), ncol = 1L))
    vapply(levels[-1L], function(level) as.numeric(x == level), numeric(length(ids)))
  })
  mat <- if (length(cols)) do.call(cbind, cols) else matrix(numeric(), length(ids), 0L)
  known_rows <- !is.na(idx) & rowSums(!is.finite(mat)) == 0L
  informative <- apply(mat[training & known_rows, , drop = FALSE], 2L, function(x) {
    length(unique(x)) > 1L
  })
  mat <- mat[, informative, drop = FALSE]
  if (any(!known_rows)) {
    mat <- cbind(mat, ifelse(known_rows, 0, NA_real_))
  }
  list(matrix = mat, missing = setdiff(ids, cov_ids), extra = setdiff(cov_ids, ids))
}

refine_align_kinship <- function(kinship, ids) {
  if (is.null(kinship)) return(NULL)
  stop_if_not(is.matrix(kinship) && is.numeric(kinship) && nrow(kinship) == ncol(kinship) &&
                !is.null(rownames(kinship)) && !is.null(colnames(kinship)) &&
                !anyNA(rownames(kinship)) && !anyNA(colnames(kinship)) &&
                !anyDuplicated(rownames(kinship)) && !anyDuplicated(colnames(kinship)) &&
                setequal(rownames(kinship), colnames(kinship)) && all(ids %in% rownames(kinship)) &&
                all(is.finite(kinship)),
              "Kinship must be a finite, square numeric matrix with unique IDs covering the genotype cohort.")
  kinship <- kinship[ids, ids, drop = FALSE]
  stop_if_not(isSymmetric(kinship, tol = 1e-8),
              "The aligned kinship matrix must be symmetric.")
  values <- eigen(kinship, symmetric = TRUE, only.values = TRUE)$values
  stop_if_not(min(values) >= -max(1, max(values)) * 1e-8 && max(values) > 0,
              "The kinship matrix must be positive semidefinite with positive variance.")
  kinship
}

refine_equivalence_p <- function(contrast, margin) {
  if (is.null(margin) || !is.finite(contrast[["estimate"]]) ||
      !is.finite(contrast[["se"]]) || !is.finite(contrast[["df"]]) ||
      contrast[["df"]] <= 0) return(NA_real_)
  estimate <- contrast[["estimate"]]
  se <- contrast[["se"]]
  if (se == 0) return(as.numeric(abs(estimate) >= margin))
  max(
    stats::pt((estimate + margin) / se, contrast[["df"]], lower.tail = FALSE),
    stats::pt((estimate - margin) / se, contrast[["df"]])
  )
}

refine_merge_audit <- function(states, ids, eligible, blocked, protected,
                               variant_models, hap_model, min_n, alpha, control) {
  empty <- data.table::data.table(
    i = integer(), j = integer(), group1 = character(), group2 = character(),
    n_group1 = integer(), n_group2 = integer(), n_compared = integer(),
    n_differing = integer(), differing_blocks = character(), estimate = numeric(),
    lower = numeric(), upper = numeric(), p_value = numeric(), can_merge = logical(),
    reason = character()
  )
  margin <- control$equivalence_margin
  if (length(eligible) < 2L || is.null(margin)) return(empty)
  block_ids <- colnames(states)
  rows <- lapply(seq_len(length(eligible) - 1L), function(index) {
    i <- eligible[index]
    js <- eligible[seq.int(index + 1L, length(eligible))]
    data.table::rbindlist(lapply(js, function(j) {
      compared <- !is.na(states[i, ]) & !is.na(states[j, ])
      differing <- which(compared & states[i, ] != states[j, ])
      if (length(differing) > control$max_diff_blocks) return(NULL)
      diff_ids <- block_ids[differing]
      contrast <- refine_difference(hap_model, ids[j], ids[i], level = 1 - 2 * alpha)
      p <- refine_equivalence_p(contrast, margin)
      low_effect <- function(k) {
        fit <- variant_models[[block_ids[k]]]
        if (is.null(fit)) return(FALSE)
        left <- as.character(states[i, k])
        right <- as.character(states[j, k])
        if (!isTRUE(all(fit$counts[c(left, right)] >= min_n))) return(FALSE)
        effect <- refine_difference(fit, right, left, level = 1 - 2 * alpha)
        bound <- control$variant_effect_margin
        isTRUE(effect[["lower"]] > -bound && effect[["upper"]] < bound)
      }
      reason <- if (any(diff_ids %in% protected)) {
        "protected_genotype"
      } else if (blocked[i] || blocked[j]) {
        "unresolved_genotype"
      } else if (!all(vapply(differing, low_effect, logical(1L)))) {
        "variant_effect_uncertain"
      } else if (is.na(p)) {
        "haplotype_effect_unestimable"
      } else if (p >= alpha) {
        "outside_equivalence"
      } else {
        "equivalent"
      }
      data.table::data.table(
        i = i, j = j, group1 = ids[i], group2 = ids[j],
        n_group1 = unname(hap_model$counts[ids[i]]),
        n_group2 = unname(hap_model$counts[ids[j]]),
        n_compared = sum(compared), n_differing = length(differing),
        differing_blocks = paste(diff_ids, collapse = ";"),
        estimate = contrast[["estimate"]], lower = contrast[["lower"]], upper = contrast[["upper"]],
        p_value = p, can_merge = reason == "equivalent", reason = reason
      )
    }))
  })
  data.table::rbindlist(c(list(empty), rows))
}

refine_merge_clusters <- function(n, audit) {
  cluster <- seq_len(n)
  edges <- audit[audit$can_merge, ]
  if (nrow(edges) == 0L) return(cluster)
  edges <- edges[order(edges$n_differing, edges$p_value, edges$i, edges$j), ]
  edge_keys <- paste(edges$i, edges$j, sep = ":")
  for (row in seq_len(nrow(edges))) {
    first <- cluster[edges$i[row]]
    second <- cluster[edges$j[row]]
    if (first == second) next
    left <- which(cluster == first)
    right <- which(cluster == second)
    compatible <- vapply(left, function(i) {
      all(paste(pmin(i, right), pmax(i, right), sep = ":") %in% edge_keys)
    }, logical(1L))
    if (all(compatible)) cluster[c(left, right)] <- min(first, second)
  }
  cluster
}

refine_validate_merges <- function(audit, cluster, y, pattern, covariates,
                                  holdout, kinship, min_n, alpha, margin) {
  edges <- audit[audit$can_merge & cluster[audit$i] == cluster[audit$j], ]
  if (!any(holdout) || !nrow(edges)) return(NULL)
  y[!holdout] <- NA_real_
  fit <- refine_model(y, pattern, covariates, kinship)
  n1 <- unname(fit$counts[edges$group1])
  n2 <- unname(fit$counts[edges$group2])
  results <- lapply(seq_len(nrow(edges)), function(i) {
    contrast <- refine_difference(fit, edges$group2[i], edges$group1[i], level = 1 - 2 * alpha)
    p <- if (isTRUE(n1[i] >= min_n && n2[i] >= min_n)) {
      refine_equivalence_p(contrast, margin)
    } else {
      NA_real_
    }
    data.table::data.table(
      group1 = edges$group1[i], group2 = edges$group2[i], n_group1 = n1[i], n_group2 = n2[i],
      estimate = contrast[["estimate"]], lower = contrast[["lower"]], upper = contrast[["upper"]],
      p_value = p
    )
  })
  results <- data.table::rbindlist(results)
  results[["p_adj"]] <- stats::p.adjust(ifelse(is.na(results$p_value), 1, results$p_value), "holm")
  results[["equivalent"]] <- results$p_adj < alpha
  results
}

refine_reference <- function(reference, pattern_ids, original, pattern, n, min_n) {
  if (is.null(reference)) {
    if (max(n) < min_n) return(NA_character_)
    return(pattern_ids[which.max(n)])
  }
  if (!reference %in% pattern_ids) {
    matches <- unique(pattern[original == reference])
    stop_if_not(length(matches) == 1L,
                "The reference must identify one original genotype pattern.")
    reference <- matches[[1L]]
  }
  stop_if_not(n[match(reference, pattern_ids)] >= min_n,
              "The reference pattern has too few training phenotypes.")
  reference
}

refine_hap_ranking <- function(model, group, y, model_rows, reference, direction, min_n, holdout) {
  levels <- sort(unique(group), method = "radix")
  indices <- match(group, levels)
  comparison <- t(vapply(levels, function(level) {
    refine_difference(model, level, reference)
  }, c(estimate = 0, se = 0, lower = 0, upper = 0, df = 0)))
  size <- unname(model$counts[levels])
  size[is.na(size)] <- 0L
  score <- if (direction == "higher") comparison[, "lower"] else -comparison[, "upper"]
  score[size < min_n | !is.finite(score)] <- NA_real_
  out <- data.table::data.table(
    hap_id = levels, n_carrier = tabulate(indices, length(levels)),
    n_phenotype = tabulate(indices[!is.na(y)], length(levels)),
    n_training = tabulate(indices[model_rows], length(levels)), n_model = size,
    n_validation = tabulate(indices[holdout & !is.na(y)], length(levels)),
    adjusted_mean = unname(model$means[levels]),
    reference_hap = reference, estimate = comparison[, "estimate"],
    lower = comparison[, "lower"], upper = comparison[, "upper"], rank_score = score
  )
  out <- out[order(is.na(out$rank_score), -out$rank_score, -out$n_model, out$hap_id), ]
  out[["rank"]] <- ifelse(is.na(out$rank_score), NA_integer_, seq_len(nrow(out)))
  out
}

refine_display_tables <- function(hap, raw, samples, map, pattern_ids, group_for_pattern, mixed_label) {
  labels <- matrix(c("0/0", "0/1", "1/1")[raw$dosage + 1L],
                   nrow = nrow(raw$dosage), dimnames = dimnames(raw$dosage))
  make_hap <- function(membership, levels) {
    indices <- split(seq_along(membership), factor(membership, levels = levels))
    tab <- data.table::data.table(
      hap_id = levels, sample_n = lengths(indices),
      samples = vapply(indices, function(i) paste(samples$sample_id[i], collapse = ";"), character(1L))
    )
    for (column in colnames(labels)) {
      tab[[column]] <- vapply(indices, function(i) {
        values <- unique(labels[i, column])
        if (length(values) == 1L) return(values[[1L]])
        mixed_label
      }, character(1L))
    }
    object <- hap
    object$haplotypes <- tab
    object$variants <- raw$variants
    object$sample_haplotypes <- data.table::data.table(
      sample_id = samples$sample_id, hap_id = membership,
      original_hap_id = samples$original_hap_id, hap_pattern = samples$pattern_id,
      non_missing_variant_n = rowSums(!is.na(raw$dosage))
    )
    object$meta$genotype_mode <- "dosage"
    object$meta$haplotype_n <- nrow(tab)
    object$meta$sample_n <- nrow(samples)
    object
  }
  refined_hap <- make_hap(samples$hap_id, unique(group_for_pattern))
  refined_hap$meta$refined <- TRUE
  refined_hap$meta$original_haplotype_n <- nrow(map)
  refined_hap$haplotypes[["original_hap_n"]] <- as.integer(table(factor(
    group_for_pattern, levels = refined_hap$haplotypes$hap_id
  )))
  refined_hap$haplotypes[["original_hap_ids"]] <- vapply(
    refined_hap$haplotypes$hap_id,
    function(id) paste(map$hap_id[map$refined_hap_id == id], collapse = ";"),
    character(1L)
  )
  list(refined_hap = refined_hap, pattern_hap = make_hap(samples$pattern_id, pattern_ids))
}

prepare_refined_variant_plot <- function(x, collapse_refined, variant_view,
                                         show_group_label, group_label_prefix) {
  hap <- x$original_hap
  variants <- data.table::copy(x$variant_map)
  if (identical(variant_view, "original")) {
    variants[["representative_id"]] <- variants$variant_id
    variants[["block_n"]] <- 1L
  }
  representatives <- unique(variants$representative_id)
  sites <- hap$variants[match(variants$variant_id, hap$variants$variant_id), ]
  block_members <- split(seq_len(nrow(variants)),
                         factor(variants$representative_id, levels = representatives))
  samples <- x$sample_map
  dosage <- x$genotype_dosage[samples$sample_id, variants$variant_id, drop = FALSE]

  pattern_labels <- vapply(split(samples$original_hap_id, samples$pattern_id), function(ids) {
    paste(sort(unique(ids), method = "radix"), collapse = "+")
  }, character(1L))
  split_labels <- duplicated(pattern_labels) | duplicated(pattern_labels, fromLast = TRUE)
  pattern_labels[split_labels] <- paste(pattern_labels[split_labels],
                                       names(pattern_labels)[split_labels], sep = "/")

  membership <- if (isTRUE(collapse_refined)) samples$hap_id else samples$pattern_id
  members <- split(seq_len(nrow(samples)), membership)
  ids <- names(members)
  labels <- stats::setNames(ids, ids)
  if (!isTRUE(collapse_refined)) {
    labels <- pattern_labels[ids]
    groups <- samples$hap_id[match(ids, samples$pattern_id)]
    if (identical(variant_view, "blocks") && isTRUE(show_group_label)) {
      labels[] <- paste(paste0(group_label_prefix, groups), labels, sep = " / ")
    }
    order_index <- order(match(groups, x$haplotype_rank$hap_id), labels)
  } else {
    order_index <- order(match(ids, x$haplotype_rank$hap_id), ids)
  }

  haplotypes <- data.table::data.table(
    hap_id = ids, sample_n = lengths(members),
    samples = vapply(members, function(rows) paste(samples$sample_id[rows], collapse = ";"), character(1L))
  )
  ref <- format_hap_allele(sites$ref)
  alt <- format_hap_allele(sites$alt)
  alleles <- cbind(ref, paste(ref, alt, sep = "/"), alt)
  allele_sequence <- function(row, columns) {
    tokens <- alleles[cbind(columns, dosage[row, columns] + 1L)]
    if (length(tokens) == 1L) return(tokens)
    paste(tokens, collapse = "-")
  }
  reference_sequences <- vapply(block_members[lengths(block_members) > 1L], function(columns) {
    c(REF = paste(ref[columns], collapse = "-"), ALT = paste(alt[columns], collapse = "-"))
  }, c(REF = "", ALT = ""))
  block_view <- identical(variant_view, "blocks")
  states <- if (block_view) {
    matrix(NA_character_, length(ids), length(representatives), dimnames = list(ids, representatives))
  }
  for (j in seq_along(representatives)) {
    columns <- block_members[[j]]
    # Exact redundancy gives one sequence per representative dosage.
    calls <- lapply(members, function(rows) unique(dosage[rows, representatives[j]]))
    haplotypes[[representatives[j]]] <- vapply(seq_along(members), function(i) {
      rows <- members[[i]]
      if (length(calls[[i]]) == 1L) return(allele_sequence(rows[1L], columns))
      # A varying cell names every constituent pattern and its sequence.
      rows <- rows[!duplicated(samples$pattern_id[rows])]
      member_names <- unname(pattern_labels[samples$pattern_id[rows]])
      member_alleles <- vapply(rows, allele_sequence, character(1L), columns = columns)
      member_order <- order(member_names, method = "radix")
      paste(paste0(member_names[member_order], ": ", member_alleles[member_order]), collapse = "\n")
    }, character(1L))
    if (block_view) {
      states[, j] <- vapply(calls, function(values) {
        if (length(values) > 1L) return(x$parameters$mixed_label)
        if (is.na(values)) return("Unknown")
        c("0/0", "0/1", "1/1")[values + 1L]
      }, character(1L))
    }
  }
  hap$haplotypes <- haplotypes
  hap$meta$genotype_mode <- "string"
  attr(hap, "gtr_refined_plot_layout") <- list(
    variants = variants, row_labels = labels[order_index], states = states,
    reference_sequences = reference_sequences
  )
  hap
}
