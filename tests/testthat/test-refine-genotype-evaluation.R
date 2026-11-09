test_that("combination ranking uses the fixed cohort mean and keeps carrier aliases", {
  genotype <- c(rep(0L, 12L), rep(2L, 12L), 2L)
  dosage <- cbind(rV002_1 = genotype, rV001_1 = genotype)
  y <- c(rep(c(0, 1), 6L), rep(c(10, 11), 6L), NA_real_)
  ids <- sprintf("S%02d", seq_along(y))
  covariates <- matrix(numeric(), nrow = length(y), ncol = 0L)

  higher <- refine_combination_rank(dosage, y, covariates, NULL, "higher", 6L, 2L, ids)
  expect_equal(nrow(higher), 2L)
  expect_equal(higher$combination[1L], "rV001_1=1/1")
  expect_equal(higher$estimate[1L], 5)
  expect_equal(higher$score[1L], higher$lower[1L])
  expect_equal(higher$sample_n[1L], 12L)
  expect_equal(higher$phenotype_n[1L], 12L)
  expect_equal(higher$carrier_n[1L], 13L)
  expect_equal(higher$carrier_ids[[1L]], ids[13:25])
  expect_equal(higher$model_carrier_ids[[1L]], ids[13:24])
  expect_setequal(higher$aliases[[1L]], c("rV002_1=1/1", "rV001_1=1/1;rV002_1=1/1"))
  expect_equal(higher$rank, 1:2)

  lower <- refine_combination_rank(dosage[, 2:1], y, covariates, NULL, "lower", 6L, 2L, ids)
  expect_equal(lower$combination[1L], "rV001_1=0/0")
  expect_equal(lower$score[1L], -lower$upper[1L])
  perm <- rev(seq_along(y))
  expect_identical(higher, refine_combination_rank(dosage[perm, 2:1], y[perm], covariates[perm, , drop = FALSE],
                                                  NULL, "higher", 6L, 2L, ids[perm]))
})

test_that("disjoint complete cohorts receive only subset ranks", {
  genotype <- rep(c(0L, 2L), each = 6L)
  dosage <- cbind(rV001_1 = c(genotype, rep(NA_integer_, 12L)),
                  rV002_1 = c(rep(NA_integer_, 12L), genotype))
  y <- c(seq_len(12L), 100 + seq_len(12L))
  ids <- sprintf("S%02d", seq_along(y))

  result <- refine_combination_rank(dosage, y, NULL, NULL, "higher", 6L, 2L, ids)
  expect_identical(attr(result, "comparison_status"), "subset_cohorts")
  expect_equal(attr(result, "common_n"), 0L)
  expect_equal(result$comparison_n, rep(12L, 4L))
  expect_true(all(is.na(result$rank)))
  expect_false(any(result$comparable))
  expect_equal(result$subset_rank, rep(1:2, 2L))
  expect_equal(result$estimate, rep(c(3, -3), 2L))
})

test_that("observed two-locus states are compared with the best single on the same rows", {
  a <- rep(c(0L, 0L, 2L, 2L), each = 6L)
  b <- rep(c(0L, 2L, 0L, 2L), each = 6L)
  dosage <- cbind(rV001_1 = a, rV002_1 = b)
  y <- 10 * (a == b) + rep(c(-1, 0, 1), 8L)
  ids <- sprintf("S%02d", seq_along(y))

  result <- refine_combination_rank(dosage, y, NULL, NULL, "higher", 6L, 2L, ids)
  pair <- result[result$order == 2L, ]
  expect_equal(nrow(pair), 4L)
  expect_true(all(pair$single_delta_r2 > 0.95))
  expect_equal(pair$best_single, rep("rV001_1", 4L))
  expect_setequal(pair$state, c("0/0|0/0", "0/0|1/1", "1/1|0/0", "1/1|1/1"))
})

test_that("holdout error uses common rows and unseen patterns are unknown", {
  pattern <- c(rep(c("A", "B"), each = 6L), "A", "B", "A", "C")
  group <- c(rep("R", 14L), NA_character_, "R")
  y <- c(rep(c(0, 10), each = 6L), 1, 9, 50, 100)
  holdout <- seq_along(y) > 12L

  result <- refine_holdout(y, pattern, group, NULL, holdout, min_coverage = 0.5)
  expect_equal(result$n_holdout, 4L)
  expect_equal(result$n_common, 2L)
  expect_equal(result$original_coverage, 0.75)
  expect_equal(result$refined_coverage, 0.5)
  expect_equal(result$common_coverage, 0.5)
  expect_equal(result$original_rmse, 1)
  expect_equal(result$refined_rmse, 4)
  expect_equal(result$rmse_increase, 3)
  expect_true(result$exploratory)
})

test_that("paired holdout uncertainty can support a prespecified prediction limit", {
  pattern <- rep(c("A", "B"), each = 12L)
  group <- paste0("R", pattern)
  holdout <- rep(c(rep(FALSE, 6L), rep(TRUE, 6L)), 2L)
  y <- rep(c(0, 10), each = 12L) + rep(c(-1, 0, 1), 8L)
  covariates <- cbind(x = rep(c(0, 1), length.out = length(y)))

  result <- refine_holdout(y, pattern, group, covariates, holdout, max_rmse_increase = 0)
  expect_equal(result$rmse_increase, 0)
  expect_equal(unname(result$rmse_increase_ci), c(0, 0))
  expect_identical(result$status, "passed")
  expect_false(result$exploratory)
})

test_that("empty panels return a typed table and explicit comparison status", {
  result <- refine_combination_rank(matrix(nrow = 12L, ncol = 0L), seq_len(12L), NULL,
                                    NULL, "higher", 6L, 2L, as.character(seq_len(12L)))
  expect_equal(nrow(result), 0L)
  expect_identical(attr(result, "comparison_status"), "no_candidates")
  expect_type(result$carrier_ids, "list")
})

test_that("constant training outcomes still support fixed-effect point predictions", {
  y <- c(rep(10, 12L), rep(9, 10L))
  pattern <- rep("A", length(y))
  holdout <- seq_along(y) > 12L

  result <- refine_holdout(y, pattern, pattern, NULL, holdout, max_rmse_increase = 0)
  expect_equal(result$original_coverage, 1)
  expect_equal(result$original_rmse, 1)
  expect_equal(result$refined_rmse, 1)
  expect_equal(result$predictions$refined_prediction, rep(10, 10L))
})

test_that("kinship predictions do not use independent-line RMSE intervals", {
  pattern <- rep(c("A", "B"), each = 12L)
  holdout <- rep(c(rep(FALSE, 6L), rep(TRUE, 6L)), 2L)
  y <- rep(c(0, 10), each = 12L) + rep(c(-1, 0, 1), 8L)

  result <- refine_holdout(y, pattern, pattern, NULL, holdout, kinship = diag(length(y)),
                           max_rmse_increase = 0)
  expect_equal(result$common_coverage, 1)
  expect_equal(result$rmse_increase, 0)
  expect_true(all(is.na(result$rmse_increase_ci)))
  expect_true(result$exploratory)
})
test_that("discovery aliases retain distinct carrier lists outside the model cohort", {
  dosage <- cbind(
    rV001_1 = c(rep(0, 6L), rep(2, 6L), 0),
    rV002_1 = c(rep(0, 6L), rep(2, 6L), 2)
  )
  ids <- sprintf("Line%02d", seq_len(nrow(dosage)))
  result <- refine_combination_rank(
    dosage, c(0:5, 10:15, NA_real_), NULL, NULL, "higher", 6L, 1L, ids
  )
  expect_equal(nrow(result), 4L)
  alternate <- result[result$state == "1/1", ]
  expect_setequal(alternate$carrier_n, c(6L, 7L))
  expect_true("Line13" %in% alternate$carrier_ids[[match("rV002_1", alternate$subset)]])
})
