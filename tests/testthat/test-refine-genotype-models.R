test_that("refinement standardizes group means and contrasts over complete cases", {
  group <- rep(c("A", "B"), each = 6L)
  z <- c(1:6, 3:8)
  y <- 1 + 2 * (group == "B") + 3 * z + rep(c(1, -1, 0), 4L)
  y[1L] <- NA_real_
  model <- refine_model(y, group, cbind(z = z))

  keep <- !is.na(y)
  data <- data.frame(y = y[keep], group = group[keep], z = z[keep])
  full <- stats::lm(y ~ z + group, data)
  null <- stats::lm(y ~ z, data)
  means <- stats::predict(full, data.frame(z = mean(data$z), group = c("A", "B")))
  expected_r2 <- 1 - stats::deviance(full) / stats::deviance(null)
  expect_equal(model$counts, c(A = 5L, B = 6L))
  expect_equal(model$rows, keep)
  expect_equal(unname(model$means), unname(means))
  expect_equal(unname(model$covariate_coefficients), unname(stats::coef(full)["z"]))
  expect_equal(model$p_value, stats::anova(null, full)[2L, "Pr(>F)"])
  expect_equal(model$partial_r2, expected_r2)

  difference <- refine_difference(model, "B", "A")
  expected <- summary(full)$coefficients["groupB", ]
  expect_equal(unname(difference[1:2]), unname(expected[1:2]))
  expect_equal(unname(difference[3:4]), as.numeric(stats::confint(full, "groupB")))
  expect_equal(refine_difference(model, "B", c(B = 0.75, A = 0.25))[1:4],
               difference[1:4] / 4)
  expect_equal(refine_difference(model, "B")[["estimate"]], unname(means[2L]))
  expect_equal(refine_difference(model, "B", "B")[1:4],
               c(estimate = 0, se = 0, lower = 0, upper = 0))
})

test_that("a deterministic fit retains point predictions without intervals", {
  group <- rep(c("A", "B"), each = 4L)
  z <- rep(1:4, 2L)
  y <- 2 * z + 4 * (group == "B")
  model <- refine_model(y, group, cbind(z = z))

  expect_equal(model$status, "zero_residual_variance")
  expect_equal(model$means, c(A = 5, B = 9))
  expect_equal(unname(model$covariate_coefficients), 2)
  expect_true(all(is.na(model$covariance)))
  expect_true(is.na(refine_difference(model, "B", "A")[["lower"]]))
})

test_that("singleton and rank-deficient group designs report their limits", {
  group <- c("B", rep("A", 4L))
  y <- c(NA_real_, 1, 2, 4, 3)
  model <- refine_model(y, group, matrix(numeric(), 5L, 0L))
  expect_equal(model$counts, c(A = 4L, B = 0L))
  expect_equal(model$means, c(A = 2.5, B = NA_real_))
  expect_equal(model$df, 3L)
  expect_true(is.na(model$p_value))
  expect_true(is.na(refine_difference(model, "B", "A")[["estimate"]]))

  group <- rep(c("A", "B"), each = 4L)
  model <- refine_model(seq_along(group), group, cbind(group == "B") * 1)
  expect_equal(model$status, "rank_deficient")
  expect_equal(model$counts, c(A = 4L, B = 4L))
  expect_true(all(is.na(model$means)))
})

test_that("kinship REML uses a positive family variance in a blocked cohort", {
  family <- rep(1:8, each = 4L)
  group <- rep(c("A", "A", "B", "B"), 8L)
  y <- rep(c(-4, 2, -1, 3, 6, -6, -3, 3), each = 4L) +
    2 * (group == "B") + rep(c(-0.3, 0.3, 0.2, -0.2), 8L)
  covariates <- matrix(numeric(), length(y), 0L)
  kinship <- outer(family, family, "==") * 1
  model <- refine_model(y, group, covariates, kinship)
  linear <- refine_model(y, group, covariates)

  expect_equal(model$method, "kinship_reml")
  expect_equal(model$status, "ok")
  expect_gt(model$kinship_fraction, 0.5)
  expect_equal(model$means, c(A = 0, B = 2), tolerance = 1e-7)
  expect_lt(refine_difference(model, "B", "A")[["se"]],
            refine_difference(linear, "B", "A")[["se"]])
  expect_lt(model$p_value, 0.001)
  expect_true(is.na(model$partial_r2))

  identity <- refine_model(y, group, covariates, diag(length(y)))
  expect_equal(identity$means, linear$means)
  expect_equal(identity$covariance, linear$covariance)
  expect_equal(identity$method, "kinship_reml")

  dosage <- cbind(rV001 = as.integer(group == "B") * 2L,
                  rV002 = rep(c(0L, 2L, 0L, 2L), 8L))
  variants <- refine_variant_models(dosage, y, covariates, kinship, character(),
                                    NULL, 0.05, 6L, 2L)
  expect_equal(variants$marginal_models[["rV001"]]$covariance, model$covariance,
               ignore_attr = TRUE)
  expect_equal(variants$models[["rV002"]]$method, "kinship_reml")

  invalid <- kinship
  invalid[1L, 1L] <- -1
  expect_equal(refine_model(y, group, covariates, invalid)$status, "invalid_kinship")
})

test_that("variant screening preserves heterozygote evidence and manual protection", {
  grid <- expand.grid(a = 0:2, b = c(0L, 2L), c = c(0L, 2L), repeat_id = 1:8)
  dosage <- as.matrix(grid[, c("a", "b", "c")])
  colnames(dosage) <- c("rV001", "rV002", "rV003")
  y <- 4 * (grid$a == 1L) + rep(c(-0.4, 0.4, -0.2, 0.2), each = 12L, length.out = nrow(grid))
  covariates <- matrix(numeric(), nrow(grid), 0L)
  variants <- refine_variant_models(dosage, y, covariates, NULL, "rV003", 1, 0.05, 6L, 2L)
  rank <- variants$rank[match(colnames(dosage), variants$rank$block_id), ]

  expect_equal(variants$panel, c("rV003", "rV001"))
  expect_equal(rank$protected, c(TRUE, FALSE, TRUE))
  expect_equal(rank$evidence, c("conditional", "insufficient", "insufficient"))
  expect_lt(rank$conditional_p_adj[1L], 0.05)
  expect_equal(rank$conditional_p_adj,
               stats::p.adjust(rank$conditional_p, "holm", n = ncol(dosage)))
  expect_equal(rank$marginal_p_adj,
               stats::p.adjust(rank$marginal_p, "holm", n = ncol(dosage)))
  expect_equal(variants$models[["rV002"]]$status, "ok")
  expect_equal(variants$models[["rV002"]]$counts, c(`0` = 48L, `2` = 48L))
  expect_equal(refine_difference(variants$models[["rV001"]], "1", "0")[["estimate"]], 4)

  flipped <- refine_variant_models(dosage, -y, covariates, NULL, "rV003", 1, 0.05, 6L, 2L)
  expect_equal(flipped$rank$block_id, variants$rank$block_id)
  expect_equal(flipped$rank$protected, variants$rank$protected)
  permuted <- refine_variant_models(dosage[, 3:1], y, covariates, NULL, "rV003", 1, 0.05, 6L, 2L)
  expect_equal(permuted$rank, variants$rank)
  expect_error(refine_variant_models(dosage, y, covariates, NULL,
                                     c("rV002", "rV003"), 1, 0.05, 6L, 1L),
               "exceed")
})

test_that("variant evidence is unknown with rare or collinear states", {
  state <- rep(c(0L, 2L), each = 12L)
  dosage <- cbind(rV001 = state, rV002 = state, rV003 = c(1L, rep(0L, 23L)))
  y <- 2 * state + rep(c(-1, 1), 12L)
  covariates <- matrix(numeric(), length(y), 0L)
  variants <- refine_variant_models(dosage, y, covariates, NULL, c("rV001", "rV002"),
                                    1, 0.05, 6L, 2L)
  rank <- variants$rank[match(colnames(dosage), variants$rank$block_id), ]

  expect_true(all(is.na(rank$conditional_p_adj)))
  expect_equal(rank$protected, c(TRUE, TRUE, FALSE))
  expect_true(is.na(rank$marginal_p_adj[3L]))
  expect_equal(variants$models[["rV001"]]$status, "rank_deficient")
})

test_that("conditional fitting uses the target's complete marker background", {
  dosage <- cbind(rV001 = rep(c(0L, 1L, 2L), each = 8L),
                  rV002 = rep(c(0L, 2L), length.out = 24L))
  dosage[dosage[, "rV001"] == 1L, "rV002"] <- NA_integer_
  y <- 3 * dosage[, "rV001"] + rep(c(-0.3, 0.3, -0.4, 0.4), 6L)
  y[24L] <- NA_real_
  covariates <- matrix(numeric(), length(y), 0L)
  variants <- refine_variant_models(dosage, y, covariates, NULL, "rV002", NULL, 0.05, 6L, 2L)
  model <- variants$models[["rV001"]]

  expect_equal(model$rows, is.finite(y) & !is.na(dosage[, "rV002"]))
  expect_equal(model$counts, c(`0` = 8L, `1` = 0L, `2` = 7L))
  expect_true(is.na(model$means[["1"]]))
  expect_true(is.na(variants$rank$conditional_p_adj[match("rV001", variants$rank$block_id)]))
  expect_equal(variants$rank$block_id[variants$rank$protected], "rV002")
})

test_that("variant ranking has typed empty output when no blocks are available", {
  empty <- matrix(numeric(), 8L, 0L)
  variants <- refine_variant_models(empty, 1:8, empty, NULL, character(), NULL, 0.05, 6L, 2L)

  expect_equal(variants$rank$block_id, character())
  expect_equal(variants$rank$protected, logical())
  expect_equal(variants$rank$conditional_p_adj, numeric())
  expect_equal(nrow(variants$contrasts), 0L)
  expect_length(variants$models, 0L)
})
