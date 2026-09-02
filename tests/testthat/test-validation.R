# Validation tests
# Confirm that the package errors clearly on invalid input, and that
# reliability_estimator actually takes effect when a custom function is
# supplied. (regex patterns now properly escape literal parentheses.)

test_that("Custom reliability_estimator is actually called", {
  study_df <- make_factor_data(n_per_group = 100, k = 3, effect = 0.5,
                                seed = 11)
  ref_df <- study_df

  # Custom estimator: always returns 0.8 regardless of items
  custom_estimator <- function(items) 0.8

  result <- compute_ltg_smd(
    study_df, ref_df,
    group_var = "condition",
    items = paste0("item", 1:3),
    reliability_estimator = custom_estimator,
    group_levels = c(reference = "a", focal = "b")
  )

  # rho_g should be exactly 0.8 in both groups
  expect_equal(unname(result$rho_g["focal"]), 0.8, tolerance = 1e-12)
  expect_equal(unname(result$rho_g["reference"]), 0.8, tolerance = 1e-12)

  # Reliability factor should be (0.8 * 0.8)^{1/4} = 0.8^{1/2}
  expect_equal(unname(result$factors["reliability_factor"]),
               sqrt(0.8), tolerance = 1e-12)
})


test_that("Custom estimator with item-dependent value is propagated", {
  study_df <- make_factor_data(n_per_group = 100, k = 2, effect = 0.5,
                                seed = 22)
  ref_df <- study_df

  # Estimator depending on number of columns
  varying_estimator <- function(items) {
    if (ncol(as.matrix(items)) == 2) 0.6 else 0.9
  }

  result <- compute_ltg_smd(
    study_df, ref_df,
    group_var = "condition",
    items = c("item1", "item2"),
    reliability_estimator = varying_estimator,
    group_levels = c(reference = "a", focal = "b")
  )

  # Should be 0.6 (we supplied 2 items)
  expect_equal(unname(result$rho_g["focal"]), 0.6, tolerance = 1e-12)
})


test_that("Unknown reliability_estimator string throws an error", {
  df <- data.frame(condition = rep(c("a","b"), 10), x = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    reliability_estimator = "omega",
                    rho_external = NULL,
                    group_levels = c(reference = "a", focal = "b")),
    "reliability_estimator"
  )
})


test_that("rho_external below 0 throws an error", {
  df <- data.frame(condition = rep(c("a","b"), 10), x = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    rho_external = c(focal = -0.2, reference = 0.8),
                    group_levels = c(reference = "a", focal = "b")),
    "must lie"
  )
})


test_that("rho_external above 1 throws an error", {
  df <- data.frame(condition = rep(c("a","b"), 10), x = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    rho_external = c(focal = 1.2, reference = 0.8),
                    group_levels = c(reference = "a", focal = "b")),
    "must lie"
  )
})


test_that("rho_external with missing names throws an error", {
  df <- data.frame(condition = rep(c("a","b"), 10), x = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    rho_external = c(0.7, 0.8),  # unnamed
                    group_levels = c(reference = "a", focal = "b")),
    "named"
  )
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    rho_external = c(focal = 0.7),  # missing reference
                    group_levels = c(reference = "a", focal = "b")),
    "focal|reference"
  )
})


test_that("Negative alpha estimate triggers reliability error", {
  # Construct data with anticorrelated items, which yields negative alpha
  set.seed(99)
  n <- 100
  base_vec <- rnorm(2 * n)
  df <- data.frame(
    condition = rep(c("a","b"), each = n),
    item1 = base_vec,
    item2 = -base_vec + rnorm(2 * n, sd = 0.1)
  )
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition",
                    items = c("item1", "item2"),
                    group_levels = c(reference = "a", focal = "b")),
    "must lie|reliability"
  )
})


test_that("Zero variance in one group throws an error", {
  df <- data.frame(condition = rep(c("a","b"), each = 10))
  df$x <- c(rep(0, 10), rnorm(10))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition", score_var = "x",
                    rho_external = c(focal = 0.8, reference = 0.8),
                    group_levels = c(reference = "a", focal = "b")),
    "variance"
  )
})


test_that("Missing item column throws an informative error", {
  df <- data.frame(condition = rep(c("a","b"), each = 10),
                   item1 = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition",
                    items = c("item1", "item2_does_not_exist"),
                    group_levels = c(reference = "a", focal = "b")),
    "not found"
  )
})


test_that("Missing group_var or score_var throws informative errors", {
  df <- data.frame(condition = rep(c("a","b"), 10), x = rnorm(20))
  expect_error(
    compute_ltg_smd(df, df, group_var = "nonexistent_group",
                    score_var = "x",
                    rho_external = c(focal = 0.8, reference = 0.8)),
    "group_var"
  )
  expect_error(
    compute_ltg_smd(df, df, group_var = "condition",
                    score_var = "nonexistent_score",
                    rho_external = c(focal = 0.8, reference = 0.8)),
    "score_var"
  )
})


test_that("More than two groups throws an error", {
  df <- data.frame(grp = rep(c("a","b","c"), each = 10), x = rnorm(30))
  expect_error(
    compute_ltg_smd(df, df, group_var = "grp", score_var = "x",
                    rho_external = c(focal = 0.8, reference = 0.8))
  )
})
