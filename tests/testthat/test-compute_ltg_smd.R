test_that("compute_ltg_smd returns correct structure", {
  study_df <- make_factor_data(n_per_group = 100, k = 4, effect = 0.4,
                                seed = 20240501)
  ref_df   <- make_factor_data(n_per_group = 200, k = 4, effect = 0.4,
                                seed = 20240502)

  result <- compute_ltg_smd(
    study_data     = study_df,
    reference_data = ref_df,
    group_var      = "condition",
    score_var      = "outcome",
    items          = paste0("item", 1:4),
    group_levels   = c(reference = "a", focal = "b")
  )

  expect_s3_class(result, "ltg_smd")
  expect_named(result$estimates,
    c("hedges_g", "welch_smd", "observed_geometric",
      "external_observed", "ltg_smd"))
  expect_true(all(is.finite(result$estimates)))
  expect_named(result$factors,
    c("reliability_factor", "study_to_target_factor",
      "predicted_ratio_observed_to_LTG"))
  expect_true(result$se_analytic > 0)
  # Sanity: reliability should be reasonable with the factor-model data
  expect_true(all(result$rho_g > 0.5 & result$rho_g < 1))
})

test_that("coef_alpha matches manual formula", {
  set.seed(1)
  k <- 5; n <- 200
  items <- matrix(rnorm(k * n), ncol = k)
  manual <- (k / (k - 1)) *
    (1 - sum(apply(items, 2, var)) / var(rowSums(items)))
  expect_equal(coef_alpha(items), manual, tolerance = 1e-10)
})

test_that("LTG-SMD agrees with observed geometric SMD when reliability = 1", {
  study_df <- make_factor_data(n_per_group = 200, seed = 99)
  ref_df   <- make_factor_data(n_per_group = 500, seed = 100)

  # Override reliability to 1.0
  result <- compute_ltg_smd(
    study_data = study_df, reference_data = ref_df,
    group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 1.0, reference = 1.0),
    group_levels = c(reference = "a", focal = "b")
  )
  # Under rho = 1, LTG = external observed exactly
  expect_equal(unname(result$estimates["ltg_smd"]),
               unname(result$estimates["external_observed"]),
               tolerance = 1e-10)
})

test_that("Reliability factor matches decomposition when c_g = 1", {
  # Use the SAME data as both study and reference -> c_g should equal 1
  data_set <- make_factor_data(n_per_group = 1000, k = 3, effect = 0.5,
                                seed = 42)

  result <- compute_ltg_smd(
    study_data = data_set, reference_data = data_set,
    group_var = "condition", items = paste0("item", 1:3),
    group_levels = c(reference = "a", focal = "b")
  )

  # c_g should be 1.0 (within numerical precision)
  expect_equal(unname(result$c_g["focal"]), 1.0, tolerance = 1e-10)
  expect_equal(unname(result$c_g["reference"]), 1.0, tolerance = 1e-10)
  # Study-to-target factor should be 1.0
  expect_equal(unname(result$factors["study_to_target_factor"]), 1.0,
               tolerance = 1e-10)
})

test_that("compute_ltg_smd errors on missing group_var", {
  df <- data.frame(x = rnorm(20), grp = rep(c("a","b"), 10))
  expect_error(
    compute_ltg_smd(df, df, group_var = "nonexistent", score_var = "x",
                    rho_external = c(focal = 0.8, reference = 0.8))
  )
})

test_that("compute_ltg_smd errors on more than two groups", {
  df <- data.frame(x = rnorm(30), grp = rep(c("a","b","c"), 10))
  expect_error(
    compute_ltg_smd(df, df, group_var = "grp", score_var = "x",
                    rho_external = c(focal = 0.8, reference = 0.8))
  )
})
