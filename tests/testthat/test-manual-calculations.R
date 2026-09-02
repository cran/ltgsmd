# Hand-computed verification tests
# These tests build small datasets where every quantity can be checked
# against an independent manual computation. If any of these fail, the
# core math is wrong.

test_that("compute_ltg_smd matches hand-computed values on a tiny example", {
  # Build a deterministic study sample with known means and variances
  study_df <- data.frame(
    condition = c(rep("ctrl", 4), rep("trt", 4)),
    outcome = c(1, 2, 3, 4,    # control: mean = 2.5, var = 5/3
                3, 4, 5, 6)    # treatment: mean = 4.5, var = 5/3
  )

  # Reference sample (different distribution)
  ref_df <- data.frame(
    condition = c(rep("ctrl", 6), rep("trt", 6)),
    outcome = c(0, 1, 2, 3, 4, 5,        # ctrl: mean = 2.5, var = 3.5
                0, 1, 2, 3, 4, 5)        # trt:  mean = 2.5, var = 3.5
  )

  # Use fixed reliability so the comparison is deterministic
  result <- compute_ltg_smd(
    study_df, ref_df,
    group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt"),
    se_method = "normal"
  )

  # Hand computation:
  # study mean diff: 4.5 - 2.5 = 2.0
  # ref var (both groups): 3.5
  # A_g = 0.8 * 3.5 = 2.8
  # D = (2.8 * 2.8)^(1/4) = 2.8^(1/2) = sqrt(2.8) = 1.67332...
  # LTG-SMD = 2.0 / 1.67332 = 1.19523...
  expect_equal(unname(result$estimates["ltg_smd"]),
               2.0 / sqrt(2.8), tolerance = 1e-10)

  # Hand computation of Hedges's g:
  # var_study (each group): 5/3
  # pooled SD = sqrt(5/3) (since both vars equal)
  # Cohen's d = 2.0 / sqrt(5/3) = 2.0 * sqrt(3/5) ~ 1.54919
  # J(8) = 1 - 3/(4*8-9) = 1 - 3/23 = 20/23 ~ 0.86957
  # Hedges's g ~ 1.54919 * 0.86957 = 1.34712
  expect_equal(unname(result$estimates["hedges_g"]),
               (1 - 3/23) * 2.0 / sqrt(5/3), tolerance = 1e-10)

  # External observed (in ref sample, both vars 3.5):
  # 2.0 / (3.5 * 3.5)^(1/4) = 2.0 / sqrt(3.5)
  expect_equal(unname(result$estimates["external_observed"]),
               2.0 / sqrt(3.5), tolerance = 1e-10)
})


test_that("Reliability factor matches (rho_1 rho_0)^{1/4} exactly", {
  set.seed(7)
  n <- 100
  df <- data.frame(condition = rep(c("a","b"), each = n))
  df$x <- rnorm(2 * n) + ifelse(df$condition == "b", 0.3, 0)

  result <- compute_ltg_smd(
    df, df, group_var = "condition", score_var = "x",
    rho_external = c(focal = 0.7, reference = 0.6),
    group_levels = c(reference = "a", focal = "b")
  )

  # Hand computation
  expected_rel_factor <- (0.7 * 0.6)^(1/4)
  expect_equal(unname(result$factors["reliability_factor"]),
               expected_rel_factor, tolerance = 1e-12)

  # When study_data == reference_data, c_g = 1 exactly, so
  # study_to_target_factor = 1, and predicted ratio = reliability factor
  expect_equal(unname(result$factors["study_to_target_factor"]),
               1.0, tolerance = 1e-12)
  expect_equal(unname(result$factors["predicted_ratio_observed_to_LTG"]),
               expected_rel_factor, tolerance = 1e-12)
})


test_that("Study-to-target factor matches 1/sqrt(c_focal * c_ref)", {
  # Build a study sample with study variance HALF the reference variance
  set.seed(1)
  n_study <- 200; n_ref <- 200
  # Reference: var = 4
  ref_df <- data.frame(condition = rep(c("a","b"), each = n_ref))
  ref_df$x <- 2 * rnorm(2 * n_ref) + ifelse(ref_df$condition == "b", 0.5, 0)
  # Study: var = 1 (sd = 1, so half of ref sd = 2)
  study_df <- data.frame(condition = rep(c("a","b"), each = n_study))
  study_df$x <- rnorm(2 * n_study) + ifelse(study_df$condition == "b", 0.5, 0)

  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "x",
    rho_external = c(focal = 0.9, reference = 0.9),
    group_levels = c(reference = "a", focal = "b")
  )

  # c_g should be approximately 0.5 (study sd is half of reference sd)
  # study_to_target_factor = 1 / sqrt(c_focal * c_ref)
  expected <- 1 / sqrt(result$c_g[["focal"]] * result$c_g[["reference"]])
  expect_equal(unname(result$factors["study_to_target_factor"]),
               unname(expected), tolerance = 1e-12)
  # Roughly equals 2 because c_g ~ 0.5
  expect_true(unname(result$factors["study_to_target_factor"]) > 1.5 &&
              unname(result$factors["study_to_target_factor"]) < 2.5)
})
