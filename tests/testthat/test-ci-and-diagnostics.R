test_that("ltg_smd_ci returns analytic interval with correct structure", {
  study_df <- make_factor_data(n_per_group = 100, k = 4, effect = 0.5,
                                seed = 2)
  ref_df   <- make_factor_data(n_per_group = 300, k = 4, effect = 0.5,
                                seed = 3)
  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "outcome",
    items = paste0("item", 1:4),
    group_levels = c(reference = "a", focal = "b")
  )
  ci <- ltg_smd_ci(result, method = "analytic")
  expect_s3_class(ci, "ltg_smd_ci")
  expect_true(!is.null(ci$analytic))
  expect_true(is.null(ci$bootstrap))
  expect_true(ci$analytic[["lower"]] <= ci$analytic[["estimate"]])
  expect_true(ci$analytic[["upper"]] >= ci$analytic[["estimate"]])
})


test_that("ltg_smd_ci bootstrap interval works with small B", {
  skip_if_not_installed("boot")
  study_df <- make_factor_data(n_per_group = 60, k = 4, effect = 0.5,
                                seed = 33)
  ref_df   <- make_factor_data(n_per_group = 100, k = 4, effect = 0.5,
                                seed = 34)
  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "outcome",
    items = paste0("item", 1:4),
    group_levels = c(reference = "a", focal = "b")
  )
  ci <- ltg_smd_ci(
    result,
    method = c("analytic", "bootstrap"),
    boot_type = "bc",
    B = 100, seed = 42,
    study_data = study_df, reference_data = ref_df,
    group_var = "condition", score_var = "outcome",
    items = paste0("item", 1:4),
    group_levels = c(reference = "a", focal = "b")
  )
  expect_true(!is.null(ci$bootstrap))
  expect_true(ci$bootstrap[["lower"]] <= ci$bootstrap[["estimate"]])
  expect_true(ci$bootstrap[["upper"]] >= ci$bootstrap[["estimate"]])
  expect_equal(ci$boot_details$type, "bc")
})


test_that("denominator_diagnostics returns 13-row data frame", {
  study_df <- make_factor_data(n_per_group = 100, k = 4, seed = 4)
  ref_df   <- make_factor_data(n_per_group = 300, k = 4, seed = 5)
  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "outcome",
    items = paste0("item", 1:4),
    group_levels = c(reference = "a", focal = "b")
  )
  diag <- denominator_diagnostics(result)
  expect_s3_class(diag, "denominator_diagnostics")
  expect_equal(nrow(diag), 13)
  expect_named(diag,
    c("Quantity", "Group_1_focal", "Group_0_reference", "Combined"))
})


test_that("denominator_sensitivity returns six denominator rows", {
  study_df <- make_factor_data(n_per_group = 100, k = 4, seed = 6)
  ref_df   <- make_factor_data(n_per_group = 300, k = 4, seed = 7)
  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "outcome",
    items = paste0("item", 1:4),
    group_levels = c(reference = "a", focal = "b")
  )
  sens <- denominator_sensitivity(result)
  expect_s3_class(sens, "denominator_sensitivity")
  expect_equal(nrow(sens), 6)
  expect_true(all(c("geometric","arithmetic","harmonic","root_mean_sq",
                    "minimum","maximum") %in% sens$denominator))
})


test_that("All six denominators agree when group variances are equal", {
  # Use rho_external with equal-variance data
  df <- make_factor_data(n_per_group = 1000, k = 4, effect = 0.5,
                          seed = 99)
  result <- compute_ltg_smd(
    df, df, group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.7, reference = 0.7),
    group_levels = c(reference = "a", focal = "b")
  )
  sens <- denominator_sensitivity(result)
  estimates <- sens$estimate
  # All six should be within ~5% when group variances ~equal
  expect_true(max(estimates) / min(estimates) < 1.05)
})
