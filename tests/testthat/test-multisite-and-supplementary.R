# Tests for multisite external reference path, sensitivity_reference,
# and export_supplementary

make_multisite_factor_data <- function(seed = 1, n_per_cell = 60,
                                        n_sites = 5, k = 4,
                                        effect = 0.4) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sites)) {
    for (cond in c("ctrl", "trt")) {
      eff_val <- ifelse(cond == "trt", effect, 0)
      latent <- rnorm(n_per_cell)
      block <- data.frame(
        site = paste0("S", s),
        condition = cond,
        stringsAsFactors = FALSE
      )
      block <- block[rep(1, n_per_cell), ]
      for (i in seq_len(k)) {
        eps <- rnorm(n_per_cell, sd = 0.5)
        block[[paste0("item", i)]] <- latent + eff_val + eps
      }
      rows[[paste(s, cond, sep = "_")]] <- block
    }
  }
  do.call(rbind, rows)
}


test_that("multisite_ltg_smd works with external reference + score_var = NULL", {
  skip_if_not_installed("metafor")
  dat <- make_multisite_factor_data(seed = 7, n_per_cell = 60)
  exp_data <- dat[dat$site %in% c("S1", "S2", "S3"), ]
  ext_ref  <- dat[dat$site %in% c("S4", "S5"), ]

  result <- multisite_ltg_smd(
    data               = exp_data,
    site_var           = "site",
    group_var          = "condition",
    score_var          = NULL,
    items              = paste0("item", 1:4),
    group_levels       = c(reference = "ctrl", focal = "trt"),
    reference_strategy = "external",
    reference_data     = ext_ref,
    min_n_per_group    = 30
  )

  expect_s3_class(result, "multisite_ltg_smd")
  expect_true(nrow(result$site_table) >= 2)
  expect_true(all(is.finite(result$site_table$ltg_smd)))
  expect_equal(result$reference_strategy, "external")
})


test_that("multisite_ltg_smd se_ltg_method documentation is exposed", {
  skip_if_not_installed("metafor")
  dat <- make_multisite_factor_data(seed = 8, n_per_cell = 60)
  result <- multisite_ltg_smd(
    data = dat, site_var = "site", group_var = "condition",
    items = paste0("item", 1:4),
    group_levels = c(reference = "ctrl", focal = "trt"),
    reference_strategy = "pooled_across_sites",
    min_n_per_group = 30
  )
  expect_true(!is.null(result$se_ltg_method))
  expect_true(grepl("fixed", result$se_ltg_method))
})


test_that("sensitivity_reference compares multiple references", {
  # Primary reference matches study; alternatives have different scale
  study_df <- make_factor_data(n_per_group = 200, k = 3, effect = 0.5,
                                seed = 15, sigma_err = 0.5)
  ref_primary <- study_df

  # Alternative ref 1: narrower (halved SD via reduced loading)
  ref_alt1 <- make_factor_data(n_per_group = 200, k = 3, effect = 0.5,
                                seed = 16, loading = 0.5,
                                sigma_err = 0.25)

  # Alternative ref 2: wider (doubled SD via increased loading)
  ref_alt2 <- make_factor_data(n_per_group = 200, k = 3, effect = 0.5,
                                seed = 17, loading = 2.0,
                                sigma_err = 1.0)

  result <- compute_ltg_smd(
    study_df, ref_primary, group_var = "condition",
    items = paste0("item", 1:3),
    group_levels = c(reference = "a", focal = "b")
  )

  sens <- sensitivity_reference(
    object = result,
    alternative_references = list(narrow = ref_alt1, wide = ref_alt2),
    study_data = study_df, group_var = "condition",
    items = paste0("item", 1:3),
    group_levels = c(reference = "a", focal = "b")
  )

  expect_equal(nrow(sens), 3)
  expect_true(all(c("primary", "narrow", "wide") %in% sens$reference))
  ltg_primary <- sens$ltg_smd[sens$reference == "primary"]
  ltg_narrow  <- sens$ltg_smd[sens$reference == "narrow"]
  ltg_wide    <- sens$ltg_smd[sens$reference == "wide"]
  # Narrower reference => larger absolute LTG-SMD; wider => smaller
  expect_true(abs(ltg_narrow) > abs(ltg_primary))
  expect_true(abs(ltg_wide)   < abs(ltg_primary))
})


test_that("export_supplementary returns expected structure and optionally writes a file", {
  study_df <- make_factor_data(n_per_group = 80, k = 3, effect = 0.4,
                                seed = 25)
  ref_df <- study_df

  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition",
    items = paste0("item", 1:3),
    group_levels = c(reference = "a", focal = "b")
  )

  out <- export_supplementary(result)
  expect_true(is.list(out))
  expect_true(all(c("estimates","diagnostics","sensitivity","ltg_smd")
                  %in% names(out)))

  tf <- tempfile(fileext = ".md")
  on.exit(unlink(tf), add = TRUE)
  invisible(export_supplementary(result, file = tf))
  expect_true(file.exists(tf))
  content <- readLines(tf, warn = FALSE)
  txt <- paste(content, collapse = "\n")
  expect_true(grepl("Effect-size estimates", txt))
  expect_true(grepl("Denominator diagnostics", txt))
  expect_true(grepl("Denominator-sensitivity profile", txt))
})


test_that("export_supplementary includes CIs when supplied", {
  study_df <- make_factor_data(n_per_group = 80, k = 3, effect = 0.4,
                                seed = 35)
  ref_df <- study_df

  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition",
    items = paste0("item", 1:3),
    group_levels = c(reference = "a", focal = "b")
  )
  ci <- ltg_smd_ci(result, method = "analytic")

  tf <- tempfile(fileext = ".md")
  on.exit(unlink(tf), add = TRUE)
  invisible(export_supplementary(result, ci = ci, file = tf))
  expect_true(file.exists(tf))
  content <- readLines(tf, warn = FALSE)
  txt <- paste(content, collapse = "\n")
  expect_true(grepl("Confidence intervals", txt))
  expect_true(grepl("Analytic", txt))
})
