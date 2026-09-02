# Bootstrap stratification tests
# Confirm that bootstrap replicates preserve the four study x group
# cell sample sizes.

make_unequal_factor_data <- function(seed = 1) {
  set.seed(seed)
  # Generate latent factor and items first, then split into two
  # different-sized groups for unequal cell sizes
  N <- 300
  latent <- rnorm(N)
  effect_vec <- c(rep(0, 30), rep(0.5, 70),    # study: 30 ctrl, 70 trt
                  rep(0, 120), rep(0.5, 80))   # ref:   120 ctrl, 80 trt
  source <- c(rep("study", 100), rep("reference", 200))
  condition <- c(rep("ctrl", 30), rep("trt", 70),
                 rep("ctrl", 120), rep("trt", 80))
  k <- 4
  items <- matrix(NA_real_, nrow = N, ncol = k)
  for (i in seq_len(k)) {
    items[, i] <- 1.0 * latent + effect_vec + rnorm(N, sd = 0.5)
  }
  colnames(items) <- paste0("item", 1:k)
  combined <- data.frame(source = source, condition = condition,
                          items, stringsAsFactors = FALSE)
  combined$outcome <- rowMeans(items)

  list(
    study = combined[combined$source == "study", , drop = FALSE],
    ref   = combined[combined$source == "reference", , drop = FALSE]
  )
}


test_that("Bootstrap stratifies on source x group_var", {
  skip_if_not_installed("boot")
  d <- make_unequal_factor_data(seed = 11)

  result <- compute_ltg_smd(
    d$study, d$ref, group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt")
  )

  ci <- ltg_smd_ci(
    result, method = "bootstrap", boot_type = "bc",
    B = 50, seed = 42,
    study_data = d$study, reference_data = d$ref,
    group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt")
  )

  expect_equal(ci$boot_details$strata_design, "source x group_var")

  # Stratum sizes correspond to the four cells
  ss <- as.numeric(ci$boot_details$stratum_sizes)
  expect_equal(sum(ss), 100 + 200)  # total
  expect_true(30 %in% ss)
  expect_true(70 %in% ss)
  expect_true(120 %in% ss)
  expect_true(80 %in% ss)
})


test_that("Bootstrap stratification produces a CI containing the point estimate", {
  skip_if_not_installed("boot")
  d <- make_unequal_factor_data(seed = 22)
  result <- compute_ltg_smd(
    d$study, d$ref, group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt")
  )
  ci <- ltg_smd_ci(
    result, method = "bootstrap", B = 100, seed = 1,
    study_data = d$study, reference_data = d$ref,
    group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt")
  )
  expect_true(
    ci$bootstrap[["lower"]] <= ci$bootstrap[["estimate"]] + 0.01 &&
    ci$bootstrap[["upper"]] >= ci$bootstrap[["estimate"]] - 0.01
  )
})


test_that("Warning when a stratum has fewer than 4 observations", {
  skip_if_not_installed("boot")
  set.seed(33)
  # Construct a setting with a tiny stratum (only 3 trt in reference)
  N_study <- 100; N_ref <- 103
  latent <- rnorm(N_study + N_ref)
  effect_vec <- c(rep(0, 50), rep(0.5, 50),  # study
                  rep(0, 100), rep(0.5, 3))   # ref: only 3 trt
  source <- c(rep("study", N_study), rep("reference", N_ref))
  condition <- c(rep("ctrl", 50), rep("trt", 50),
                 rep("ctrl", 100), rep("trt", 3))
  items <- sapply(1:3, function(i) latent + effect_vec + rnorm(N_study + N_ref, sd = 0.5))
  colnames(items) <- paste0("item", 1:3)
  combined <- data.frame(source = source, condition = condition,
                          items, stringsAsFactors = FALSE)
  combined$outcome <- rowMeans(items)

  study_df <- combined[combined$source == "study", , drop = FALSE]
  ref_df   <- combined[combined$source == "reference", , drop = FALSE]

  result <- compute_ltg_smd(
    study_df, ref_df, group_var = "condition", score_var = "outcome",
    rho_external = c(focal = 0.8, reference = 0.8),
    group_levels = c(reference = "ctrl", focal = "trt")
  )

  # Use pattern-matching to ensure we catch the right warning, ignoring
  # any incidental warnings from boot::boot itself. Note: the warning
  # message uses the plural "strata", so the pattern must too.
  expect_warning(
    ltg_smd_ci(result, method = "bootstrap", B = 30, seed = 1,
               study_data = study_df, reference_data = ref_df,
               group_var = "condition", score_var = "outcome",
               rho_external = c(focal = 0.8, reference = 0.8),
               group_levels = c(reference = "ctrl", focal = "trt")),
    regexp = "strata"
  )
})
