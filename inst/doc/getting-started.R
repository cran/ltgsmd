## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)
library(ltgsmd)
set.seed(20240501)

## ----single-study-example, eval = FALSE---------------------------------------
# # Simulated example
# n <- 60; k <- 4
# study_df <- data.frame(condition = rep(c("control", "treatment"), each = n))
# for (i in 1:k) {
#   eff <- ifelse(study_df$condition == "treatment", 0.5, 0)
#   study_df[[paste0("item", i)]] <- rnorm(2 * n, mean = eff)
# }
# study_df$outcome <- rowMeans(study_df[, paste0("item", 1:k)])
# 
# # Holdout reference
# ref_df <- data.frame(condition = rep(c("control", "treatment"), each = n))
# for (i in 1:k) {
#   eff <- ifelse(ref_df$condition == "treatment", 0.5, 0)
#   ref_df[[paste0("item", i)]] <- rnorm(2 * n, mean = eff)
# }
# ref_df$outcome <- rowMeans(ref_df[, paste0("item", 1:k)])
# 
# result <- compute_ltg_smd(
#   study_data     = study_df,
#   reference_data = ref_df,
#   group_var      = "condition",
#   score_var      = "outcome",
#   items          = paste0("item", 1:k),
#   group_levels   = c(reference = "control", focal = "treatment")
# )
# 
# print(result)

## ----single-study-ci, eval = FALSE--------------------------------------------
# ci <- ltg_smd_ci(
#   result,
#   method         = c("analytic", "bootstrap"),
#   boot_type      = "bc",
#   B              = 2000,
#   seed           = 20240501,
#   study_data     = study_df,
#   reference_data = ref_df,
#   group_var      = "condition",
#   score_var      = "outcome",
#   items          = paste0("item", 1:k),
#   group_levels   = c(reference = "control", focal = "treatment")
# )
# print(ci)

## ----single-study-diag, eval = FALSE------------------------------------------
# diag <- denominator_diagnostics(result)
# print(diag)

## ----single-study-sens, eval = FALSE------------------------------------------
# sens <- denominator_sensitivity(result)
# print(sens)

## ----external-ref-example, eval = FALSE---------------------------------------
# # Suppose 'study_df' is a focused US 18-30 subsample and 'norm_df' is
# # the rest of an Open Psychometrics dataset serving as the broader
# # reference distribution.
# 
# result_external <- compute_ltg_smd(
#   study_data     = study_df,
#   reference_data = norm_df,
#   group_var      = "gender",
#   score_var      = "neuroticism_score",
#   items          = c("N1", "N2", "N3", "N4"),
#   group_levels   = c(reference = "Male", focal = "Female")
# )
# print(result_external)

## ----multisite-example, eval = FALSE------------------------------------------
# ml2_result <- multisite_ltg_smd(
#   data               = ml2_data,
#   site_var           = "Source.Global",
#   group_var          = "condition",
#   score_var          = "SWB",
#   items              = paste0("and_item_", 1:25),
#   reference_strategy = "pooled_across_sites",
#   min_n_per_group    = 50
# )
# print(ml2_result)

## ----ref-sensitivity, eval = FALSE--------------------------------------------
# # Compare three plausible references for an analysis
# sens_ref <- sensitivity_reference(
#   object                  = result,
#   alternative_references  = list(
#     norm_alt1 = ref_df_alt1,
#     norm_alt2 = ref_df_alt2
#   ),
#   study_data   = study_df,
#   group_var    = "condition",
#   score_var    = "outcome",
#   items        = paste0("item", 1:k),
#   group_levels = c(reference = "control", focal = "treatment")
# )
# print(sens_ref)

## ----supplementary, eval = FALSE----------------------------------------------
# export_supplementary(result, ci = ci, file = "supplementary_appendix.md")

