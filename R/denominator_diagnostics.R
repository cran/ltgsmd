#' Denominator diagnostics for an LTG-SMD analysis
#'
#' Produces the diagnostic table recommended in Section 8.3 of the
#' companion paper (Table F.1 format). The table reports group-level
#' sample sizes, observed standard deviations, reference reliabilities,
#' true-score standard deviations, study-to-reference ratios c_g,
#' Hedges's g, the LTG-SMD, the reliability factor
#' \eqn{(\rho_1 \rho_0)^{1/4}}, the study-to-target factor
#' \eqn{(c_1 c_0)^{-1/2}}, the predicted ratio of the within-study
#' observed geometric SMD to the LTG-SMD (exact under the
#' decomposition), and the observed ratio of Hedges's g to the LTG-SMD
#' (which approximately matches the predicted ratio with small
#' deviations due to the pooled-versus-geometric denominator and
#' Hedges's J correction).
#'
#' @param x An "ltg_smd" object.
#' @return A data frame in a long format suitable for direct reporting.
#'   The data frame has class "denominator_diagnostics".
#' @export
#' @examples
#' set.seed(2026)
#' gen_items <- function(n, shift = 0) {
#'   true <- rnorm(n)
#'   data.frame(
#'     item1 = true + rnorm(n, sd = 0.6) + shift,
#'     item2 = true + rnorm(n, sd = 0.6) + shift,
#'     item3 = true + rnorm(n, sd = 0.6) + shift,
#'     item4 = true + rnorm(n, sd = 0.6) + shift
#'   )
#' }
#' study_df <- rbind(
#'   cbind(condition = "control",   gen_items(15)),
#'   cbind(condition = "treatment", gen_items(15, shift = 0.5))
#' )
#' reference_df <- rbind(
#'   cbind(condition = "control",   gen_items(30)),
#'   cbind(condition = "treatment", gen_items(30, shift = 0.5))
#' )
#' result <- compute_ltg_smd(study_df, reference_df,
#'   group_var = "condition", items = c("item1", "item2", "item3", "item4"),
#'   group_levels = c(reference = "control", focal = "treatment"))
#' diag <- denominator_diagnostics(result)
#' print(diag)
denominator_diagnostics <- function(x) {
  if (!inherits(x, "ltg_smd")) {
    stop("x must be an 'ltg_smd' object.")
  }

  # Extract group-level info
  s <- x$study
  r <- x$reference
  focal_idx <- which(s$group == "1_focal")
  ref_idx   <- which(s$group == "0_reference")

  n_focal <- s$n[focal_idx]
  n_ref   <- s$n[ref_idx]
  m_focal <- r$n[focal_idx]
  m_ref   <- r$n[ref_idx]
  sd_focal_study <- sqrt(s$var[focal_idx])
  sd_ref_study   <- sqrt(s$var[ref_idx])
  sd_focal_R     <- sqrt(r$var[focal_idx])
  sd_ref_R       <- sqrt(r$var[ref_idx])
  rho_focal <- x$rho_g[["focal"]]
  rho_ref   <- x$rho_g[["reference"]]
  A_focal_sd <- sqrt(x$A_g[["focal"]])
  A_ref_sd   <- sqrt(x$A_g[["reference"]])
  c_focal   <- x$c_g[["focal"]]
  c_ref     <- x$c_g[["reference"]]

  combined_n_study <- n_focal + n_ref
  combined_m_ref   <- m_focal + m_ref
  combined_sd_study <- sqrt(
    ((n_focal - 1) * s$var[focal_idx] + (n_ref - 1) * s$var[ref_idx]) /
      (n_focal + n_ref - 2)
  )
  combined_sd_R <- (r$var[focal_idx] * r$var[ref_idx])^(1/4)
  combined_rho  <- sqrt(rho_focal * rho_ref)
  combined_A_sd <- (x$A_g[["focal"]] * x$A_g[["reference"]])^(1/4)
  combined_c    <- sqrt(c_focal * c_ref)

  hedges_g <- unname(x$estimates["hedges_g"])
  ltg_smd  <- unname(x$estimates["ltg_smd"])
  rel_factor <- unname(x$factors["reliability_factor"])
  s2t_factor <- unname(x$factors["study_to_target_factor"])
  pred_ratio_obs <- unname(x$factors["predicted_ratio_observed_to_LTG"])
  observed_ratio_g_LTG <- hedges_g / ltg_smd

  diag <- data.frame(
    Quantity = c(
      "Study sample size n_g",
      "Reference sample size m_g",
      "Study observed SD s_g",
      "Reference observed SD s_{g,R}",
      "Reference reliability rho_g",
      "Reference true-score SD A_g^{1/2}",
      "Study-to-target ratio c_g",
      "Hedges's g",
      "LTG-SMD",
      "Reliability factor (rho_1 rho_0)^{1/4}",
      "Study-to-target factor (c_1 c_0)^{-1/2}",
      "Predicted ratio observed-SMD / LTG",
      "Observed ratio Hedges's g / LTG"
    ),
    Group_1_focal = c(
      sprintf("%d", n_focal),
      sprintf("%d", m_focal),
      sprintf("%.3f", sd_focal_study),
      sprintf("%.3f", sd_focal_R),
      sprintf("%.3f", rho_focal),
      sprintf("%.3f", A_focal_sd),
      sprintf("%.3f", c_focal),
      "--",
      "--",
      "--",
      "--",
      "--",
      "--"
    ),
    Group_0_reference = c(
      sprintf("%d", n_ref),
      sprintf("%d", m_ref),
      sprintf("%.3f", sd_ref_study),
      sprintf("%.3f", sd_ref_R),
      sprintf("%.3f", rho_ref),
      sprintf("%.3f", A_ref_sd),
      sprintf("%.3f", c_ref),
      "--",
      "--",
      "--",
      "--",
      "--",
      "--"
    ),
    Combined = c(
      sprintf("%d", combined_n_study),
      sprintf("%d", combined_m_ref),
      sprintf("%.3f", combined_sd_study),
      sprintf("%.3f", combined_sd_R),
      sprintf("%.3f", combined_rho),
      sprintf("%.3f", combined_A_sd),
      sprintf("%.3f", combined_c),
      sprintf("%.3f", hedges_g),
      sprintf("%.3f", ltg_smd),
      sprintf("%.3f", rel_factor),
      sprintf("%.3f", s2t_factor),
      sprintf("%.3f", pred_ratio_obs),
      sprintf("%.3f", observed_ratio_g_LTG)
    ),
    stringsAsFactors = FALSE
  )

  attr(diag, "note") <- paste(
    "Predicted ratio is exact for the within-study observed geometric SMD",
    "divided by the LTG-SMD. The observed ratio uses Hedges's g as the",
    "numerator and matches the predicted ratio approximately, with small",
    "deviations due to the pooled-vs-geometric denominator and Hedges's",
    "small-sample correction."
  )

  class(diag) <- c("denominator_diagnostics", "data.frame")
  diag
}


#' @export
print.denominator_diagnostics <- function(x, ...) {
  cat("Denominator diagnostics (Section 8.3, Table F.1 format)\n")
  cat("======================================================\n\n")
  print.data.frame(x, row.names = FALSE)
  note <- attr(x, "note")
  if (!is.null(note)) {
    cat("\nNote: ", note, "\n", sep = "")
  }
  invisible(x)
}
