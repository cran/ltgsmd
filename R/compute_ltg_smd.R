#' Compute the LTG-SMD and related effect-size estimates
#'
#' Computes the latent true-score and target-population anchored geometric
#' standardized mean difference (LTG-SMD) along with conventional
#' comparators: Hedges's g, Welch-type SMD, observed geometric SMD (within
#' study), and external observed geometric SMD (in the reference sample).
#'
#' @param study_data A data frame containing the study sample, with one row
#'   per participant.
#' @param reference_data A data frame containing the reference sample. May
#'   be the same as `study_data` if an internal (holdout) reference is used;
#'   see Section 8.2 of the paper.
#' @param group_var Character. Name of the grouping variable in both data
#'   frames.
#' @param score_var Character. Name of the precomputed composite score
#'   variable. If NULL, the composite is built by averaging the items in
#'   `items`.
#' @param items Optional character vector of item names. Required if
#'   `score_var` is NULL, or if reliability is to be estimated from item-
#'   level data.
#' @param group_levels Named character vector with names "reference" (group
#'   0) and "focal" (group 1) specifying the levels of `group_var`. If NULL,
#'   the two unique values are sorted alphabetically with the first taken
#'   as reference.
#' @param reliability_estimator Either the string "alpha" (default; uses
#'   [coef_alpha()]) or a function taking a matrix/data frame of items and
#'   returning a numeric scalar in (0, 1]. Custom estimators allow plugging
#'   in `psych::omega()`, CFA-based reliability, or any user-supplied
#'   procedure.
#' @param rho_external Optional named numeric vector with names "focal" and
#'   "reference" specifying externally supplied reliability estimates. If
#'   provided, these override the internal estimator. Each value must lie
#'   in (0, 1].
#' @param se_method Character. "fourth_moment" (default; uses sample
#'   fourth moments as described in Section 4.2 and Supplementary Section C
#'   of the paper) or "normal" (uses the simpler 2*sigma^4/(m-1)
#'   approximation that is valid under normality).
#' @param verbose Logical. If TRUE, prints intermediate quantities.
#'
#' @return An object of class "ltg_smd" containing:
#' \describe{
#'   \item{point}{Numeric. The LTG-SMD point estimate.}
#'   \item{se_analytic}{Numeric. The analytic delta-method standard error
#'     (reliability-fixed approximation; see Section 4.2).}
#'   \item{se_method}{Character. Which SE method was used.}
#'   \item{estimates}{Numeric vector with elements hedges_g, welch_smd,
#'     observed_geometric, external_observed, ltg_smd.}
#'   \item{study}{Data frame of group-level study-sample statistics.}
#'   \item{reference}{Data frame of group-level reference-sample statistics.}
#'   \item{c_g}{Numeric vector of length 2: study-to-reference SD ratios
#'     (c_focal, c_reference).}
#'   \item{rho_g}{Numeric vector of length 2: reliability estimates
#'     (rho_focal, rho_reference) in the reference sample.}
#'   \item{factors}{Named numeric vector with reliability_factor,
#'     study_to_target_factor, predicted_ratio_observed_to_LTG.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The LTG-SMD estimand (Section 2.3 of the companion paper) is
#' \deqn{\delta_{\mathrm{LTG}} = \frac{\mu_{T1} - \mu_{T0}}
#'   {(\sigma^2_{T1,R}\sigma^2_{T0,R})^{1/4}}}
#' where \eqn{\mu_{Tg}} are true-score means and
#' \eqn{\sigma^2_{Tg,R}} are group-specific true-score variances in the
#' target reference population \eqn{R}. The plug-in estimator substitutes
#' observed study-sample means for the true-score means (assuming
#' \eqn{\mu_{Yg} = \mu_{Tg}}) and reference-sample reliabilities times
#' reference-sample observed variances for the true-score variances.
#'
#' Hedges's small-sample correction factor \eqn{J(N) = 1 - 3/(4N - 9)} is
#' applied to the conventional Hedges's g estimate but not to the LTG-SMD
#' or to the observed geometric SMDs, consistent with conventional practice.
#'
#' The analytic SE under `se_method = "fourth_moment"` uses
#' \eqn{\widehat{\mathrm{Var}}(s^2) = (\hat\mu_4 - \frac{m-3}{m-1}(s^2)^2)/m},
#' which is consistent under arbitrary distributions of the score with
#' finite fourth moments and reduces to the normal-theory approximation
#' \eqn{2\sigma^4/(m-1)} when the score is normal. Reliability is treated
#' as fixed at its point estimate (the reliability-fixed approximation of
#' Section 4.2); a reliability-propagated alternative is implemented by
#' the bootstrap in [ltg_smd_ci()].
#'
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
#'
#' result <- compute_ltg_smd(
#'   study_data     = study_df,
#'   reference_data = reference_df,
#'   group_var      = "condition",
#'   items          = c("item1", "item2", "item3", "item4"),
#'   group_levels   = c(reference = "control", focal = "treatment")
#' )
#' print(result)
#'
#' # Custom reliability estimator: use psych::omega instead of alpha
#' # (illustrated on a larger synthetic sample, since omega needs more
#' # data than alpha to fit stably)
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   gen_items6 <- function(n, shift = 0) {
#'     true <- rnorm(n)
#'     items <- as.data.frame(replicate(6, true + rnorm(n, sd = 0.6),
#'                                       simplify = FALSE))
#'     names(items) <- paste0("item", 1:6)
#'     items[] <- lapply(items, `+`, shift)
#'     items
#'   }
#'   study_big <- rbind(
#'     cbind(condition = "control",   gen_items6(60)),
#'     cbind(condition = "treatment", gen_items6(60, shift = 0.5))
#'   )
#'   reference_big <- rbind(
#'     cbind(condition = "control",   gen_items6(100)),
#'     cbind(condition = "treatment", gen_items6(100, shift = 0.5))
#'   )
#'   my_omega <- function(items) psych::omega(items, plot = FALSE)$omega.tot
#'   result_omega <- compute_ltg_smd(
#'     study_big, reference_big,
#'     group_var = "condition",
#'     items = paste0("item", 1:6),
#'     reliability_estimator = my_omega,
#'     group_levels = c(reference = "control", focal = "treatment")
#'   )
#'   print(result_omega)
#' }
compute_ltg_smd <- function(study_data,
                             reference_data,
                             group_var,
                             score_var = NULL,
                             items = NULL,
                             group_levels = NULL,
                             reliability_estimator = "alpha",
                             rho_external = NULL,
                             se_method = c("fourth_moment", "normal"),
                             verbose = FALSE) {

  se_method <- match.arg(se_method)

  # ---- Input validation ----
  if (!is.data.frame(study_data)) {
    stop("study_data must be a data frame.")
  }
  if (!is.data.frame(reference_data)) {
    stop("reference_data must be a data frame.")
  }
  if (!group_var %in% names(study_data)) {
    stop("group_var '", group_var, "' not found in study_data.")
  }
  if (!group_var %in% names(reference_data)) {
    stop("group_var '", group_var, "' not found in reference_data.")
  }
  check_rho_external(rho_external)

  # ---- Build composite if score_var not provided ----
  if (is.null(score_var)) {
    if (is.null(items)) {
      stop("Either score_var or items must be supplied.")
    }
    score_var <- ".ltgsmd_composite"
    study_data[[score_var]] <- build_composite(study_data, items, "mean")
    reference_data[[score_var]] <- build_composite(reference_data, items, "mean")
  } else {
    if (!score_var %in% names(study_data)) {
      stop("score_var '", score_var, "' not found in study_data.")
    }
    if (!score_var %in% names(reference_data)) {
      stop("score_var '", score_var, "' not found in reference_data.")
    }
  }

  # ---- Group-level summaries ----
  # Means come from the study sample; the reliability estimator is not
  # invoked on the study sample because the LTG-SMD denominator depends
  # only on the reference-sample variances and reliabilities.
  study_summary <- group_summary(
    study_data, group_var, score_var, items = NULL,
    reliability_estimator = reliability_estimator,
    group_levels = group_levels
  )

  # Reference sample: variances + reliability (unless rho_external supplied)
  ref_summary <- group_summary(
    reference_data, group_var, score_var,
    items = if (is.null(rho_external)) items else NULL,
    reliability_estimator = reliability_estimator,
    group_levels = group_levels
  )

  # ---- Override or validate reliabilities ----
  if (!is.null(rho_external)) {
    ref_summary$alpha[ref_summary$group == "1_focal"] <-
      rho_external[["focal"]]
    ref_summary$alpha[ref_summary$group == "0_reference"] <-
      rho_external[["reference"]]
  }

  rho_focal <- ref_summary$alpha[ref_summary$group == "1_focal"]
  rho_ref   <- ref_summary$alpha[ref_summary$group == "0_reference"]

  if (is.na(rho_focal) || is.na(rho_ref)) {
    stop("Reliability estimates are NA. Provide items + a working ",
         "reliability_estimator, or supply rho_external directly.")
  }
  check_reliability(rho_focal, "Focal-group reliability")
  check_reliability(rho_ref,   "Reference-group reliability")

  # ---- Extract group-specific quantities ----
  # Group 1 = focal, Group 0 = reference
  m_focal <- study_summary$mean[study_summary$group == "1_focal"]
  m_ref   <- study_summary$mean[study_summary$group == "0_reference"]
  s2_focal_study <- study_summary$var[study_summary$group == "1_focal"]
  s2_ref_study   <- study_summary$var[study_summary$group == "0_reference"]
  n_focal <- study_summary$n[study_summary$group == "1_focal"]
  n_ref   <- study_summary$n[study_summary$group == "0_reference"]

  s2_focal_R <- ref_summary$var[ref_summary$group == "1_focal"]
  s2_ref_R   <- ref_summary$var[ref_summary$group == "0_reference"]
  mu4_focal_R <- ref_summary$mu4[ref_summary$group == "1_focal"]
  mu4_ref_R   <- ref_summary$mu4[ref_summary$group == "0_reference"]
  m_focal_R <- ref_summary$n[ref_summary$group == "1_focal"]
  m_ref_R   <- ref_summary$n[ref_summary$group == "0_reference"]

  # ---- Effect-size estimates ----

  delta_mean <- m_focal - m_ref

  # Pooled within-study SD (Cohen's d, Hedges's g)
  s_pooled_study <- sqrt(
    ((n_focal - 1) * s2_focal_study + (n_ref - 1) * s2_ref_study) /
      (n_focal + n_ref - 2)
  )
  cohens_d <- delta_mean / s_pooled_study
  N <- n_focal + n_ref
  hedges_g <- hedges_J(N) * cohens_d

  # Welch-type SMD
  s_welch <- sqrt((s2_focal_study + s2_ref_study) / 2)
  welch_smd <- delta_mean / s_welch

  # Observed geometric SMD (within study)
  s_geom_study <- (s2_focal_study * s2_ref_study)^(1/4)
  observed_geometric <- delta_mean / s_geom_study

  # External observed geometric SMD (in reference sample)
  s_geom_R <- (s2_focal_R * s2_ref_R)^(1/4)
  external_observed <- delta_mean / s_geom_R

  # LTG-SMD: true-score geometric mean denominator
  A_focal <- rho_focal * s2_focal_R          # true-score variance, focal
  A_ref   <- rho_ref   * s2_ref_R            # true-score variance, reference
  D <- (A_focal * A_ref)^(1/4)               # geom mean of true-score SDs
  ltg_smd <- delta_mean / D

  # ---- Decomposition factors (Section 3) ----

  c_focal <- sqrt(s2_focal_study) / sqrt(s2_focal_R)
  c_ref   <- sqrt(s2_ref_study)   / sqrt(s2_ref_R)
  reliability_factor <- (rho_focal * rho_ref)^(1/4)
  study_to_target_factor <- 1 / sqrt(c_focal * c_ref)
  predicted_ratio_obs_to_LTG <- reliability_factor * study_to_target_factor

  # ---- Analytic delta-method SE (reliability-fixed) ----
  # Supplementary Section C.4, equation (C.2). Reliability is treated as fixed
  # at its point estimate.

  # Variance of group means
  var_m_focal <- s2_focal_study / n_focal
  var_m_ref   <- s2_ref_study / n_ref

  # Variance of reference-sample variances
  if (se_method == "fourth_moment") {
    var_s2_focal_R <- var_s2_plugin(mu4_focal_R, s2_focal_R, m_focal_R)
    var_s2_ref_R   <- var_s2_plugin(mu4_ref_R,   s2_ref_R,   m_ref_R)
  } else {
    var_s2_focal_R <- 2 * s2_focal_R^2 / max(m_focal_R - 1, 1)
    var_s2_ref_R   <- 2 * s2_ref_R^2   / max(m_ref_R - 1, 1)
  }
  # Variance of A_g = rho_g * s2_g, treating rho_g as fixed
  var_A_focal <- rho_focal^2 * var_s2_focal_R
  var_A_ref   <- rho_ref^2   * var_s2_ref_R

  # Gradient of delta_LTG with respect to (m_focal, m_ref, A_focal, A_ref)
  dD_dAfocal <- -delta_mean * 0.25 * A_focal^(-5/4) * A_ref^(-1/4)
  dD_dAref   <- -delta_mean * 0.25 * A_focal^(-1/4) * A_ref^(-5/4)
  dD_dmfocal <-  1 / D
  dD_dmref   <- -1 / D

  var_LTG <- dD_dmfocal^2 * var_m_focal +
             dD_dmref^2   * var_m_ref +
             dD_dAfocal^2 * var_A_focal +
             dD_dAref^2   * var_A_ref
  se_analytic <- sqrt(var_LTG)

  # ---- Pack results ----

  estimates <- c(
    hedges_g           = hedges_g,
    welch_smd          = welch_smd,
    observed_geometric = observed_geometric,
    external_observed  = external_observed,
    ltg_smd            = ltg_smd
  )

  factors <- c(
    reliability_factor              = reliability_factor,
    study_to_target_factor          = study_to_target_factor,
    predicted_ratio_observed_to_LTG = predicted_ratio_obs_to_LTG
  )

  out <- list(
    point        = ltg_smd,
    se_analytic  = se_analytic,
    se_method    = se_method,
    estimates    = estimates,
    study        = study_summary,
    reference    = ref_summary,
    c_g          = c(focal = c_focal, reference = c_ref),
    rho_g        = c(focal = rho_focal, reference = rho_ref),
    A_g          = c(focal = A_focal, reference = A_ref),
    D            = D,
    factors      = factors,
    score_var    = score_var,
    items        = items,
    reliability_estimator = reliability_estimator,
    call         = match.call()
  )
  class(out) <- "ltg_smd"

  if (verbose) print(out)

  out
}


#' @export
print.ltg_smd <- function(x, digits = 3, ...) {
  cat("LTG-SMD analysis\n")
  cat("================\n\n")
  cat("Point estimates:\n")
  print(round(x$estimates, digits))
  cat("\n")
  cat("LTG-SMD point estimate:    ", round(x$point, digits), "\n")
  cat("Analytic SE (rho-fixed, ", x$se_method, "): ",
      round(x$se_analytic, digits), "\n", sep = "")
  cat("Reference geometric true-score SD (D): ",
      round(x$D, digits), "\n\n")
  cat("Decomposition factors (Section 3):\n")
  print(round(x$factors, digits))
  cat("\n")
  cat("Study-to-reference SD ratios c_g:\n")
  print(round(x$c_g, digits))
  cat("\n")
  cat("Reference-sample reliabilities rho_g:\n")
  print(round(x$rho_g, digits))
  invisible(x)
}


#' @export
summary.ltg_smd <- function(object, ...) {
  out <- list(
    estimates   = object$estimates,
    point       = object$point,
    se_analytic = object$se_analytic,
    se_method   = object$se_method,
    factors     = object$factors,
    c_g         = object$c_g,
    rho_g       = object$rho_g,
    study       = object$study,
    reference   = object$reference,
    call        = object$call
  )
  class(out) <- "summary.ltg_smd"
  out
}


#' @export
print.summary.ltg_smd <- function(x, digits = 3, ...) {
  cat("LTG-SMD analysis summary\n")
  cat("========================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Study-sample group summaries:\n")
  print(x$study, digits = digits)
  cat("\nReference-sample group summaries:\n")
  print(x$reference, digits = digits)
  cat("\n")
  cat("Effect-size estimates:\n")
  print(round(x$estimates, digits))
  cat("\nLTG-SMD: ", round(x$point, digits),
      "  (analytic SE = ", round(x$se_analytic, digits),
      ", method = ", x$se_method, ")\n", sep = "")
  cat("\nDecomposition factors:\n")
  print(round(x$factors, digits))
  invisible(x)
}
