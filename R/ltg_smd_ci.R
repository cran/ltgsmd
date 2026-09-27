#' Compute confidence intervals for the LTG-SMD
#'
#' Computes one or more confidence intervals for the LTG-SMD, including the
#' analytic delta-method interval (reliability-fixed; Section 4.2 of the
#' paper) and the bias-corrected (BC) or bias-corrected-and-accelerated
#' (BCa) nonparametric bootstrap interval.
#'
#' @param object An object of class "ltg_smd" returned by
#'   [compute_ltg_smd()].
#' @param method Character vector with elements from c("analytic",
#'   "bootstrap"). Both can be requested.
#' @param boot_type Character. Type of bootstrap interval: "bc"
#'   (bias-corrected, default), "bca" (bias-corrected and accelerated),
#'   or "percentile" (uncorrected percentile interval).
#' @param B Integer. Number of bootstrap replicates. Default 2000.
#' @param level Numeric. Confidence level, default 0.95.
#' @param seed Optional integer for reproducibility of the bootstrap.
#' @param parallel Character; passed to [boot::boot()]. Default "no".
#' @param ncpus Integer; passed to [boot::boot()] if parallel != "no".
#' @param study_data Data frame. Required for bootstrap (the original study
#'   sample passed to compute_ltg_smd).
#' @param reference_data Data frame. Required for bootstrap.
#' @param group_var,score_var,items,group_levels,reliability_estimator
#'   Same as in [compute_ltg_smd()]; required when bootstrap is requested.
#' @param rho_external Optional; see [compute_ltg_smd()].
#'
#' @return An object of class "ltg_smd_ci" containing:
#' \describe{
#'   \item{analytic}{Named numeric vector with elements estimate, lower,
#'     upper, se, level.}
#'   \item{bootstrap}{Named numeric vector with bootstrap interval, or
#'     NULL if bootstrap not requested.}
#'   \item{boot_details}{If bootstrap requested, list with B, type, seed,
#'     z0, acceleration, strata_design, and the bootstrap distribution.}
#' }
#'
#' @details
#' **Stratification.** The bootstrap resamples within the four strata
#' defined by the cross-classification of source (study vs reference) and
#' experimental group (focal vs reference). This preserves the sample
#' sizes of all four study x group cells across replicates, which is the
#' appropriate stratification for a two-group LTG-SMD analysis. Earlier
#' versions of this package stratified only on source; that design
#' allowed the group-1/group-0 ratio within each source to vary across
#' replicates, producing unstable confidence intervals under small
#' samples or unequal group sizes.
#'
#' **Reliability propagation.** Because the bootstrap recomputes the
#' entire LTG-SMD plug-in within each replicate, item-level reliability
#' is re-estimated each time (when `items` is supplied). This propagates
#' uncertainty in the reliability estimator through to the interval. By
#' contrast, the analytic interval treats reliability as fixed at its
#' point estimate.
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
#' result <- compute_ltg_smd(study_df, reference_df,
#'   group_var = "condition", items = c("item1", "item2", "item3", "item4"),
#'   group_levels = c(reference = "control", focal = "treatment"))
#'
#' # B is kept small here so the example runs quickly; in practice use a
#' # larger B (e.g. 2000, the default) for stable bootstrap intervals.
#' ci <- ltg_smd_ci(result,
#'   method = c("analytic", "bootstrap"),
#'   B = 200, seed = 20260814,
#'   study_data = study_df, reference_data = reference_df,
#'   group_var = "condition", items = c("item1", "item2", "item3", "item4"),
#'   group_levels = c(reference = "control", focal = "treatment"))
#' print(ci)
ltg_smd_ci <- function(object,
                        method = c("analytic", "bootstrap"),
                        boot_type = c("bc", "bca", "percentile"),
                        B = 2000,
                        level = 0.95,
                        seed = NULL,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = 1L,
                        study_data = NULL,
                        reference_data = NULL,
                        group_var = NULL,
                        score_var = NULL,
                        items = NULL,
                        group_levels = NULL,
                        reliability_estimator = "alpha",
                        rho_external = NULL) {

  if (!inherits(object, "ltg_smd")) {
    stop("object must be an 'ltg_smd' object from compute_ltg_smd().")
  }

  method <- match.arg(method, choices = c("analytic", "bootstrap"),
                     several.ok = TRUE)
  boot_type <- match.arg(boot_type)
  parallel <- match.arg(parallel)
  alpha_lev <- 1 - level
  z <- stats::qnorm(1 - alpha_lev / 2)

  check_rho_external(rho_external)

  analytic_ci <- NULL
  if ("analytic" %in% method) {
    analytic_ci <- c(
      estimate = unname(object$point),
      lower    = unname(object$point - z * object$se_analytic),
      upper    = unname(object$point + z * object$se_analytic),
      se       = unname(object$se_analytic),
      level    = level
    )
  }

  bootstrap_ci <- NULL
  boot_details <- NULL
  if ("bootstrap" %in% method) {
    if (is.null(study_data) || is.null(reference_data) || is.null(group_var)) {
      stop("Bootstrap requires study_data, reference_data, and group_var ",
           "(plus score_var or items) to be supplied.")
    }
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package 'boot' is required for bootstrap intervals. ",
           "Install with install.packages('boot').")
    }

    if (!is.null(seed)) set.seed(seed)

    # Construct combined data with two tagging columns: source (study /
    # reference) and the group variable. Bootstrap resampling stratifies
    # on the interaction of these two columns, which preserves all four
    # source x group cell sample sizes across replicates.
    study_data$.source <- "study"
    reference_data$.source <- "reference"
    combined <- rbind(study_data, reference_data)
    combined$.stratum <- interaction(
      combined$.source, combined[[group_var]], drop = TRUE
    )
    strata <- factor(combined$.stratum)

    n_stratum <- table(strata)
    if (any(n_stratum < 4)) {
      warning("One or more source x group strata has fewer than 4 ",
              "observations; bootstrap intervals may be unstable.")
    }

    boot_stat <- function(d, i) {
      sub <- d[i, , drop = FALSE]
      study_sub <- sub[sub$.source == "study", , drop = FALSE]
      ref_sub   <- sub[sub$.source == "reference", , drop = FALSE]
      if (nrow(study_sub) < 4 || nrow(ref_sub) < 4) return(NA_real_)
      # Drop the .source and .stratum columns before passing
      study_sub$.source <- NULL
      study_sub$.stratum <- NULL
      ref_sub$.source <- NULL
      ref_sub$.stratum <- NULL
      out <- tryCatch(
        compute_ltg_smd(
          study_data     = study_sub,
          reference_data = ref_sub,
          group_var      = group_var,
          score_var      = score_var,
          items          = items,
          group_levels   = group_levels,
          reliability_estimator = reliability_estimator,
          rho_external   = rho_external,
          verbose        = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(out)) return(NA_real_)
      out$point
    }

    b_out <- boot::boot(
      data       = combined,
      statistic  = boot_stat,
      R          = B,
      strata     = strata,
      parallel   = parallel,
      ncpus      = ncpus
    )

    t_orig <- as.numeric(b_out$t0)
    t_boot <- as.numeric(b_out$t)
    t_boot <- t_boot[is.finite(t_boot)]
    B_eff <- length(t_boot)

    if (B_eff < 10) {
      stop("Fewer than 10 bootstrap replicates returned finite values. ",
           "Check sample sizes and reliability estimator.")
    }

    # BC factor
    p_below <- mean(t_boot < t_orig)
    if (p_below %in% c(0, 1)) {
      z0 <- 0
    } else {
      z0 <- stats::qnorm(p_below)
    }

    # Acceleration via jackknife (only if BCa)
    accel <- NA_real_
    if (boot_type == "bca") {
      n_total <- nrow(combined)
      jk <- vapply(seq_len(n_total), function(idx) {
        sub <- combined[-idx, , drop = FALSE]
        study_sub <- sub[sub$.source == "study", , drop = FALSE]
        ref_sub   <- sub[sub$.source == "reference", , drop = FALSE]
        if (nrow(study_sub) < 4 || nrow(ref_sub) < 4) return(NA_real_)
        study_sub$.source <- NULL
        study_sub$.stratum <- NULL
        ref_sub$.source <- NULL
        ref_sub$.stratum <- NULL
        out <- tryCatch(
          compute_ltg_smd(
            study_data     = study_sub,
            reference_data = ref_sub,
            group_var      = group_var,
            score_var      = score_var,
            items          = items,
            group_levels   = group_levels,
            reliability_estimator = reliability_estimator,
            rho_external   = rho_external,
            verbose        = FALSE
          ),
          error = function(e) NULL
        )
        if (is.null(out)) NA_real_ else out$point
      }, numeric(1))
      jk <- jk[is.finite(jk)]
      jk_mean <- mean(jk)
      num <- sum((jk_mean - jk)^3)
      den <- 6 * (sum((jk_mean - jk)^2))^(3/2)
      accel <- if (den > 0) num / den else 0
    }

    # Compute percentile cutoffs depending on type
    if (boot_type == "percentile") {
      alpha_lo <- alpha_lev / 2
      alpha_hi <- 1 - alpha_lev / 2
    } else if (boot_type == "bc") {
      alpha_lo <- stats::pnorm(2 * z0 + stats::qnorm(alpha_lev / 2))
      alpha_hi <- stats::pnorm(2 * z0 + stats::qnorm(1 - alpha_lev / 2))
    } else { # bca
      a <- accel
      za <- stats::qnorm(alpha_lev / 2)
      zb <- stats::qnorm(1 - alpha_lev / 2)
      alpha_lo <- stats::pnorm(z0 + (z0 + za) / (1 - a * (z0 + za)))
      alpha_hi <- stats::pnorm(z0 + (z0 + zb) / (1 - a * (z0 + zb)))
    }

    boot_lower <- stats::quantile(t_boot, alpha_lo, type = 7)
    boot_upper <- stats::quantile(t_boot, alpha_hi, type = 7)

    bootstrap_ci <- c(
      estimate = t_orig,
      lower    = unname(boot_lower),
      upper    = unname(boot_upper),
      se       = stats::sd(t_boot),
      level    = level
    )

    boot_details <- list(
      B               = B,
      B_effective     = B_eff,
      type            = boot_type,
      seed            = seed,
      z0              = z0,
      acceleration    = accel,
      strata_design   = "source x group_var",
      stratum_sizes   = n_stratum,
      boot_dist       = t_boot
    )
  }

  out <- list(
    analytic     = analytic_ci,
    bootstrap    = bootstrap_ci,
    boot_details = boot_details,
    method       = method,
    level        = level
  )
  class(out) <- "ltg_smd_ci"
  out
}


#' @export
print.ltg_smd_ci <- function(x, digits = 3, ...) {
  cat("LTG-SMD confidence intervals (level = ", x$level, ")\n", sep = "")
  cat("================================================\n\n")
  if (!is.null(x$analytic)) {
    cat("Analytic delta-method (rho-fixed):\n")
    cat("  Estimate: ", round(x$analytic[["estimate"]], digits), "\n", sep = "")
    cat("  SE:       ", round(x$analytic[["se"]], digits), "\n", sep = "")
    cat("  ",
        format(100 * x$level), "% CI: [",
        round(x$analytic[["lower"]], digits), ", ",
        round(x$analytic[["upper"]], digits), "]\n\n", sep = "")
  }
  if (!is.null(x$bootstrap)) {
    cat("Nonparametric bootstrap (", x$boot_details$type, ", B = ",
        x$boot_details$B, ", stratified on ",
        x$boot_details$strata_design, "):\n", sep = "")
    cat("  Estimate: ", round(x$bootstrap[["estimate"]], digits), "\n", sep = "")
    cat("  SE:       ", round(x$bootstrap[["se"]], digits), "\n", sep = "")
    cat("  ",
        format(100 * x$level), "% CI: [",
        round(x$bootstrap[["lower"]], digits), ", ",
        round(x$bootstrap[["upper"]], digits), "]\n", sep = "")
    if (!is.null(x$boot_details$z0)) {
      cat("  Bias-correction factor z0: ",
          round(x$boot_details$z0, digits), "\n", sep = "")
    }
    if (!is.na(x$boot_details$acceleration)) {
      cat("  Acceleration: ",
          round(x$boot_details$acceleration, digits), "\n", sep = "")
    }
  }
  invisible(x)
}
