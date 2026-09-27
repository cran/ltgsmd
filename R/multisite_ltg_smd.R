#' Multi-site LTG-SMD analysis with cross-site reference anchoring
#'
#' Implements the multi-site workflow described in Section 7.3 and
#' illustrated in Example 3 (Section 6.3) of the paper (Many Labs 2
#' Anderson.1).
#' Each site is treated as a study sample, the cross-site pooled
#' distribution (or an external sample) is treated as the target
#' reference, and site-level Hedges's g and LTG-SMD estimates are
#' returned along with a meta-analytic comparison from [metafor::rma()].
#'
#' @param data A data frame with one row per participant containing
#'   group, site, item-level, and (optionally) composite-score columns.
#' @param site_var Character. Name of the site grouping variable.
#' @param group_var Character. Name of the experimental-condition
#'   variable.
#' @param score_var Optional character. Name of a precomputed composite
#'   score column. If NULL, built from `items` (in both `data` and, if
#'   supplied, `reference_data`).
#' @param items Optional character vector of item column names.
#' @param group_levels Optional named character vector with names
#'   "reference" and "focal".
#' @param reliability_estimator Either "alpha" (default) or a function;
#'   see [compute_ltg_smd()].
#' @param reference_strategy One of "pooled_across_sites" (default;
#'   uses all retained sites as the reference) or "external" (requires
#'   `reference_data`).
#' @param reference_data Optional data frame; required when
#'   `reference_strategy = "external"`.
#' @param min_n_per_group Integer. Minimum per-group sample size for a
#'   site to be retained. Defaults to 50.
#' @param rho_external Optional named numeric vector with names "focal"
#'   and "reference" to override item-level reliability estimation.
#'
#' @return An object of class "multisite_ltg_smd" containing:
#' \describe{
#'   \item{site_table}{Data frame of site-level estimates: site, n_focal,
#'     n_reference, hedges_g, ltg_smd, se_hedges, se_ltg.}
#'   \item{reference_summary}{List of pooled reference quantities
#'     (variances, reliabilities, true-score geometric SD).}
#'   \item{meta_hedges}{Output of [metafor::rma()] for Hedges's g.}
#'   \item{meta_ltg}{Output of [metafor::rma()] for the LTG-SMD.}
#'   \item{meta_comparison}{Data frame summarizing the two meta-analytic
#'     pools side by side.}
#'   \item{se_ltg_method}{Character. Notes that the site-level SE for
#'     the LTG-SMD treats the reference denominator as fixed.}
#' }
#'
#' @details
#' The site-level standard error for the LTG-SMD is approximated as
#' \deqn{\widehat{\mathrm{SE}}_{\mathrm{LTG, site}} =
#'   \sqrt{s^2_{1,\text{site}}/n_{1,\text{site}} +
#'         s^2_{0,\text{site}}/n_{0,\text{site}}} \,/\, D,}
#' treating the cross-site (or external) reference denominator \eqn{D}
#' as fixed. This approximation is reasonable when the reference sample
#' is large relative to each site, which is typical in coordinated
#' multi-site projects, but it understates uncertainty when the
#' reference is small. For a full reliability-propagated interval at
#' the site level, use [ltg_smd_ci()] with `method = "bootstrap"` on
#' each site separately.
#'
#' @export
#' @examples
#' set.seed(2026)
#' gen_site <- function(site_name, n_per_cond = 15, shift = 0.3) {
#'   true_ctrl <- rnorm(n_per_cond)
#'   true_trt  <- rnorm(n_per_cond)
#'   ctrl <- data.frame(
#'     site = site_name, condition = "control",
#'     item1 = true_ctrl + rnorm(n_per_cond, sd = 0.6),
#'     item2 = true_ctrl + rnorm(n_per_cond, sd = 0.6),
#'     item3 = true_ctrl + rnorm(n_per_cond, sd = 0.6),
#'     item4 = true_ctrl + rnorm(n_per_cond, sd = 0.6)
#'   )
#'   trt <- data.frame(
#'     site = site_name, condition = "treatment",
#'     item1 = true_trt + rnorm(n_per_cond, sd = 0.6) + shift,
#'     item2 = true_trt + rnorm(n_per_cond, sd = 0.6) + shift,
#'     item3 = true_trt + rnorm(n_per_cond, sd = 0.6) + shift,
#'     item4 = true_trt + rnorm(n_per_cond, sd = 0.6) + shift
#'   )
#'   rbind(ctrl, trt)
#' }
#' ml_data <- do.call(rbind, lapply(
#'   c("SiteA", "SiteB", "SiteC", "SiteD"), gen_site
#' ))
#'
#' result <- multisite_ltg_smd(
#'   data = ml_data,
#'   site_var = "site",
#'   group_var = "condition",
#'   items = paste0("item", 1:4),
#'   reference_strategy = "pooled_across_sites",
#'   min_n_per_group = 10
#' )
#' print(result)
multisite_ltg_smd <- function(data,
                               site_var,
                               group_var,
                               score_var = NULL,
                               items = NULL,
                               group_levels = NULL,
                               reliability_estimator = "alpha",
                               reference_strategy = c("pooled_across_sites",
                                                      "external"),
                               reference_data = NULL,
                               min_n_per_group = 50L,
                               rho_external = NULL) {

  reference_strategy <- match.arg(reference_strategy)
  check_rho_external(rho_external)

  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!site_var %in% names(data)) {
    stop("site_var '", site_var, "' not found in data.")
  }
  if (!group_var %in% names(data)) {
    stop("group_var '", group_var, "' not found in data.")
  }
  if (reference_strategy == "external") {
    if (is.null(reference_data)) {
      stop("reference_strategy = 'external' requires reference_data.")
    }
    if (!is.data.frame(reference_data)) {
      stop("reference_data must be a data frame.")
    }
    if (!group_var %in% names(reference_data)) {
      stop("group_var '", group_var, "' not found in reference_data.")
    }
  }

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for meta-analytic pooling.")
  }

  # ---- Build composite if not supplied ----
  # IMPORTANT: build for both data and reference_data when external.
  if (is.null(score_var)) {
    if (is.null(items)) {
      stop("Either score_var or items must be supplied.")
    }
    score_var <- ".ltgsmd_composite"
    data[[score_var]] <- build_composite(data, items, "mean")
    if (!is.null(reference_data)) {
      reference_data[[score_var]] <-
        build_composite(reference_data, items, "mean")
    }
  } else {
    if (!score_var %in% names(data)) {
      stop("score_var '", score_var, "' not found in data.")
    }
    if (!is.null(reference_data) &&
        !score_var %in% names(reference_data)) {
      stop("score_var '", score_var, "' not found in reference_data.")
    }
  }

  # ---- Determine group levels ----
  if (is.null(group_levels)) {
    levels_use <- sort(unique(data[[group_var]]))
    if (length(levels_use) != 2) {
      stop("group_var must have exactly two unique values, or supply ",
           "group_levels.")
    }
    group_levels <- c(reference = levels_use[1], focal = levels_use[2])
  }
  ref_label <- group_levels[["reference"]]
  foc_label <- group_levels[["focal"]]

  # ---- Filter sites by sample size in each condition ----
  site_n <- stats::aggregate(
    list(n = data[[score_var]]),
    by = list(site = data[[site_var]], grp = data[[group_var]]),
    FUN = function(z) sum(!is.na(z))
  )
  ref_sites <- site_n$site[site_n$grp == ref_label & site_n$n >= min_n_per_group]
  foc_sites <- site_n$site[site_n$grp == foc_label & site_n$n >= min_n_per_group]
  retained_sites <- intersect(ref_sites, foc_sites)

  if (length(retained_sites) < 2) {
    stop("Fewer than two sites pass the min_n_per_group filter. ",
         "Reduce min_n_per_group or check the data.")
  }

  # ---- Build the reference ----
  if (reference_strategy == "pooled_across_sites") {
    reference_data <- data[data[[site_var]] %in% retained_sites, , drop = FALSE]
  }

  # ---- Compute reference-level quantities once ----
  ref_summary <- group_summary(
    reference_data, group_var, score_var,
    items = if (is.null(rho_external)) items else NULL,
    reliability_estimator = reliability_estimator,
    group_levels = group_levels
  )

  if (!is.null(rho_external)) {
    ref_summary$alpha[ref_summary$group == "1_focal"] <- rho_external[["focal"]]
    ref_summary$alpha[ref_summary$group == "0_reference"] <- rho_external[["reference"]]
  }

  rho_focal <- ref_summary$alpha[ref_summary$group == "1_focal"]
  rho_ref   <- ref_summary$alpha[ref_summary$group == "0_reference"]
  check_reliability(rho_focal, "Focal-group reliability in reference")
  check_reliability(rho_ref,   "Reference-group reliability in reference")

  s2_focal_R <- ref_summary$var[ref_summary$group == "1_focal"]
  s2_ref_R   <- ref_summary$var[ref_summary$group == "0_reference"]
  A_focal <- rho_focal * s2_focal_R
  A_ref   <- rho_ref   * s2_ref_R
  D <- (A_focal * A_ref)^(1/4)

  # ---- Compute per-site effects ----
  site_rows <- list()
  for (s in retained_sites) {
    sub <- data[data[[site_var]] == s, , drop = FALSE]
    focal <- sub[sub[[group_var]] == foc_label, , drop = FALSE]
    refg  <- sub[sub[[group_var]] == ref_label, , drop = FALSE]
    n_focal <- sum(!is.na(focal[[score_var]]))
    n_ref   <- sum(!is.na(refg[[score_var]]))
    if (n_focal < 2 || n_ref < 2) next

    m_focal <- mean(focal[[score_var]], na.rm = TRUE)
    m_ref   <- mean(refg[[score_var]], na.rm = TRUE)
    s2_focal <- stats::var(focal[[score_var]], na.rm = TRUE)
    s2_ref   <- stats::var(refg[[score_var]], na.rm = TRUE)

    if (!is.finite(s2_focal) || s2_focal <= 0 ||
        !is.finite(s2_ref)   || s2_ref   <= 0) {
      next
    }

    delta_mean <- m_focal - m_ref
    s_pooled <- sqrt(
      ((n_focal - 1) * s2_focal + (n_ref - 1) * s2_ref) /
        (n_focal + n_ref - 2)
    )
    cohens_d <- delta_mean / s_pooled
    N <- n_focal + n_ref
    g <- hedges_J(N) * cohens_d
    ltg <- delta_mean / D

    # Approximate within-site SEs
    se_g <- sqrt(N / (n_focal * n_ref) + g^2 / (2 * N))
    se_ltg <- sqrt(s2_focal / n_focal + s2_ref / n_ref) / D

    site_rows[[as.character(s)]] <- data.frame(
      site = s,
      n_focal = n_focal,
      n_reference = n_ref,
      mean_diff = delta_mean,
      sd_pooled = s_pooled,
      hedges_g = g,
      ltg_smd = ltg,
      se_g = se_g,
      se_ltg = se_ltg,
      stringsAsFactors = FALSE
    )
  }

  if (length(site_rows) < 2) {
    stop("Fewer than two sites produced usable estimates. Check the ",
         "data or reduce min_n_per_group.")
  }

  site_table <- do.call(rbind, site_rows)
  rownames(site_table) <- NULL

  # ---- Meta-analysis via metafor ----
  meta_g <- metafor::rma(yi = site_table$hedges_g, sei = site_table$se_g,
                          method = "REML")
  meta_l <- metafor::rma(yi = site_table$ltg_smd, sei = site_table$se_ltg,
                          method = "REML")

  meta_comparison <- data.frame(
    metric = c("pooled_estimate", "ci_lower", "ci_upper", "se_pooled",
               "tau2", "I2_percent", "Q_stat", "Q_df", "Q_pvalue", "k_sites"),
    hedges_g = c(meta_g$beta, meta_g$ci.lb, meta_g$ci.ub, meta_g$se,
                 meta_g$tau2, meta_g$I2, meta_g$QE, meta_g$k - 1,
                 meta_g$QEp, meta_g$k),
    ltg_smd = c(meta_l$beta, meta_l$ci.lb, meta_l$ci.ub, meta_l$se,
                meta_l$tau2, meta_l$I2, meta_l$QE, meta_l$k - 1,
                meta_l$QEp, meta_l$k),
    stringsAsFactors = FALSE
  )

  ref_info <- list(
    rho_focal = rho_focal,
    rho_reference = rho_ref,
    var_focal_R = s2_focal_R,
    var_reference_R = s2_ref_R,
    true_score_var_focal = A_focal,
    true_score_var_reference = A_ref,
    true_score_geometric_SD = D,
    n_retained_sites = length(retained_sites),
    n_total_participants = sum(site_table$n_focal + site_table$n_reference)
  )

  out <- list(
    site_table        = site_table,
    reference_summary = ref_info,
    meta_hedges       = meta_g,
    meta_ltg          = meta_l,
    meta_comparison   = meta_comparison,
    retained_sites    = retained_sites,
    min_n_per_group   = min_n_per_group,
    reference_strategy = reference_strategy,
    se_ltg_method     = paste(
      "Site-level SE for LTG-SMD treats the cross-site (or external)",
      "reference denominator D as fixed (study-only contribution).",
      "Use ltg_smd_ci(method = 'bootstrap') for each site separately to",
      "propagate uncertainty in the reference denominator."
    ),
    call              = match.call()
  )
  class(out) <- "multisite_ltg_smd"
  out
}


#' @export
print.multisite_ltg_smd <- function(x, digits = 3, ...) {
  cat("Multi-site LTG-SMD analysis\n")
  cat("===========================\n\n")
  cat("Sites retained: ", x$reference_summary$n_retained_sites,
      " (min_n_per_group = ", x$min_n_per_group, ")\n", sep = "")
  cat("Total participants: ", x$reference_summary$n_total_participants, "\n", sep = "")
  cat("Reference strategy: ", x$reference_strategy, "\n", sep = "")
  cat("\nReference true-score geometric SD (D): ",
      round(x$reference_summary$true_score_geometric_SD, digits), "\n", sep = "")
  cat("Reference reliabilities: focal = ",
      round(x$reference_summary$rho_focal, digits),
      ", reference = ",
      round(x$reference_summary$rho_reference, digits), "\n\n", sep = "")

  cat("Meta-analytic comparison:\n")
  print.data.frame(
    data.frame(
      metric = x$meta_comparison$metric,
      hedges_g = round(x$meta_comparison$hedges_g, digits),
      ltg_smd  = round(x$meta_comparison$ltg_smd, digits),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )

  cat("\nFirst rows of site-level estimates:\n")
  st <- utils::head(x$site_table, 8)
  st_print <- data.frame(
    site = st$site,
    n_focal = st$n_focal,
    n_reference = st$n_reference,
    hedges_g = round(st$hedges_g, digits),
    ltg_smd = round(st$ltg_smd, digits),
    stringsAsFactors = FALSE
  )
  print.data.frame(st_print, row.names = FALSE)

  if (nrow(x$site_table) > 8) {
    cat("(", nrow(x$site_table) - 8, " more rows in site_table)\n", sep = "")
  }

  cat("\nNote on SE: ", x$se_ltg_method, "\n", sep = "")

  invisible(x)
}
