#' Compute coefficient alpha (Cronbach's alpha)
#'
#' Computes Cronbach's coefficient alpha from item-level data. Uses
#' listwise-complete cases by default. The implementation follows the
#' standard formula
#' \deqn{\alpha = \frac{k}{k-1}\left(1 - \frac{\sum_i s_i^2}{s_T^2}\right)}
#' where \eqn{k} is the number of items, \eqn{s_i^2} are item variances,
#' and \eqn{s_T^2} is the variance of the total score.
#'
#' @param items A data frame, matrix, or list-like object whose columns are
#'   item-level responses (one row per respondent).
#' @param na.rm Logical. If TRUE (default), use complete cases only.
#' @return A numeric scalar, the value of coefficient alpha.
#' @export
#' @examples
#' set.seed(1)
#' items <- matrix(rnorm(400), ncol = 4)
#' coef_alpha(items)
coef_alpha <- function(items, na.rm = TRUE) {
  if (is.list(items) && !is.data.frame(items)) {
    items <- do.call(cbind, items)
  }
  items <- as.matrix(items)

  if (na.rm) {
    items <- items[stats::complete.cases(items), , drop = FALSE]
  }

  k <- ncol(items)
  if (k < 2) {
    stop("coef_alpha requires at least two items.")
  }

  item_vars <- apply(items, 2, stats::var)
  total_var <- stats::var(rowSums(items))

  if (total_var <= 0) {
    return(NA_real_)
  }

  alpha <- (k / (k - 1)) * (1 - sum(item_vars) / total_var)
  as.numeric(alpha)
}


#' Compute the fourth central moment of a numeric vector
#'
#' Helper for the analytic delta-method variance of \eqn{s^2} in the
#' LTG-SMD asymptotic variance derivation (see Section 4 of the paper).
#'
#' @param x A numeric vector.
#' @param na.rm Logical. If TRUE (default), missing values are removed.
#' @return Numeric scalar, the sample fourth central moment.
#' @noRd
fourth_moment <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  m <- mean(x)
  mean((x - m)^4)
}


#' Variance of s^2 using the sample fourth-moment plug-in
#'
#' Implements the variance approximation
#' \deqn{\widehat{\mathrm{Var}}(s^2) = \frac{\hat\mu_4 -
#'   \frac{m-3}{m-1}(s^2)^2}{m}}
#' which reduces to \eqn{2 \sigma^4 / (m-1)} under normality. This is the
#' plug-in used in Section 6 of the companion paper.
#'
#' @param mu4 Sample fourth central moment.
#' @param s2 Sample variance.
#' @param m Sample size.
#' @return Estimated variance of s^2.
#' @noRd
var_s2_plugin <- function(mu4, s2, m) {
  if (m < 4) return(NA_real_)
  val <- (mu4 - ((m - 3) / (m - 1)) * s2^2) / m
  # The estimator can occasionally be negative for small m under heavy
  # skew; fall back to the normal-approximation lower bound in that case.
  normal_approx <- 2 * s2^2 / (m - 1)
  if (!is.finite(val) || val <= 0) {
    return(normal_approx)
  }
  val
}


#' Compute the small-sample correction factor for Hedges's g
#'
#' Computes the J correction \eqn{J(N) = 1 - 3 / (4N - 9)} used to convert
#' Cohen's d into Hedges's g, where N is the total sample size of the two
#' groups (Hedges, 1981).
#'
#' @param N Total sample size (sum of group sizes).
#' @return Numeric scalar, the J correction factor.
#' @noRd
hedges_J <- function(N) {
  if (N <= 2) return(NA_real_)
  1 - 3 / (4 * N - 9)
}


#' Build a composite score from item-level data
#'
#' Computes the row-mean (or row-sum) of a set of items. Used internally
#' by [compute_ltg_smd()] when a precomputed composite score is not
#' supplied.
#'
#' @param data A data frame containing the items as columns.
#' @param items Character vector of item column names.
#' @param method Either "mean" (default) or "sum".
#' @param min_valid Minimum number of valid items required per respondent.
#'   Respondents with fewer valid items receive NA. Defaults to all items.
#' @return Numeric vector of composite scores.
#' @noRd
build_composite <- function(data, items, method = c("mean", "sum"),
                            min_valid = NULL) {
  method <- match.arg(method)
  missing_cols <- setdiff(items, names(data))
  if (length(missing_cols) > 0) {
    stop("Item columns not found in data: ",
         paste(missing_cols, collapse = ", "))
  }
  M <- as.matrix(data[, items, drop = FALSE])

  if (is.null(min_valid)) min_valid <- length(items)
  n_valid <- rowSums(!is.na(M))

  out <- rep(NA_real_, nrow(M))
  ok <- n_valid >= min_valid

  if (method == "mean") {
    out[ok] <- rowMeans(M[ok, , drop = FALSE], na.rm = TRUE)
  } else {
    out[ok] <- rowSums(M[ok, , drop = FALSE], na.rm = TRUE)
  }
  out
}


#' Validate a reliability estimate
#'
#' Internal helper that enforces \eqn{0 < \rho \le 1}. Stops with an
#' informative error if the value is missing, non-finite, non-positive,
#' or exceeds 1.
#'
#' @param rho Numeric scalar.
#' @param name Character. Used in the error message to identify the source.
#' @return Invisibly, the validated value.
#' @noRd
check_reliability <- function(rho, name = "reliability estimate") {
  if (length(rho) != 1L || is.na(rho) || !is.finite(rho)) {
    stop(sprintf("%s is missing or non-finite (got: %s).",
                 name, format(rho)))
  }
  if (rho <= 0 || rho > 1) {
    stop(sprintf(paste0(
      "%s must lie in (0, 1] (got: %s). Negative, zero, or above-1 ",
      "coefficient alpha is symptomatic of items that do not consistently ",
      "measure a common construct, or of strong sampling-induced item ",
      "covariance distortions. The LTG-SMD framework requires a defensible ",
      "reliability estimator with a value in (0, 1]; consider using ",
      "omega, a CFA-based reliability, or supplying rho_external from ",
      "a separate calibration sample."),
      name, format(rho, digits = 4)))
  }
  invisible(rho)
}


#' Validate an external reliability specification
#'
#' Checks that `rho_external` is a named numeric vector with the
#' required names "focal" and "reference", and that each value lies in
#' (0, 1]. Returns invisibly.
#'
#' @param rho_external Named numeric vector or NULL.
#' @noRd
check_rho_external <- function(rho_external) {
  if (is.null(rho_external)) return(invisible(NULL))
  if (!is.numeric(rho_external) || is.null(names(rho_external))) {
    stop("rho_external must be a named numeric vector.")
  }
  if (!all(c("focal", "reference") %in% names(rho_external))) {
    stop("rho_external must contain names 'focal' and 'reference'.")
  }
  check_reliability(rho_external[["focal"]],
                    "rho_external['focal']")
  check_reliability(rho_external[["reference"]],
                    "rho_external['reference']")
  invisible(rho_external)
}


#' Resolve a reliability_estimator argument to a function
#'
#' Internal helper: converts the user-facing reliability_estimator
#' argument (either "alpha" or a function) into a callable function that
#' takes a matrix of items and returns a numeric scalar.
#'
#' @param reliability_estimator Character "alpha" or a function.
#' @return A function.
#' @noRd
resolve_reliability_estimator <- function(reliability_estimator) {
  if (is.function(reliability_estimator)) {
    return(reliability_estimator)
  }
  if (is.character(reliability_estimator) &&
      length(reliability_estimator) == 1L) {
    if (reliability_estimator == "alpha") {
      return(coef_alpha)
    }
    stop("Unknown reliability_estimator string: '", reliability_estimator,
         "'. Use 'alpha' or supply a function.")
  }
  stop("reliability_estimator must be the string 'alpha' or a function ",
       "that accepts a matrix or data frame of items and returns a ",
       "numeric scalar.")
}


#' Compute group-specific summary statistics
#'
#' Helper that returns a data frame of group sizes, means, observed
#' variances, fourth central moments, and (if items are provided)
#' reliability estimates from the supplied estimator.
#'
#' @param data A data frame.
#' @param group_var Name of the grouping variable.
#' @param score_var Name of the composite score variable.
#' @param items Optional character vector of item names for reliability.
#' @param reliability_estimator Function or "alpha". Used to compute
#'   reliability from the items matrix.
#' @param group_levels Optional named character vector with names
#'   "reference" and "focal" specifying the levels of group_var.
#' @return A data frame with one row per group: group, label, n, mean,
#'   var, mu4, alpha.
#' @noRd
group_summary <- function(data, group_var, score_var, items = NULL,
                          reliability_estimator = "alpha",
                          group_levels = NULL) {
  if (!is.null(group_levels)) {
    if (!all(c("reference", "focal") %in% names(group_levels))) {
      stop("group_levels must be a named character vector with names ",
           "'reference' and 'focal'.")
    }
    ref_label <- group_levels[["reference"]]
    foc_label <- group_levels[["focal"]]
    if (!all(c(ref_label, foc_label) %in% data[[group_var]])) {
      stop("Specified group_levels not found in group_var.")
    }
    levels_use <- c(ref_label, foc_label)
  } else {
    levels_use <- sort(unique(data[[group_var]]))
    if (length(levels_use) != 2) {
      stop("group_var must have exactly two unique values, or you must ",
           "supply group_levels.")
    }
  }

  rel_fun <- resolve_reliability_estimator(reliability_estimator)

  out <- data.frame(
    group = c("0_reference", "1_focal"),
    label = levels_use,
    n = NA_integer_,
    mean = NA_real_,
    var = NA_real_,
    mu4 = NA_real_,
    alpha = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(levels_use)) {
    sub <- data[data[[group_var]] == levels_use[i], , drop = FALSE]
    score <- sub[[score_var]]
    out$n[i] <- sum(!is.na(score))
    if (out$n[i] < 2) {
      stop("Group '", levels_use[i], "' has fewer than 2 valid ",
           "observations for score_var '", score_var, "'.")
    }
    out$mean[i] <- mean(score, na.rm = TRUE)
    out$var[i]  <- stats::var(score, na.rm = TRUE)
    if (!is.finite(out$var[i]) || out$var[i] <= 0) {
      stop("Group '", levels_use[i], "' has zero or non-finite ",
           "variance for score_var '", score_var,
           "'. The LTG-SMD framework requires positive variability ",
           "in both groups.")
    }
    out$mu4[i] <- fourth_moment(score, na.rm = TRUE)
    if (!is.null(items) && length(items) >= 2) {
      missing_cols <- setdiff(items, names(sub))
      if (length(missing_cols) > 0) {
        stop("Item columns not found in data: ",
             paste(missing_cols, collapse = ", "))
      }
      out$alpha[i] <- rel_fun(sub[, items, drop = FALSE])
    }
  }
  out
}
