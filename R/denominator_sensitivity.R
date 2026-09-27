#' Denominator-sensitivity profile across six symmetric denominators
#'
#' Computes the LTG-SMD estimate under each of six symmetric denominator
#' rules applied to the group-specific reference-sample true-score
#' standard deviations: geometric (default), arithmetic, harmonic, root
#' mean square, minimum, and maximum. This profile mirrors the six
#' symmetric denominators compared in Simulation D (Section 5.3 of the
#' companion paper; Supplementary Section D.4) and is recommended as
#' supplementary material (Supplementary Section F.2.4), particularly when
#' group-specific reference true-score variances differ markedly.
#'
#' @param x An "ltg_smd" object.
#' @return A data frame with one row per denominator rule, columns
#'   `denominator`, `D`, `estimate`. Class
#'   "denominator_sensitivity".
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
#' sens <- denominator_sensitivity(result)
#' print(sens)
denominator_sensitivity <- function(x) {
  if (!inherits(x, "ltg_smd")) {
    stop("x must be an 'ltg_smd' object.")
  }

  A_focal <- x$A_g[["focal"]]
  A_ref   <- x$A_g[["reference"]]
  sd_focal <- sqrt(A_focal)
  sd_ref   <- sqrt(A_ref)
  delta_mean <- x$study$mean[x$study$group == "1_focal"] -
                x$study$mean[x$study$group == "0_reference"]

  denoms <- list(
    geometric    = sqrt(sd_focal * sd_ref),       # GM of SDs = (A1 A0)^{1/4}
    arithmetic   = (sd_focal + sd_ref) / 2,
    harmonic     = 2 * sd_focal * sd_ref / (sd_focal + sd_ref),
    root_mean_sq = sqrt((sd_focal^2 + sd_ref^2) / 2),
    minimum      = pmin(sd_focal, sd_ref),
    maximum      = pmax(sd_focal, sd_ref)
  )

  out <- data.frame(
    denominator = names(denoms),
    D = unlist(denoms),
    estimate = delta_mean / unlist(denoms),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "note") <- paste(
    "All six denominators are symmetric in the group labels (label",
    "invariant). The geometric, arithmetic, harmonic, and RMS",
    "denominators satisfy HM <= GM <= AM <= RMS (with equality when",
    "group SDs are equal). The minimum and maximum anchor the effect",
    "to one extreme of the two group-specific SDs. Different denominators",
    "define different estimands; see Section 5.3 (Simulation D) of the paper."
  )
  class(out) <- c("denominator_sensitivity", "data.frame")
  out
}


#' @export
print.denominator_sensitivity <- function(x, digits = 3, ...) {
  cat("Denominator-sensitivity profile (six symmetric denominators)\n")
  cat("==========================================================\n\n")
  df <- x
  df$D <- round(df$D, digits)
  df$estimate <- round(df$estimate, digits)
  print.data.frame(df, row.names = FALSE)
  note <- attr(x, "note")
  if (!is.null(note)) {
    cat("\nNote: ", note, "\n", sep = "")
  }
  invisible(x)
}


#' Sensitivity of the LTG-SMD to the choice of target reference
#'
#' Recomputes the LTG-SMD under each of several user-supplied alternative
#' reference data frames. Useful when more than one target reference is
#' defensible (Section 8.3 of the companion paper).
#'
#' @param object An "ltg_smd" object computed with a primary reference.
#' @param alternative_references A named list of reference data frames.
#'   Names are used as labels in the output.
#' @param study_data Original study data (required to recompute).
#' @param group_var,score_var,items,group_levels,rho_external
#'   Arguments to pass through to [compute_ltg_smd()].
#' @return A data frame with columns `reference`, `ltg_smd`, `c_focal`,
#'   `c_reference`, `rho_focal`, `rho_reference`, including the primary
#'   reference in the first row.
#' @export
sensitivity_reference <- function(object,
                                   alternative_references,
                                   study_data,
                                   group_var,
                                   score_var = NULL,
                                   items = NULL,
                                   group_levels = NULL,
                                   rho_external = NULL) {
  if (!inherits(object, "ltg_smd")) {
    stop("object must be an 'ltg_smd' object.")
  }
  if (!is.list(alternative_references) ||
      is.null(names(alternative_references))) {
    stop("alternative_references must be a named list of data frames.")
  }

  results <- list()

  # Row for the primary reference (from the original object)
  results[["primary"]] <- data.frame(
    reference     = "primary",
    ltg_smd       = unname(object$point),
    c_focal       = unname(object$c_g[["focal"]]),
    c_reference   = unname(object$c_g[["reference"]]),
    rho_focal     = unname(object$rho_g[["focal"]]),
    rho_reference = unname(object$rho_g[["reference"]]),
    stringsAsFactors = FALSE
  )

  for (nm in names(alternative_references)) {
    alt <- alternative_references[[nm]]
    new_obj <- compute_ltg_smd(
      study_data     = study_data,
      reference_data = alt,
      group_var      = group_var,
      score_var      = score_var,
      items          = items,
      group_levels   = group_levels,
      rho_external   = rho_external
    )
    results[[nm]] <- data.frame(
      reference     = nm,
      ltg_smd       = unname(new_obj$point),
      c_focal       = unname(new_obj$c_g[["focal"]]),
      c_reference   = unname(new_obj$c_g[["reference"]]),
      rho_focal     = unname(new_obj$rho_g[["focal"]]),
      rho_reference = unname(new_obj$rho_g[["reference"]]),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
