#' ltgsmd: Latent True-Score and Target-Population Anchored Geometric SMD
#'
#' Implements the LTG-SMD framework for two-group effect-size analysis.
#' The LTG-SMD is a standardized mean difference whose denominator is the
#' geometric mean of group-specific true-score standard deviations in an
#' explicitly chosen target reference population. The package supports
#' single-study analyses with internal holdout references, analyses with
#' external norm samples, and multi-site coordinated-replication
#' analyses with cross-site pooled references.
#'
#' @section Core functions:
#' \describe{
#'   \item{[compute_ltg_smd()]}{Compute the LTG-SMD and four conventional
#'     comparators (Hedges's g, Welch-type SMD, observed geometric SMD,
#'     external observed SMD) from study and reference data.}
#'   \item{[ltg_smd_ci()]}{Analytic delta-method and bias-corrected (BC) or
#'     BCa nonparametric bootstrap confidence intervals.}
#'   \item{[denominator_diagnostics()]}{Table F.1 diagnostics for reporting
#'     (Section 8.3 of the paper).}
#'   \item{[denominator_sensitivity()]}{Profile across six symmetric
#'     denominators (Simulation D; Section 5.3 of the paper).}
#'   \item{[sensitivity_reference()]}{Sensitivity to the choice of target
#'     reference.}
#'   \item{[multisite_ltg_smd()]}{Multi-site analysis with cross-site
#'     pooled reference and meta-analytic comparison.}
#'   \item{[export_supplementary()]}{Bundle diagnostics, CIs, and
#'     sensitivity profiles for a supplementary appendix.}
#' }
#'
#' @section Companion paper:
#' Nakamura, D. (2026). The denominator chooses the estimand: A
#' target-population true-score framework for standardized mean
#' differences. *Psychological Methods*. Advance online publication.
#' \doi{10.1037/met0000875}
#'
#' Run `citation("ltgsmd")` for this reference in text and BibTeX form.
#' See the package vignette
#' `vignette("getting-started", package = "ltgsmd")` for worked examples
#' corresponding to Sections 6.1, 6.2, and 6.3 of the paper.
#'
#' @docType package
#' @name ltgsmd-package
#' @aliases ltgsmd
"_PACKAGE"
