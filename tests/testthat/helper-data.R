# Helper functions for generating test data with realistic factor
# structure. testthat sources files starting with "helper-" before
# running tests, so make_factor_data() is available in every test.
#
# Why this matters:
#   The first release of v0.2.0 used `rnorm()` to generate items
#   independently. That gave alpha ~ 0 (often negative), which the new
#   reliability range check correctly rejects. The fix is to generate
#   items with a shared latent factor, which is also what real
#   psychometric data look like.

#' Generate two-group test data with a shared latent factor
#'
#' Produces a data frame with `group` (a two-level factor) and `k`
#' item columns, where each item is generated as
#'     item_ij = loading * latent_i + group_effect_i + eps_ij
#' with latent ~ N(0,1), eps ~ N(0, sigma_err), and group_effect
#' equal to `effect` for the focal group, 0 otherwise.
#'
#' With loading = 1 and sigma_err = 0.5, four items give alpha ~ 0.94.
#' With loading = 1 and sigma_err = 1.0, four items give alpha ~ 0.80.
#'
#' @param n_per_group Integer.
#' @param k Number of items.
#' @param effect Group mean shift on each item.
#' @param loading Common loading (default 1).
#' @param sigma_err Item-specific noise SD (default 0.5; alpha ~ 0.94 at k=4).
#' @param group_levels Length-2 character vector of group labels.
#' @param group_name Name for the grouping column.
#' @param item_prefix Prefix for item column names.
#' @param compute_composite If TRUE, also add a row-mean composite
#'   column named "outcome".
#' @param seed Optional integer seed.
#' @return A data frame.
make_factor_data <- function(n_per_group = 100,
                              k = 4,
                              effect = 0.5,
                              loading = 1,
                              sigma_err = 0.5,
                              group_levels = c("a", "b"),
                              group_name = "condition",
                              item_prefix = "item",
                              compute_composite = TRUE,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- 2L * n_per_group
  group <- rep(group_levels, each = n_per_group)
  latent <- stats::rnorm(N)
  effect_vec <- ifelse(group == group_levels[2], effect, 0)

  df <- data.frame(g = group, stringsAsFactors = FALSE)
  names(df)[1] <- group_name

  for (i in seq_len(k)) {
    eps <- stats::rnorm(N, sd = sigma_err)
    df[[paste0(item_prefix, i)]] <-
      loading * latent + effect_vec + eps
  }

  if (compute_composite) {
    items <- paste0(item_prefix, seq_len(k))
    df$outcome <- rowMeans(df[, items, drop = FALSE])
  }

  df
}
