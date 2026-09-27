# ltgsmd: Latent True-Score and Target-Population Anchored Geometric SMD

Companion R package to the methodological paper:

> Nakamura, D. (2026). The denominator chooses the estimand: A target-population true-score framework for standardized mean differences. *Psychological Methods*. Advance online publication. https://doi.org/10.1037/met0000875

## Installation

```r
# From CRAN:
install.packages("ltgsmd")

# Development version from OSF:
# https://doi.org/10.17605/OSF.IO/KW9R6
```

## Quick start

```r
library(ltgsmd)

# Single study with an internal holdout reference
result <- compute_ltg_smd(
  study_data     = study_df,
  reference_data = reference_df,
  group_var      = "condition",
  score_var      = "outcome",
  items          = c("item1", "item2", "item3", "item4"),
  group_levels   = c(reference = "control", focal = "treatment")
)
print(result)

# Confidence intervals
ci <- ltg_smd_ci(result, method = c("analytic", "bootstrap"),
                  B = 2000, seed = 20240501,
                  study_data = study_df, reference_data = reference_df,
                  group_var = "condition", score_var = "outcome",
                  items = c("item1", "item2", "item3", "item4"),
                  group_levels = c(reference = "control", focal = "treatment"))

# Denominator diagnostics (Table F.1 format)
diag <- denominator_diagnostics(result)

# Six-denominator sensitivity profile
sens <- denominator_sensitivity(result)

# Multi-site / meta-analysis
multisite <- multisite_ltg_smd(
  data               = ml_data,
  site_var           = "site",
  group_var          = "condition",
  items              = paste0("item", 1:25),
  reference_strategy = "pooled_across_sites",
  min_n_per_group    = 50
)
```

## What the LTG-SMD is

The latent true-score and target-population anchored geometric standardized mean difference (LTG-SMD) is

$$\delta_{\mathrm{LTG}} = \frac{\mu_{T1} - \mu_{T0}}{(\sigma^2_{T1,R}\,\sigma^2_{T0,R})^{1/4}},$$

a two-group SMD whose denominator is the geometric mean of group-specific true-score standard deviations in an explicitly chosen target reference population $R$. The framework treats the denominator as an *estimand* choice, in the spirit of Lundberg, Johnson, and Stewart (2021), rather than as a technical detail.

## Functions

| Function | Purpose |
|---|---|
| `compute_ltg_smd()` | Plug-in estimator + five SMD comparators (Hedges's g, Welch, observed geometric, external observed, LTG-SMD) |
| `ltg_smd_ci()` | Analytic delta-method + BC/BCa bootstrap confidence intervals |
| `denominator_diagnostics()` | Table F.1 diagnostics for reporting (Section 8.3) |
| `denominator_sensitivity()` | Profile across six symmetric denominators (Simulation D, Section 5.3) |
| `sensitivity_reference()` | Sensitivity to alternative target reference distributions |
| `multisite_ltg_smd()` | Multi-site wrapper with cross-site reference and meta-analytic pooling |
| `export_supplementary()` | Bundle diagnostics, CIs, and sensitivity for a supplementary appendix |
| `coef_alpha()` | Cronbach's coefficient alpha (internal default reliability estimator) |

## Vignette

The package vignette `vignette("getting-started", package = "ltgsmd")` works through the three Section 6 examples from the companion paper.

## Reproducing the empirical illustrations

The three empirical illustrations in Section 6 of the paper use:

1. **Example 1**: Ottmar et al. (2025) open mathematics-learning data (Journal of Open Psychology Data, https://doi.org/10.5334/jopd.139).
2. **Example 2**: Open Psychometrics IPIP Big Five 2014 dataset.
3. **Example 3**: Many Labs 2 Anderson.1 effect, Slate 1 data from Klein et al. (2018).

R scripts that reproduce each example are archived on OSF (https://doi.org/10.17605/OSF.IO/KW9R6); see also the data availability statement in the companion paper. The scripts read the public data and call this package's functions in the order documented in the vignette.

## Citation

Please cite the companion paper:

> Nakamura, D. (2026). The denominator chooses the estimand: A target-population true-score framework for standardized mean differences. *Psychological Methods*. Advance online publication. https://doi.org/10.1037/met0000875

In R, `citation("ltgsmd")` returns this reference, together with a citation for the package itself, in text and BibTeX form. Volume and page details will be added once the article is assigned to an issue.

## License

MIT. See the `LICENSE` file.

## Dependencies

- **R** (≥ 4.1.0)
- **stats** (base)
- **boot** (suggested; required for bootstrap confidence intervals)
- **metafor** (suggested; required for `multisite_ltg_smd()`)
- **psych** (suggested; an alternative reliability estimator)
- **testthat** (suggested; for unit tests)
- **knitr**, **rmarkdown** (suggested; for building the vignette)

## Status

This is version 0.2.3. The companion paper has been published in *Psychological Methods* (advance online publication, 2026), and the API is stable.
