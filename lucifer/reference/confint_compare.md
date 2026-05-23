# Compare credible and confidence intervals across fits

Creates a forest plot comparing intervals from multiple fit objects.
This is useful for visually assessing the agreement between Bayesian
credible intervals and frequentist confidence intervals, or between
different estimation methods.

## Usage

``` r
confint_compare(fits, true_values = NULL, prob = 0.95, parm = NULL, col = NULL)
```

## Arguments

- fits:

  a named list of fitted model objects (any class accepted by
  [`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)).

- true_values:

  optional named numeric vector of true parameter values, overlaid as
  reference points.

- prob:

  probability level for intervals, default 0.95.

- parm:

  character vector of parameter names to include. If `NULL`, all
  parameters common to all fits are shown.

- col:

  optional named list of color overrides passed to `.plot_colors()`.

## Value

A `ggplot2` object (invisible).

## Details

Implementation of `confint_compare`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- list(
    Laplace = laplace_fit,
    MCMC    = mcmc_fit
)
confint_compare(fits, true_values = c(beta1 = 1.0, beta2 = -0.5))
} # }
```
