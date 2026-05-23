# Plot method for SDE model objects

Produces diagnostic visualizations of the SDE model before fitting,
including observed data time series, phase portraits, and prior
predictive trajectories.

## Usage

``` r
# S3 method for class 'sde_model'
plot(
  x,
  type = c("data", "phase", "prior_predictive"),
  n.paths = 10,
  col = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `sde_model`.

- type:

  Character string specifying the plot type: `"data"` (default) for
  observed time series, `"phase"` for 2D phase portrait (multivariate
  data only), `"prior_predictive"` for simulated trajectories from the
  prior.

- n.paths:

  Integer number of prior predictive trajectories to simulate. Default
  10.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments passed to plotting functions.

## Value

The ggplot2 object, invisibly.

## Details

Implementation of `plot.sde_model`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md),
[`compile.sde_model`](https://robustecologies.github.io/lucifer/reference/compile.sde_model.md),
[`print.sde_model`](https://robustecologies.github.io/lucifer/reference/print.sde_model.md),
[`simulate.sde_model`](https://robustecologies.github.io/lucifer/reference/simulate.sde_model.md),
[`summary.sde_model`](https://robustecologies.github.io/lucifer/reference/summary.sde_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.sde_model
} # }
```
