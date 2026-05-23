# Plot method for SDE fit objects

Produces SDE-specific visualizations of a fitted model, including
posterior predictive trajectories, phase portraits, residuals, and
filtered state estimates.

## Usage

``` r
# S3 method for class 'sde_fit'
plot(
  x,
  type = c("trajectory", "predictive", "phase", "residuals", "posterior", "filtered",
    "interaction"),
  n.paths = 20,
  col = NULL,
  components = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `sde_fit`.

- type:

  Character string specifying the plot type: `"trajectory"` (default)
  for observed data with posterior mean filtered trajectory and
  observation-level 95% credible bands, `"predictive"` for posterior
  replicated observations overlaid on data (filtered states plus
  observation noise), `"phase"` for phase portrait with posterior
  predictive trajectories, `"residuals"` for one-step-ahead prediction
  residuals (data minus filtered states), `"posterior"` for standard
  MCMC trace/density plots (delegates to `plot.demonoid`), `"filtered"`
  for filtered state estimates.

- n.paths:

  Integer number of posterior predictive trajectories. Default 20.

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- components:

  Optional character or integer vector selecting which state components
  to plot for multivariate SDE models. Defaults to `NULL` (plot all
  state components).

- ...:

  Additional arguments passed to downstream plot methods.

## Value

The ggplot2 object, invisibly. For `type = "posterior"`, returns
`invisible(NULL)`.

## Details

Produces diagnostic plots of a stochastic differential equation fit
produced by
[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Iacus, S. M. (2008). *Simulation and Inference for Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-0-387-75839-8](https://doi.org/10.1007/978-0-387-75839-8)

## See also

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md),
[`summary.sde_fit`](https://robustecologies.github.io/lucifer/reference/summary.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and fit OU
set.seed(42)
n <- 100; dt <- 0.1
x <- numeric(n); x[1] <- 5
for (i in 2:n) x[i] <- x[i-1] + 2*(5-x[i-1])*dt + 0.8*sqrt(dt)*rnorm(1)
y <- x + rnorm(n, 0, 0.3)
times <- seq(0, by = dt, length.out = n)

sde <- SDE(data = y, times = times, family = "ou")
fit <- SDE.fit(sde, Algorithm = "AMWG", Iterations = 30000,
    Thinning = 10)

# Posterior filtered trajectory with observation-level CI
plot(fit, type = "trajectory")

# Posterior replicated observations (PPC)
plot(fit, type = "predictive", n.paths = 30)

# One-step-ahead prediction residuals
plot(fit, type = "residuals")

# Filtered state estimates
plot(fit, type = "filtered")
} # }
```
