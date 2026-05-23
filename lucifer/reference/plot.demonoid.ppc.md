# Plots of posterior predictive checks for demonoid

Renders one of 40 posterior predictive check (PPC) plot styles for an
object of class `demonoid.ppc` using ggplot2 and the RElab theme. The
function is the canonical entry point for visual model criticism in
lucifer and is shared (via the internal `.plot_ppc()` engine) with the
analogous methods for `laplace.ppc`, `pmc.ppc`, `vb.ppc`, and
`iterquad.ppc` fits.

## Usage

``` r
# S3 method for class 'demonoid.ppc'
plot(x, Style = NULL, Data = NULL, Rows = NULL, PDF = FALSE, ...)
```

## Arguments

- x:

  A required object of class `demonoid.ppc` returned by
  [`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md).

- Style:

  Character string selecting the plot style. Defaults to `NULL`, which
  is treated as `"Density"`. The full set of supported styles is listed
  in *Details*.

- Data:

  The data list originally supplied to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
  Required for styles that need access to covariates (`"Covariates"`,
  `"Covariates, Categorical DV"`), to the response matrix
  (`"Fitted, Multivariate, *"`, `"Spatial*"`, `"Space-Time *"`, the
  multivariate `"Time-Series"` variants, and `"DW, Multivariate, C"`),
  or to spatial coordinates (`"Spatial"`, `"Spatial Uncertainty"`).
  Ignored otherwise.

- Rows:

  Optional integer vector of row indices in `x$y`/`x$yhat` to be
  plotted. Defaults to all rows. Several styles ignore this argument
  because they aggregate over rows (e.g. `"Density Overlay"`,
  `"ECDF Overlay"`, `"Stat"`, `"LOO-PIT"`, `"Mardia"`).

- PDF:

  Logical. If `TRUE`, the plots are saved to a PDF file instead of being
  printed to the active graphics device. The file name is style-specific
  and is written in the current working directory.

- ...:

  Additional arguments forwarded to the internal engine. The most useful
  are:

  `Group`

  :   Optional integer/factor grouping variable used by
      `Style = "Violin Grouped"` to split the posterior predictive
      distribution by group.

  `stat_fun`, `stat_fun2`

  :   Functions used by `Style = "Stat"` and `"Stat 2D"` to compute
      scalar test statistics over each replicated dataset (defaults:
      `mean` and [`stats::sd`](https://rdrr.io/r/stats/sd.html)).

  `loo`

  :   A `lucifer_loo` object used by `Style = "LOO-PIT"` to compute
      leave-one-out PIT values via importance weights.

  `forecast_start`

  :   Integer. For the time-series styles, the index that separates
      training observations from forecast horizon, drawn as a vertical
      reference line.

  `col`

  :   Character vector of colour overrides used by all styles. When
      `NULL` the RElab contrasting palette is used.

## Value

See Details.

## Details

**Plot styles by category.** The 40 styles fall into seven thematic
groups. Styles whose name ends in `"Multivariate, C"` facet by data
column (one panel per series), while `"Multivariate, R"` variants facet
by row.

*Univariate predictive distributions.* `"Density"` (default) plots the
posterior predictive density of `yhat[Rows, ]` with the observed
`y[Rows]` as a dashed reference line. `"Density Overlay"` and
`"Histogram Overlay"` superimpose multiple replicated datasets onto the
observed distribution. `"ECDF"` and `"ECDF Overlay"` are the
empirical-cumulative analogues.

*Discrepancy and goodness-of-fit.* `"Stat"` (univariate test statistic)
and `"Stat 2D"` (bivariate test statistics) compare the distribution of
a user-supplied scalar function applied to each replicated dataset
against the same statistic computed on `y`. `"Jarque-Bera"` and `"DW"`
(Durbin-Watson) compute classical normality and serial-correlation test
statistics on the residuals, respectively, with multivariate column
variants. `"Mardia"` reports Mardia's multivariate skewness and
kurtosis.

*Fitted-vs-observed.* `"Fitted"` plots posterior predictive medians
against observed values with 95% predictive intervals, a dashed 1:1
identity line, a LOESS smoother, and points coloured by absolute
residual. `"Fitted, Multivariate, C"` produces one such panel per
response column; `"Fitted, Multivariate, R"` per row.
`"Scatter Average"` plots the posterior predictive mean against observed
values, while `"Error Scatter"` and `"Error Histogram"` display the
residuals.

*Residual diagnostics.* `"Residuals"`, `"Residual Density"` and the
`"*, Multivariate, C"` / `"*, Multivariate, R"` variants render
time-ordered residuals, residual densities, and per-series breakdowns.

*Predictive intervals and counts.* `"Intervals"` draws a caterpillar of
per-observation 50% and 95% predictive intervals. `"Ribbon"` stacks them
into a continuous ribbon ordered by index. `"Bars"`, `"Rootogram"` and
`"Predictive Quantiles"` target count outcomes (the rootogram is the
canonical PPC for Poisson and negative binomial models).

*Covariate-aware checks.* `"Covariates"` and
`"Covariates, Categorical DV"` plot the posterior predictive mean (or
per-class probability) against each covariate column in `Data$X`, with a
LOESS smoother. These styles require `Data` so the covariate matrix can
be retrieved.

*Spatial and time-series checks.* `"Spatial"` and
`"Spatial Uncertainty"` require `Data$longitude` and `Data$latitude` (or
`Data$X`/`Data$Y` coordinates) and render the posterior predictive mean
and the width of the predictive interval as colour-coded site maps.
`"Space-Time by Space"` and `"Space-Time by Time"` build the analogous
animations as a ribbon-style spaghetti plot. `"Time-Series"` renders
observed and replicated time series in a single panel; the multivariate
variants facet by series. The `forecast_start` argument adds a vertical
line separating in-sample fit from out-of-sample forecast.

*Calibration.* `"LOO-PIT"` computes leave-one-out probability integral
transforms via PSIS importance weights and renders an ECDF against the
uniform reference, the gold standard for assessing whether a Bayesian
predictive distribution is calibrated.

*Grouped diagnostics.* `"Violin Grouped"` requires the `Group` argument
and splits the posterior predictive distribution into one violin per
group, useful for hierarchical models.

**Data requirements.** The styles that strictly require `Data` are: all
`"*, Multivariate, *"` variants; `"Covariates"`,
`"Covariates, Categorical DV"`; `"Spatial"`, `"Spatial Uncertainty"`;
`"Space-Time by Space"`, `"Space-Time by Time"`;
`"Time-Series, Multivariate, C"`, `"Time-Series, Multivariate, R"`. The
remaining styles operate purely on the `y` and `yhat` stored in `x`.

**The `"Fitted"` family in detail.** The `"Fitted"` style and its
multivariate variants are the workhorse for visual model criticism. They
render observed values on the x-axis and posterior predictive medians on
the y-axis, with vertical bars spanning the 2.5% and 97.5% quantiles,
points coloured by \\\|\hat{y} - y\_{obs}\|\\, a dashed 1:1 identity
line, and a LOESS smoother in colour. The plot region is forced
equal-aspect via `coord_cartesian()` so deviation from the identity is
geometrically meaningful, and a grey caption documents the visual
encoding. `"Fitted, Multivariate, C"` requires `Data$Y` and produces one
panel per series; `"Fitted, Multivariate, R"` produces one panel per
row.

## References

Gelman, A., Meng, X.-L., and Stern, H. (1996). "Posterior Predictive
Assessment of Model Fitness via Realized Discrepancies". *Statistica
Sinica*, 6, p. 733–807.

Durbin, J., and Watson, G.S. (1950). "Testing for Serial Correlation in
Least Squares Regression, I." *Biometrika*, 37, p. 409–428.

Jarque, C.M. and Bera, A.K. (1980). "Efficient Tests for Normality,
Homoscedasticity and Serial Independence of Regression Residuals".
*Economics Letters*, 6(3), p. 255–259.

Mardia, K.V. (1970). "Measures of Multivariate Skewness and Kurtosis
with Applications". *Biometrika*, 57(3), p. 519–530.

Vehtari, A., Gelman, A., and Gabry, J. (2017). "Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC". *Statistics
and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
[`summary.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/summary.demonoid.ppc.md),
[`ppc_dens_overlay`](https://robustecologies.github.io/lucifer/reference/ppc_dens_overlay.md),
[`ppc_stat`](https://robustecologies.github.io/lucifer/reference/ppc_stat.md),
[`ppc_rootogram`](https://robustecologies.github.io/lucifer/reference/ppc_rootogram.md),
[`plot_imputed`](https://robustecologies.github.io/lucifer/reference/plot_imputed.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## A simple linear regression with covariates
set.seed(1)
N <- 200; J <- 3
X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
beta_true <- c(1.0, 2.0, -0.5); sigma_true <- 0.7
y <- as.vector(X %*% beta_true + rnorm(N, 0, sigma_true))

Data <- list(y = y, X = X, N = N, J = J,
             parm.names = c(paste0("beta[", 1:J, "]"), "sigma"),
             mon.names = "LP")
Model <- function(parm, Data) {
    beta <- parm[1:Data$J]
    sigma <- interval(parm[Data$J + 1], 1e-6, Inf)
    parm[Data$J + 1] <- sigma
    mu <- as.vector(Data$X %*% beta)
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + sum(dnorm(beta, 0, 100, log = TRUE)) +
        dhalfcauchy(sigma, 25, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = rnorm(Data$N, mu, sigma), parm = parm)
}
fit <- lucifer(Model, Data, c(rep(0, J), 1),
               Iterations = 2000, Algorithm = "NUTS", Chains = 1)
ppc <- predict(fit, Model = Model, Data = Data)

## Default density plot of the first 16 records
plot(ppc, Style = "Density", Rows = 1:16)

## The full Fitted family: univariate, with 1:1 line, LOESS, and
## residual-magnitude colour gradient
plot(ppc, Style = "Fitted")

## Convenience wrappers
ppc_dens_overlay(ppc)         # Style = "Density Overlay"
ppc_stat(ppc, stat_fun = sd)  # Style = "Stat" with sd as test statistic

## Aggregated PPC summaries
plot(ppc, Style = "Density Overlay")
plot(ppc, Style = "Intervals")
plot(ppc, Style = "Ribbon")

## Covariate-aware diagnostic
plot(ppc, Style = "Covariates", Data = Data)
} # }
```
