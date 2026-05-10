# Posterior imputation of missing time series observations

Publication-quality visualization of Bayesian missing-data imputation
for univariate and multivariate time series. Given a fitted lucifer
model whose parameter vector includes one latent variable per missing
observation (the data-augmentation idiom of Tanner and Wong), the
function locates the imputed-value parameters in the posterior, matches
them positionally to the `NA` cells in the original data, and renders
the time series with observed values, posterior median imputations,
credible intervals, and (optionally) ground truth.

## Usage

``` r
plot_imputed(
  fit,
  data,
  Style = "default",
  indicators = NULL,
  time = NULL,
  series_names = NULL,
  ground_truth = NULL,
  level = 0.95,
  n_draws = 100,
  col = NULL,
  PDF = FALSE,
  FileName = "PPC.Imputed.pdf"
)
```

## Arguments

- fit:

  A fitted lucifer object of class `demonoid`, `vb`, `pmc`, `smc`,
  `abc`, `sbi`, `laplace`, `iterquad`, or `bayesquad`. A bare posterior
  matrix with named columns is also accepted.

- data:

  Original data containing `NA` values at the missing positions. Must be
  a numeric vector (univariate time series) or a numeric matrix / data
  frame (multivariate, rows = time, columns = series). Imputed
  parameters are matched to `NA` cells in column-major order, which is
  the convention used throughout lucifer.

- Style:

  Character string selecting the visualization style. One of:

  `"default"`

  :   Vertical credible-interval bars and hollow diamonds for the
      posterior median at each missing position (default).

  `"ribbon"`

  :   Continuous shaded ribbon spanning the credible interval, with a
      dotted median line passing through the imputed positions. Best
      when the missing fraction is concentrated in contiguous gaps.

  `"draws"`

  :   Spaghetti overlay of `n_draws` random posterior draws threading
      through the missing positions, revealing the joint correlation
      structure of nearby imputations.

- indicators:

  Optional. Character vector of parameter name patterns (matched via
  `grep`) or integer vector of column indices in the posterior matrix
  that correspond to the imputed-value parameters. When `NULL`
  (default), columns whose names contain `"imp"` or `"miss"` (case
  insensitive) are auto-detected.

- time:

  Optional numeric vector of time points. Defaults to `seq_len(N)` for
  univariate input or `seq_len(nrow(data))` for multivariate input.

- series_names:

  Optional character vector of facet labels for the multivariate case.
  Defaults to `colnames(data)` when present.

- ground_truth:

  Optional numeric vector or matrix of true values (same shape as
  `data`) used to overlay X markers at the imputed positions. Useful in
  simulation studies where the latent values are known. For univariate
  input may also be a vector of length equal to the number of missing
  observations.

- level:

  Credible interval level for the posterior summaries. Default is
  `0.95`.

- n_draws:

  Integer. Number of random posterior draws used by `Style = "draws"`.
  Default is `100`.

- col:

  Optional character vector of colors. When `NULL` the RElab contrasting
  palette is used.

- PDF:

  Logical. If `TRUE`, the plot is written to `FileName` instead of being
  printed.

- FileName:

  File path used when `PDF = TRUE`.

## Value

Invisibly returns the underlying `ggplot2` object so that downstream
code may add layers, modify scales, or save it via
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The function performs Bayesian missing-data visualization in the
data-augmentation framework. The user fits a lucifer model whose
parameter vector has been augmented with one latent variable per missing
observation; lucifer's auto-FC Gibbs sampler (or any other inference
engine) draws from the joint posterior over parameters and imputations
simultaneously. `plot_imputed()` extracts the marginal posterior of each
imputed value, summarises it with the requested credible interval, and
overlays the result on the observed time series.

The matching between imputed parameter columns and `NA` cells in `data`
is positional and follows R's column-major ordering for matrices, which
is the convention adopted throughout lucifer (e.g. by
[`as.parm.names()`](https://robustecologies.github.io/lucifer/reference/as.parm.names.md)).
When the user follows this convention the matching requires no
additional argument; otherwise the `indicators` argument can be used to
provide an explicit set of column indices.

Each `Style` highlights a different aspect of the posterior imputation.
`"default"` is the most informative single view: each missing position
carries a vertical bar for the credible interval and a hollow diamond
for the median. `"ribbon"` converts the bars into a continuous band,
which makes contiguous gaps visually prominent and is the appropriate
choice for sensor dropouts or structured missingness. `"draws"` shows
multiple posterior draws simultaneously, which reveals correlations
across adjacent imputations that the marginal credible intervals
collapse away.

For multivariate input the function facets by series and produces a
self-contained subtitle reporting the total counts of observed and
imputed cells across all series. A grey caption in the bottom-right
summarises the visual encoding of the chosen `Style`.

## References

Tanner, M.A. and Wong, W.H. (1987). The calculation of posterior
distributions by data augmentation. *Journal of the American Statistical
Association*, 82(398), 528–540.
[doi:10.1080/01621459.1987.10478458](https://doi.org/10.1080/01621459.1987.10478458)

Little, R.J.A. and Rubin, D.B. (2019). *Statistical Analysis with
Missing Data*, 3rd ed. Wiley.
[doi:10.1002/9781119482260](https://doi.org/10.1002/9781119482260)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`plot.demonoid.ppc`](https://robustecologies.github.io/lucifer/reference/plot.demonoid.ppc.md),
[`caterpillar.plot`](https://robustecologies.github.io/lucifer/reference/caterpillar.plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Univariate AR(1) with random missingness
set.seed(42)
n <- 120
phi <- 0.8; sigma <- 0.5
y_complete <- numeric(n)
for (t in 2:n) y_complete[t] <- phi * y_complete[t - 1] + rnorm(1, 0, sigma)
miss <- sort(sample(2:n, 25))
y_obs <- y_complete; y_obs[miss] <- NA

parm.names <- c("phi", "sigma", paste0("y_imp_", miss))
Data <- list(y_obs = y_obs, n = n, miss = miss, n_miss = length(miss),
             parm.names = parm.names, mon.names = "LP",
             pos.phi = 1, pos.sigma = 2, pos.imp = 3:length(parm.names))
Model <- function(parm, Data) {
    phi <- parm[Data$pos.phi]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    y <- Data$y_obs; y[Data$miss] <- parm[Data$pos.imp]
    LL <- sum(dnorm(y[-1], phi * y[-Data$n], sigma, log = TRUE))
    LP <- LL + dnorm(phi, 0, 1, log = TRUE) +
        dhalfcauchy(sigma, 5, log = TRUE)
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = c(0, phi * y[-Data$n] + rnorm(Data$n - 1, 0, sigma)),
         parm = parm)
}
iv <- c(0.5, 0.5, rep(0, length(miss)))
fit <- lucifer(Model, Data, iv, Iterations = 2000,
               Algorithm = "Gibbs", Chains = 1)

## Three styles
plot_imputed(fit, data = y_obs)                          # default
plot_imputed(fit, data = y_obs, Style = "ribbon")
plot_imputed(fit, data = y_obs, Style = "draws", n_draws = 80)

## With ground truth overlay
plot_imputed(fit, data = y_obs, ground_truth = y_complete)

## Multivariate: a 3-series VAR with missingness in two of the series
Y <- cbind(y_complete,
           y_complete + rnorm(n, 0, 0.3),
           cumsum(rnorm(n, 0, 0.2)))
Y[miss, 1] <- NA
Y[sample(n, 15), 3] <- NA
colnames(Y) <- c("temp", "humid", "pressure")
## Fit a multivariate model that imputes Y[i,j] for each NA cell, then:
# plot_imputed(fit_mv, data = Y, Style = "ribbon",
#              series_names = colnames(Y))
} # }
```
