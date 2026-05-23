# Leave-future-out cross-validation for time series models

Computes approximate leave-future-out cross-validation (LFO-CV) using
Pareto smoothed importance sampling (PSIS), as described in Burkner,
Gabry and Vehtari (2020). Standard LOO-CV is theoretically inappropriate
for time series because leaving out a single observation still allows
future data to inform past predictions; LFO-CV avoids this by evaluating
predictions sequentially in time.

## Usage

``` r
LFO(
  x,
  Model = NULL,
  Data = NULL,
  L = NULL,
  M = 1L,
  k_threshold = 0.7,
  mode = c("forward", "backward"),
  refit_args = list(),
  verbose = TRUE
)
```

## Arguments

- x:

  Either a fitted object of class `demonoid`, `vb`, or `pmc`, or an \\N
  \times S\\ matrix of pointwise log-likelihoods ordered by time. When a
  matrix is supplied and `Model`/`Data` are `NULL`, only the PSIS
  approximation is used (no refitting).

- Model:

  Model specification function. Required for refitting when Pareto k
  exceeds the threshold.

- Data:

  Data list conforming to the lucifer data specification. Observations
  must be ordered by time.

- L:

  Minimum number of observations for the initial training window.
  Defaults to `max(10, 2 * n_params)`.

- M:

  Prediction horizon. `M = 1` gives one-step-ahead predictions
  (default); `M = 4` gives four-step-ahead.

- k_threshold:

  Pareto k threshold for triggering a model refit (default 0.7,
  following Burkner et al. 2020).

- mode:

  Character, either `"forward"` (default) or `"backward"`.

- refit_args:

  A list of overrides for refitting (e.g., `Iterations`, `Algorithm`).

- verbose:

  Logical; if `TRUE`, prints progress.

## Value

An object of class `lucifer_lfo`, which is a list with components:

- elpd_lfo:

  Total expected log pointwise predictive density.

- se_elpd_lfo:

  Standard error of `elpd_lfo`.

- pointwise:

  A data.frame with columns `time`, `elpd_lfo`, `pareto_k`, and
  `refitted`.

- n_refits:

  Number of model refits triggered.

- L:

  Initial training window size used.

- M:

  Prediction horizon.

- k_threshold:

  Pareto k threshold used.

- mode:

  CV mode used.

- call:

  The matched call.

## Details

The algorithm proceeds by expanding the training window forward through
time, starting from observation L. At each step, importance weights are
computed to approximate the posterior \\p(\theta \| y\_{1:i})\\ from the
most recent refit posterior \\p(\theta \| y\_{1:i^\*})\\. When the
Pareto k diagnostic of these importance weights exceeds `k_threshold`,
the model is refitted on \\y\_{1:i}\\ to restore reliability.

The predictive log-likelihood for future observation \\y\_{i+1}\\ is
computed by deviance differencing: evaluating the model on
\\y\_{1:i+1}\\ and \\y\_{1:i}\\ at each posterior sample, then taking
the difference. This is exact for models with additive log-likelihoods.

## References

Burkner, P.-C., Gabry, J., and Vehtari, A. (2020). "Approximate
leave-future-out cross-validation for Bayesian time series models."
*Journal of Statistical Computation and Simulation*, 90(14), 2499–2523.
[doi:10.1080/00949655.2020.1783262](https://doi.org/10.1080/00949655.2020.1783262)

Vehtari, A., Gelman, A., and Gabry, J. (2017). "Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC." *Statistics
and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`Kfold`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`loo_compare`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# AR(1) model
Model <- function(parm, Data) {
    phi <- parm[1]
    sigma <- exp(parm[2])
    mu <- phi * Data$y[-Data$n]
    LL <- sum(dnorm(Data$y[-1], mu, sigma, log = TRUE))
    LP <- LL + dnorm(phi, 0, 1, log = TRUE) +
        dgamma(sigma, 1, 1, log = TRUE)
    yhat <- c(Data$y[1], mu)
    list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}
Data <- list(y = arima.sim(list(ar = 0.7), 100), n = 100,
             parm.names = c("phi", "log.sigma"), mon.names = "LP")

fit <- lucifer(Model, Data, c(0, 0), Iterations = 10000,
    Algorithm = "HARM", Specs = NULL)

result <- LFO(fit, Model, Data, L = 20, M = 1)
print(result)
summary(result)
plot(result)
} # }
```
