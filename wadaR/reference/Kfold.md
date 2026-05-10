# K-fold cross-validation for Bayesian models

Performs exact K-fold cross-validation by refitting the model K times,
each time holding out a different fold of observations and computing the
out-of-sample log predictive density for the held-out data.

## Usage

``` r
Kfold(
  Model,
  Data,
  fit = NULL,
  K = 10L,
  folds = NULL,
  fit_args = list(),
  CPUs = 1L,
  verbose = TRUE
)
```

## Arguments

- Model:

  A model specification function compatible with
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
  [`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
  or
  [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

- Data:

  A list containing the data, conforming to the standard lucifer data
  specification.

- fit:

  An existing fitted object (demonoid, vb, or pmc) used to extract
  settings and warm-start the refits. If `NULL`, settings must be
  provided via `fit_args`.

- K:

  Integer number of folds (default 10). Setting K equal to the number of
  observations gives exact leave-one-out cross-validation via refitting.

- folds:

  Optional integer vector of length N assigning each observation to a
  fold (values from 1 to K). If `NULL`, random balanced assignment is
  used.

- fit_args:

  A list of arguments passed to the inference backend, overriding
  settings extracted from `fit`. Relevant fields include `Algorithm`,
  `Iterations`, `Thinning`, `Specs`, `Method`.

- CPUs:

  Number of folds to fit in parallel (default 1, sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix or sequential fallback on Windows.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

An object of class `lucifer_kfold`, which is a list with components:

- elpd_kfold:

  Total expected log pointwise predictive density.

- se_elpd_kfold:

  Standard error of `elpd_kfold`.

- p_kfold:

  Effective number of parameters (lppd minus elpd_kfold).

- pointwise:

  A data.frame with columns `observation`, `elpd_kfold`, and `fold`.

- K:

  Number of folds.

- folds:

  Integer vector of fold assignments.

- call:

  The matched call.

## Details

K-fold cross-validation partitions the data into K non-overlapping
folds, refits the model K times, each time leaving out one fold as test
data, and evaluates the predictive log density for each held-out
observation. Unlike LOO-PSIS, K-fold CV is exact (no importance sampling
approximation) and works even when the Pareto k diagnostics indicate
that PSIS-LOO is unreliable.

The computational cost is K times the cost of fitting the model, making
it substantially more expensive than LOO-PSIS. Use LOO-PSIS first and
resort to K-fold only when diagnostics indicate problems.

The elpd is computed via deviance differencing: for each held-out
observation \\y_i\\, the pointwise log-likelihood is obtained by
evaluating the model on the training data plus \\y_i\\ and subtracting
the deviance on training data alone.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). "Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC." *Statistics
and Computing*, 27(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

## See also

[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`loo_compare`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
[`LFO`](https://robustecologies.github.io/lucifer/reference/LFO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model and data
Model <- function(parm, Data) {
    mu <- parm[1]
    sigma <- exp(parm[2])
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
        dgamma(sigma, 1, 1, log = TRUE)
    yhat <- rep(mu, Data$n)
    list(LP = LP, Dev = -2 * LL, Monitor = LP, yhat = yhat, parm = parm)
}
Data <- list(y = rnorm(50, 5, 2), n = 50,
             parm.names = c("mu", "log.sigma"), mon.names = "LP")

# Fit once
fit <- lucifer(Model, Data, Initial.Values = c(0, 0),
    Iterations = 10000, Algorithm = "HARM", Specs = NULL)

# 10-fold CV
result <- Kfold(Model, Data, fit = fit, K = 10)
print(result)
summary(result)
plot(result)
} # }
```
