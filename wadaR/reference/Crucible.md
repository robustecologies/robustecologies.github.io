# Automated Bayesian inference pipeline

Orchestrates the full lucifer inference pipeline: profiling via
[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
fitting the top-ranked methods, iterative refinement via
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
and comparison via
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md).
The result is a ranked set of converged fits with diagnostic metadata.

## Usage

``` r
Crucible(
  Model,
  Data,
  Initial.Values = NULL,
  spec = NULL,
  n_methods = 5,
  methods = NULL,
  max_rounds = 3,
  Chains = 4,
  CPUs = parallel::detectCores() - 1,
  families = c("MCMC", "VB", "Laplace", "IQ"),
  diverse = TRUE,
  reference = NULL,
  prescribe.args = list(),
  arena.args = list(),
  verbose = TRUE
)
```

## Arguments

- Model:

  A model specification function compatible with
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- Data:

  A named list of data.

- Initial.Values:

  A numeric vector of initial parameter values.

- spec:

  An optional `model_spec` object. If provided, Model, Data, and
  Initial.Values are extracted from it.

- n_methods:

  Integer. Number of top-ranked methods from
  [`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md)
  to try. Default 5. Ignored when `methods` is supplied.

- methods:

  Optional character vector of algorithm abbreviations to try,
  overriding Prescribe's automatic ranking. Each name must match one of
  the algorithms known to lucifer (e.g., `c("NUTS", "AMWG", "Slice")`).
  Prescribe is still called for profiling but the method selection is
  driven entirely by this vector.

- max_rounds:

  Integer. Maximum refinement rounds per method (1 = initial fit only).
  Default 3.

- Chains:

  Integer. Number of MCMC chains. Default 4.

- CPUs:

  Integer. CPU cores for MCMC parallelism. Default
  `parallel::detectCores() - 1`.

- families:

  Character vector of inference families to consider. Default
  `c("MCMC", "VB", "Laplace", "IQ")`. Ignored when `methods` is
  supplied.

- diverse:

  Logical. If `TRUE` (default), method selection ensures diversity
  across MCMC subcategories (gradient, ensemble, slice, geometric,
  multimodal, etc.) rather than selecting the top-N by score. This
  portfolio approach increases robustness when the profiling phase does
  not fully capture the posterior structure.

- reference:

  Reference distribution for Arena accuracy metrics. See
  [`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md).

- prescribe.args:

  Named list of extra arguments passed to
  [`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md).

- arena.args:

  Named list of extra arguments passed to
  [`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md).

- verbose:

  Logical. Default `TRUE`.

## Value

An object of class `crucible` containing:

- Prescription:

  The `prescription` object from Prescribe.

- Methods:

  Data frame of selected methods with scores.

- Results:

  Named list, one entry per method. Each contains `fit_final`,
  `consort`, `rounds` (list of fit/consort pairs), `n_rounds`,
  `appeased`, and `error`.

- Arena:

  The `arena` comparison object.

- Summary:

  Data frame with Method, Category, Rounds, Appeased, ESS.min, Minutes,
  and Rank.

- Best:

  List with `method`, `category`, `fit`, and `consort` for the
  top-ranked method.

- Config:

  List of configuration parameters.

- Call:

  The matched call.

- Minutes:

  Total elapsed time in minutes.

## Details

Crucible proceeds in six stages. First, it profiles the model using
Prescribe to rank all available inference methods. Second, it selects
the top `n_methods` methods (respecting the `families` filter) and fits
each one. Third, it evaluates each fit with Consort. Fourth, for any
non-appeased fit, it extracts Consort's structured suggestion and
re-fits, repeating up to `max_rounds` total. Fifth, it passes all final
fits to Arena for cross-method comparison. Sixth, it assembles the
results into a ranked summary.

MCMC methods are run with `Chains` parallel chains. Non-MCMC methods are
run sequentially. Each method's fitting is wrapped in error handling so
that individual failures do not abort the pipeline.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-Normalization, Folding, and Localization: An Improved Rhat
for Assessing Convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`Arena`](https://robustecologies.github.io/lucifer/reference/Arena.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnorm(mu, 0, 100, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rep(mu, Data$N)
  return(list(LP = LP, Dev = -2 * LL, Monitor = LP,
    yhat = yhat, parm = parm))
}
Data <- list(N = 50, y = rnorm(50, 5, 2),
  mon.names = "LP", parm.names = c("mu", "log.sigma"))
IV <- c(0, 0)

# Run automated pipeline
cr <- Crucible(Model, Data, IV, n_methods = 3, max_rounds = 2)

# Inspect
print(cr)
summary(cr)

# Visualize
plot(cr)
plot(cr, type = "arena")
plot(cr, type = "convergence")
} # }
```
