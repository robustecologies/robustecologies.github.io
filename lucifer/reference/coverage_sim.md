# Monte Carlo coverage simulation

Evaluates the frequentist coverage of confidence or credible intervals
produced by different estimation methods. For each simulation replicate,
data are generated from a known data-generating process, the model is
fitted with each method, intervals are extracted, and coverage of the
true parameter values is recorded.

## Usage

``` r
coverage_sim(
  dgp,
  fit_methods,
  true_parms,
  n_sims = 200,
  conf_level = 0.95,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- dgp:

  a function with no arguments that returns a list with elements
  `Model`, `Data`, and `Initial.Values`.

- fit_methods:

  a named list of fitting functions, each with signature
  `function(Model, Data, Initial.Values)`.

- true_parms:

  a named numeric vector of true parameter values.

- n_sims:

  number of simulation replicates, default 200.

- conf_level:

  confidence level for intervals, default 0.95.

- seed:

  optional random seed for reproducibility.

- verbose:

  logical; if `TRUE`, print progress updates.

## Value

An object of class `coverage_sim` with elements:

- coverage:

  matrix of empirical coverage rates (methods x parameters)

- mean_width:

  matrix of mean interval widths

- bias:

  matrix of mean bias (estimate minus true value)

- rmse:

  matrix of root mean squared error

- n_sims:

  number of completed simulations

- n_failures:

  matrix of fitting failures per method

- conf_level:

  nominal confidence level

- true_parms:

  true parameter values

## Details

This function is computationally expensive. With `n_sims = 200` and
multiple methods involving MCMC, a single call may take hours. It is
intended for methodological investigations, not routine use.

The `dgp` function must return a list with three elements: `Model` (a
lucifer-compatible model function), `Data` (a data list), and
`Initial.Values` (a numeric vector). Each call to `dgp` should produce a
fresh simulated dataset.

Each element of `fit_methods` is a function with signature
`function(Model, Data, Initial.Values)` that returns a fitted object
accepted by
[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md).

## See also

[`freq_summary`](https://robustecologies.github.io/lucifer/reference/freq_summary.md),
[`confint_compare`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Define data-generating process
dgp <- function() {
    n <- 100
    x <- rnorm(n)
    y <- 1.5 + 0.8 * x + rnorm(n, 0, 0.5)
    ## ... return list(Model, Data, Initial.Values)
}

## Define fitting methods
methods <- list(
    laplace = function(M, D, IV) LaplaceApproximation(M, IV, D)
)

## Run coverage simulation
result <- coverage_sim(dgp, methods, true_parms = c(beta0=1.5, beta1=0.8, sigma=0.5))
print(result)
summary(result)
plot(result)
} # }
```
