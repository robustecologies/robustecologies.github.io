# Compare prior and posterior distributions

Produces faceted density plots comparing the prior and posterior
distributions for selected parameters, with optional ground truth
reference lines.

## Usage

``` r
prior_vs_posterior(
  x,
  prior,
  Parms = NULL,
  ground_truth = NULL,
  n_prior = 10000,
  col = NULL,
  ...
)
```

## Arguments

- x:

  An object from any lucifer inference function (class `demonoid`, `vb`,
  `pmc`, `smc`, `abc`, `sbi`, `laplace`, `iterquad`, or `bayesquad`).

- prior:

  Either:

  - A function `rprior(n)` returning an \\n \times K\\ matrix of prior
    draws.

  - A named list of density functions, e.g.,
    `list(mu = function(x) dnorm(x, 0, 10), sigma = ...)`.

- Parms:

  Character vector of parameter names to select (uses `grep` matching).
  Defaults to all parameters (up to 9).

- ground_truth:

  Optional named numeric vector of true parameter values.

- n_prior:

  Integer; number of prior samples when `prior` is a function (default
  10000).

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the list of ggplot objects.

## Details

Implementation of `prior_vs_posterior`. Refer to the package vignettes
and the cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lucifer(Model, Data, Initial.Values, ...)

# Using a sampling function
rprior <- function(n) cbind(rnorm(n, 0, 10), rexp(n, 1))
prior_vs_posterior(fit, rprior)

# Using density functions
prior_dens <- list(
    "beta[1]" = function(x) dnorm(x, 0, 10),
    "sigma"   = function(x) dgamma(x, 1, 1)
)
prior_vs_posterior(fit, prior_dens, ground_truth = c("beta[1]" = 2))
} # }
```
