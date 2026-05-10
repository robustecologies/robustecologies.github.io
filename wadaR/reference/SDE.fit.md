# Fit an SDE model via MCMC

Convenience wrapper that bundles an SDE model specification with the
MCMC inference call. Internally calls
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
and wraps the result in an `sde_fit` S3 class that provides SDE-specific
plotting and summary methods.

## Usage

``` r
SDE.fit(
  sde,
  Algorithm = "NUTS",
  Iterations = 5000L,
  Thinning = 1L,
  Status = 100L,
  Specs = NULL,
  Chains = 1L,
  CPUs = 1L,
  ...
)
```

## Arguments

- sde:

  An object of class `sde_model`, as returned by
  [`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md) or
  one of the family constructors.

- Algorithm:

  Character string specifying the MCMC algorithm. Default `"NUTS"`.

- Iterations:

  Integer number of MCMC iterations. Default 5000.

- Thinning:

  Integer thinning interval. Default 1.

- Status:

  Integer controlling how often progress is printed (every `Status`
  iterations). Default 100.

- Specs:

  Optional list of algorithm specifications. If `NULL`, default specs
  are used.

- Chains:

  Integer number of parallel chains. Default 1.

- CPUs:

  Integer number of CPUs for parallel chains. Default 1.

- ...:

  Additional arguments passed to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Value

An object of class `sde_fit`, inheriting from `demonoid`, with
additional components:

- sde:

  The original `sde_model` object.

- posterior_trajectories:

  Reserved for posterior predictive simulations (computed lazily by
  `plot.sde_fit`).

## Details

Fit an SDE model via MCMC

Implementation of `SDE.fit`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and fit an OU process
set.seed(42)
n <- 100; dt <- 0.1
x <- numeric(n); x[1] <- 5
for (i in 2:n) x[i] <- x[i-1] + 2*(5-x[i-1])*dt + 0.8*sqrt(dt)*rnorm(1)
y <- x + rnorm(n, 0, 0.3)
times <- seq(0, by = dt, length.out = n)

sde <- SDE(data = y, times = times, family = "ou")
fit <- SDE.fit(sde, Algorithm = "AMWG", Iterations = 30000,
    Thinning = 10)

print(fit)
summary(fit)
plot(fit, type = "trajectory")
plot(fit, type = "predictive")
} # }
```
