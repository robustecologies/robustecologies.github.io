# Particle MCMC for SDE models

Convenience wrapper for pseudo-marginal MCMC using a particle filter
likelihood estimate. Recommends HARM or RWM (non-gradient algorithms)
and allows selection of the particle filter type (bootstrap, auxiliary,
or bridge). The particle filter produces an unbiased estimate of the
marginal likelihood, so the resulting MCMC chain targets the exact
posterior (Andrieu, Doucet, Holenstein 2010).

## Usage

``` r
SDE.fit.pmmh(
  sde,
  Algorithm = "HARM",
  Iterations = 10000L,
  N.particles = 500L,
  pf.type = "bootstrap",
  Thinning = 5L,
  Specs = NULL,
  ...
)
```

## Arguments

- sde:

  An object of class `sde_model` with `method = "particle"`.

- Algorithm:

  Character string. Default `"HARM"`, which is recommended for
  pseudo-marginal MCMC.

- Iterations:

  Integer number of MCMC iterations. Default 10000.

- N.particles:

  Integer number of particles. More particles reduce the variance of the
  likelihood estimate but increase cost. Default 500.

- pf.type:

  Character string: `"bootstrap"` (default), `"auxiliary"`, or
  `"bridge"`.

- Thinning:

  Integer thinning interval. Default 5.

- Specs:

  Optional list of algorithm specifications.

- ...:

  Additional arguments passed to
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Value

An object of class `sde_fit`.

## Details

Fit an SDE model via particle MCMC

Implementation of `SDE.fit.pmmh`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## References

Andrieu, C., Doucet, A., and Holenstein, R. (2010). "Particle Markov
chain Monte Carlo methods." *JRSS-B*, 72(3), p. 269–342.
[doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate predator-prey data
set.seed(42)
times <- seq(0, 10, by = 0.5)
prey <- rpois(length(times), 50)
pred <- rpois(length(times), 30)

sde <- SDE(data = cbind(prey, pred), times = times,
           family = "lotka_volterra")
fit <- SDE.fit.pmmh(sde, N.particles = 300, Iterations = 5000)

print(fit)
summary(fit)
plot(fit, type = "trajectory")
} # }
```
