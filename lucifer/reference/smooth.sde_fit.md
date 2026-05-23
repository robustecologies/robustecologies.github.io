# Forward-filtering backward-sampling smoother for SDE models

Runs a particle smoother on a fitted SDE model, producing smoothed state
estimates and sampled trajectories from the joint smoothing distribution
\\p(x\_{0:T} \| y\_{1:T}, \hat{\theta})\\.

## Usage

``` r
smooth.sde_fit(object, theta = NULL, N.particles = 500L, n.trajectories = 50L)
```

## Arguments

- object:

  An object of class `sde_fit`.

- theta:

  Optional parameter vector. If `NULL`, uses the posterior mean.

- N.particles:

  Integer number of particles. Default 500.

- n.trajectories:

  Integer number of backward-sampled trajectories. Default 50.

## Value

A list with components:

- times:

  Observation times.

- filtered_mean:

  Matrix of filtered state means.

- smoothed_mean:

  Matrix of smoothed state means.

- trajectories:

  Array (T x dim_x x n.trajectories) of sampled smoothed paths.

- loglik:

  Log-marginal-likelihood estimate.

## Details

Smooth latent states via particle smoother

Implementation of `smooth.sde_fit`. Refer to the package vignettes and
the cited references for a complete algorithmic and mathematical
description.

## See also

[`SDE.fit`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md),
[`LOO.sde_fit`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`WAIC.sde_fit`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`log_lik.sde_fit`](https://robustecologies.github.io/lucifer/reference/log_lik.md),
[`plot.sde_fit`](https://robustecologies.github.io/lucifer/reference/plot.sde_fit.md),
[`predict.sde_fit`](https://robustecologies.github.io/lucifer/reference/predict.sde_fit.md),
[`print.sde_fit`](https://robustecologies.github.io/lucifer/reference/print.sde_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving smooth.sde_fit
} # }
```
