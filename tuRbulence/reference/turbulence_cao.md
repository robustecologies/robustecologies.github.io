# Estimate embedding dimension using Cao's method

Determines the optimal embedding dimension by analyzing how false
neighbors decrease with increasing dimension.

## Usage

``` r
turbulence_cao(x, max_dim = 15L, tau = 10L, n_samples = 2000L, n_threads = 0L)
```

## Arguments

- x:

  Time series (numeric vector) or simulation object.

- max_dim:

  Maximum embedding dimension to test.

- tau:

  Time delay for embedding.

- n_samples:

  Number of reference points to sample.

- n_threads:

  Number of OpenMP threads.

## Value

Object of class "turbulence_cao" with dimension estimates.

## Details

Cao's method computes the ratio \\E_1(d)\\ which measures how the
distance to nearest neighbors changes when moving from dimension \\d\\
to \\d+1\\. For deterministic systems, \\E_1(d)\\ saturates near 1 when
the embedding dimension equals or exceeds the attractor dimension.

## References

Cao, L. (1997). Practical method for determining the minimum embedding
dimension of a scalar time series. *Physica D*, 110(1-2), 43-50.
[doi:10.1016/S0167-2789(97)00118-8](https://doi.org/10.1016/S0167-2789%2897%2900118-8)

## See also

[`print.turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_cao.md),
[`summary.turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_cao.md),
[`plot.turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_cao.md),
[`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md),
[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md),
[`turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate chaotic system
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Estimate embedding dimension
cao <- turbulence_cao(sim, max_dim = 15, tau = 10)
print(cao)

# Visualize E1* curve
plot(cao)
} # }
```
