# Plot Lyapunov exponent estimation results

Visualizes the divergence curve and linear fit used to estimate the
largest Lyapunov exponent.

## Usage

``` r
# S3 method for class 'turbulence_lyapunov'
plot(x, ...)
```

## Arguments

- x:

  Object of class "turbulence_lyapunov" from
  [`turbulence_lyapunov`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md).

- ...:

  Additional arguments (currently unused).

## Value

Returns the ggplot2 object invisibly.

## Details

The plot shows the mean log divergence of nearby trajectories over time.
The slope of the linear region (shown as a dashed red line) gives the
Lyapunov exponent estimate. Positive slopes indicate chaos.

## See also

[`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md),
[`print.turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_lyapunov.md),
[`summary.turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_lyapunov.md),
[`turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- lorenz_sim(rho = 28, n_steps = 50000, seed = 42)
lyap <- turbulence_lyapunov(sim, embed_dim = 10, tau = 10)
print(lyap)
plot(lyap)
} # }
```
