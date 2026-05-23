# Plot bifurcation diagram

Visualizes the Lyapunov exponent as a function of the control parameter,
showing transitions between stable and chaotic regimes.

## Usage

``` r
# S3 method for class 'turbulence_bifurcation'
plot(x, ...)
```

## Arguments

- x:

  Object of class "turbulence_bifurcation" from
  [`turbulence_bifurcation`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md).

- ...:

  Additional arguments (currently unused).

## Value

Returns the ggplot2 object invisibly.

## Details

The plot displays the largest Lyapunov exponent λ₁ versus the control
parameter. A horizontal dashed line at λ₁ = 0 separates chaotic (λ₁ \>
0) from stable (λ₁ \< 0) regimes. The parameter label is automatically
selected based on the system type (μ for vonkarman, η₂ for stommel,
etc.).

Bifurcation points where the system transitions from stable to chaotic
behavior can be identified where the curve crosses zero.

## See also

[`turbulence_bifurcation`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)
for the analysis function,
[`turbulence_lyapunov`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md)
for single-parameter estimation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute bifurcation diagram for von Kármán system
bif <- turbulence_bifurcation(
    system = "vonkarman",
    param_values = seq(0.1, 0.5, by = 0.05),
    n_steps = 50000,
    verbose = TRUE
)
print(bif)
summary(bif)
plot(bif)
} # }
```
