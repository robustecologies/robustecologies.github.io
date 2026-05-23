# Plot Cao's method results for embedding dimension

Visualizes the E1\* curve from Cao's method showing how the embedding
dimension saturates, with the estimated optimal dimension marked.

## Usage

``` r
# S3 method for class 'turbulence_cao'
plot(x, ...)
```

## Arguments

- x:

  Object of class "turbulence_cao" from
  [`turbulence_cao`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md).

- ...:

  Additional arguments (currently unused).

## Value

Returns the ggplot2 object invisibly.

## Details

The plot shows the E1\* ratio as a function of embedding dimension. When
E1\* approaches 1, false nearest neighbors are eliminated and the
embedding dimension is sufficient. The vertical dashed line indicates
the estimated optimal dimension where E1\* first exceeds 0.95.

For deterministic systems, E1\* should saturate clearly. For stochastic
systems or noisy data, saturation may be gradual.

## See also

[`turbulence_cao`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md)
for the estimation function.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)
cao <- turbulence_cao(sim, max_dim = 15, tau = 10)
print(cao)
plot(cao)
} # }
```
