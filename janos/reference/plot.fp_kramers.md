# Plot method for fp_kramers objects

Visualizes the Kramers analysis. The default plot shows the potential
landscape with wells, barriers, and transition arrows annotated with
rates. The `"rate_vs_epsilon"` type sweeps over noise intensities to
show how the escape rate depends on \\\varepsilon\\.

## Usage

``` r
# S3 method for class 'fp_kramers'
plot(
  x,
  type = c("summary", "rate_vs_epsilon"),
  epsilon_range = NULL,
  n_epsilon = 50L,
  ...
)
```

## Arguments

- x:

  An `fp_kramers` object.

- type:

  Character string: `"summary"` (default) or `"rate_vs_epsilon"`.

- epsilon_range:

  Numeric vector of length 2 giving the range of \\\varepsilon\\ values
  for the sweep (default: computed automatically).

- n_epsilon:

  Number of epsilon values in the sweep (default 50).

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.
