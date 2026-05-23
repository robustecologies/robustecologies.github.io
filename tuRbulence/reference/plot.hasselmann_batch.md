# Plot Hasselmann batch simulation results

Professional visualization for parameter sweep simulations of the
Hasselmann stochastic climate model.

## Usage

``` r
# S3 method for class 'hasselmann_batch'
plot(
  x,
  type = c("bifurcation", "density", "summary"),
  var = c("T", "Td", "I"),
  n_points = 2000,
  alpha = 0.15,
  ...
)
```

## Arguments

- x:

  Object of class "hasselmann_batch".

- type:

  Type of visualization: "bifurcation", "density", or "summary".

- var:

  State variable: "T" (surface temperature), "Td" (deep ocean), "I" (ice
  extent).

- n_points:

  Maximum points per parameter value.

- alpha:

  Point transparency.

- ...:

  Additional arguments.

## Value

A ggplot2 object.

## See also

[`hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_batch.md),
[`print.hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.hasselmann_batch.md),
[`summary.hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.hasselmann_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
F_seq <- seq(-1.0, 2.0, length.out = 31)
batch <- hasselmann_batch(F_seq, n_steps = 100000, n_threads = 4)

plot(batch)
plot(batch, var = "I", type = "summary")
} # }
```
