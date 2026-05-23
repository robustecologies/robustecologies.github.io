# Plot Stommel batch simulation results

Professional visualization of bifurcation diagrams and statistical
summaries for parameter sweep simulations of the Stommel thermohaline
circulation model.

## Usage

``` r
# S3 method for class 'stommel_batch'
plot(
  x,
  type = c("bifurcation", "density", "summary"),
  var = c("q", "T", "S"),
  n_points = 2000,
  alpha = 0.15,
  ...
)
```

## Arguments

- x:

  Object of class "stommel_batch" from
  [`stommel_batch`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md).

- type:

  Type of visualization: "bifurcation", "density", or "summary".

- var:

  State variable: "q" (flow strength), "T" (temperature), "S"
  (salinity).

- n_points:

  Maximum points per parameter value.

- alpha:

  Point transparency.

- ...:

  Additional arguments.

## Value

A ggplot2 object.

## See also

[`stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md),
[`print.stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.stommel_batch.md),
[`summary.stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.stommel_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
eta2_seq <- seq(0.5, 2.0, length.out = 31)
batch <- stommel_batch(eta2_seq, n_steps = 100000, n_threads = 4)

plot(batch)                      # Bifurcation diagram
plot(batch, type = "density")    # Density ridges
} # }
```
