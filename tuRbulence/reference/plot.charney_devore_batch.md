# Plot Charney-DeVore batch simulation results

Professional visualization for parameter sweep simulations of the
Charney-DeVore atmospheric blocking model.

## Usage

``` r
# S3 method for class 'charney_devore_batch'
plot(
  x,
  type = c("bifurcation", "density", "summary", "blocking"),
  var = c("x", "y", "z", "wave"),
  n_points = 2000,
  alpha = 0.15,
  ...
)
```

## Arguments

- x:

  Object of class "charney_devore_batch".

- type:

  Type of visualization: "bifurcation", "density", "summary", or
  "blocking" (blocking frequency vs forcing).

- var:

  State variable: "x" (zonal flow), "y", "z" (wave modes), "wave" (wave
  amplitude). Ignored for type = "blocking".

- n_points:

  Maximum points per parameter value.

- alpha:

  Point transparency.

- ...:

  Additional arguments.

## Value

A ggplot2 object.

## See also

[`charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md),
[`print.charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_batch.md),
[`summary.charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
F_seq <- seq(0.5, 3.0, length.out = 26)
batch <- charney_devore_batch(F_seq, n_steps = 100000, n_threads = 4)

plot(batch)
plot(batch, var = "wave", type = "summary")
plot(batch, type = "blocking")  # Blocking frequency curve
} # }
```
