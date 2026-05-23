# Plot von Kármán batch simulation results

Professional visualization of bifurcation diagrams and statistical
summaries for parameter sweep simulations of the von Kármán stochastic
Duffing model.

## Usage

``` r
# S3 method for class 'vonkarman_batch'
plot(
  x,
  type = c("bifurcation", "density", "summary"),
  var = c("x", "y"),
  n_points = 2000,
  alpha = 0.15,
  ...
)
```

## Arguments

- x:

  Object of class "vonkarman_batch" from
  [`vonkarman_batch`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_batch.md).

- type:

  Type of visualization:

  "bifurcation"

  :   Classic bifurcation diagram: state variable vs control parameter
      with point density revealing attractor structure

  "density"

  :   Ridge density plot showing probability distributions at each
      parameter value

  "summary"

  :   Statistical summary: mean ± sd ribbons with median line

- var:

  State variable to plot: "x" (Duffing position), "y" (velocity),
  "theta" (transformed observable). Default "x".

- n_points:

  Maximum points per parameter value for performance.

- alpha:

  Point transparency for bifurcation plot.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object.

## See also

[`vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_batch.md),
[`print.vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_batch.md),
[`summary.vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_batch.md),
[`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mu_seq <- seq(-0.3, 0.3, length.out = 31)
batch <- vonkarman_batch(mu_seq, n_steps = 100000, n_threads = 4)

plot(batch)                      # Bifurcation diagram
plot(batch, type = "density")    # Density ridges
plot(batch, type = "summary")    # Statistical summary
} # }
```
