# Plot a Lorenz-84 atmospheric model simulation

Visualizes the Lorenz-84 system trajectory using time series, phase
space, or 3D attractor representations.

## Usage

``` r
# S3 method for class 'lorenz84_sim'
plot(
  x,
  type = c("timeseries", "phase", "attractor"),
  var = c("x", "y", "z"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  Object of class "lorenz84_sim" from
  [`lorenz84_sim`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md).

- type:

  Type of plot: "timeseries" (default), "phase" (x-y plane), or
  "attractor" (3D plotly visualization).

- var:

  State variable to plot for time series: "x" (westerly wind), "y" (wave
  cosine), or "z" (wave sine).

- n_points:

  Maximum number of points to display (for performance).

- ...:

  Additional arguments (currently unused).

## Value

For ggplot2 plots, returns the plot object invisibly. For plotly 3D
plots, returns the plotly object invisibly.

## See also

[`lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md),
[`print.lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz84_sim.md),
[`summary.lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz84_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- lorenz84_sim(F = 8.0, G = 1.0, n_steps = 50000, seed = 42)
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")
} # }
```
