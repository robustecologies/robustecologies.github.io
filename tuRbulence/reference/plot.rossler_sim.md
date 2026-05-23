# Plot a Rössler system simulation

Visualizes the Rössler system trajectory using time series, phase space,
or 3D attractor representations.

## Usage

``` r
# S3 method for class 'rossler_sim'
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

  Object of class "rossler_sim" from
  [`rossler_sim`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md).

- type:

  Type of plot: "timeseries" (default), "phase" (x-y plane), or
  "attractor" (3D plotly visualization).

- var:

  State variable to plot for time series: "x", "y", or "z".

- n_points:

  Maximum number of points to display (for performance).

- ...:

  Additional arguments (currently unused).

## Value

For ggplot2 plots, returns the plot object invisibly. For plotly 3D
plots, returns the plotly object invisibly.

## See also

[`rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md),
[`print.rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.rossler_sim.md),
[`summary.rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.rossler_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- rossler_sim(a = 0.2, b = 0.2, c = 5.7, n_steps = 50000, seed = 42)
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")
} # }
```
