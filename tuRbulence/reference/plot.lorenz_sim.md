# Plot a Lorenz system simulation

Visualizes the Lorenz system trajectory using time series, phase space,
or 3D attractor representations.

## Usage

``` r
# S3 method for class 'lorenz_sim'
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

  Object of class "lorenz_sim" from
  [`lorenz_sim`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md).

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

[`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md),
[`print.lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz_sim.md),
[`summary.lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz_sim.md),
[`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- lorenz_sim(sigma = 10, rho = 28, n_steps = 50000, seed = 42)
plot(sim)
plot(sim, type = "phase")
plot(sim, type = "attractor")
} # }
```
