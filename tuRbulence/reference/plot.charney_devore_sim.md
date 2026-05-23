# Plot a Charney-DeVore atmospheric blocking simulation

Visualizes the Charney-DeVore model using time series, phase space, 3D
attractor, or blocking regime representations.

## Usage

``` r
# S3 method for class 'charney_devore_sim'
plot(
  x,
  type = c("timeseries", "phase", "attractor", "blocking"),
  var = c("x", "y", "z", "wave"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  Object of class "charney_devore_sim" from
  [`charney_devore_sim`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md).

- type:

  Type of plot:

  "timeseries"

  :   Time series of selected variable (default)

  "phase"

  :   x-y phase space colored by time

  "attractor"

  :   3D plotly visualization of full attractor

  "blocking"

  :   Zonal flow vs wave amplitude, colored by regime

- var:

  Variable to plot for timeseries: "x" (zonal flow, default), "y" (wave
  mode 1), "z" (wave mode 2), or "wave" (wave amplitude).

- n_points:

  Maximum number of points to plot; if NULL, plots all points. Useful
  for large simulations.

- ...:

  Additional arguments passed to internal plot functions.

## Value

Returns the plot object invisibly (ggplot2 or plotly depending on type).

## Details

The Charney-DeVore model exhibits multiple equilibria: high-index
(strong zonal flow, weak waves) and low-index (blocking: weak zonal,
amplified waves). Different plot types reveal different aspects.

The phase plot shows trajectories in the x-y plane. The 3D attractor
visualization (requires plotly) shows the full three-dimensional
structure. The blocking plot identifies regime occupation based on zonal
flow strength and wave amplitude thresholds.

## See also

[`charney_devore_sim`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md)
for simulation,
[`charney_devore_batch`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md)
for parameter sweeps.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate atmospheric blocking dynamics
sim <- charney_devore_sim(F = 1.5, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Different visualization types
plot(sim)                     # Time series of zonal flow
plot(sim, type = "phase")     # x-y phase space
plot(sim, type = "attractor") # 3D attractor (plotly)
plot(sim, type = "blocking")  # Blocking regime identification
} # }
```
