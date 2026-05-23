# Plot a Stommel box model simulation

Visualizes the Stommel thermohaline circulation model using time series,
phase space, flow regime, or density representations.

## Usage

``` r
# S3 method for class 'stommel_sim'
plot(
  x,
  type = c("timeseries", "phase", "flow", "density"),
  var = c("q", "T", "S"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  Object of class "stommel_sim" from
  [`stommel_sim`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md).

- type:

  Type of plot:

  "timeseries"

  :   Time series of selected variable (default)

  "phase"

  :   Temperature-salinity phase space colored by time

  "flow"

  :   Flow strength colored by circulation regime

  "density"

  :   Probability density of flow strength

- var:

  Variable to plot for timeseries: "q" (flow strength, default), "T"
  (temperature), or "S" (salinity).

- n_points:

  Maximum number of points to plot; if NULL, plots all points. Useful
  for large simulations.

- ...:

  Additional arguments passed to internal plot functions.

## Value

Returns the ggplot2 object invisibly.

## Details

The Stommel model exhibits bistability between thermally-dominated (q \>
0) and salinity-dominated (q \< 0) circulation regimes. Different plot
types reveal different aspects of the dynamics.

The phase plot shows trajectories in the T-S plane with the q = 0
diagonal marking the boundary between circulation regimes. The flow plot
colors points by their circulation state, making regime transitions
visible. The density plot reveals the probability distribution of flow
states, showing bimodality when both regimes are visited.

## See also

[`stommel_sim`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md)
for simulation,
[`stommel_batch`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md)
for parameter sweeps.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate bistable thermohaline circulation
sim <- stommel_sim(eta2 = 1.0, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Different visualization types
plot(sim)                    # Time series of flow strength
plot(sim, type = "phase")    # T-S phase space
plot(sim, type = "flow")     # Regime-colored flow
plot(sim, type = "density")  # Flow distribution
} # }
```
