# Plot a Hasselmann stochastic climate simulation

Visualizes the Hasselmann climate model using time series, phase space,
power spectrum, or ice-temperature relationship representations.

## Usage

``` r
# S3 method for class 'hasselmann_sim'
plot(
  x,
  type = c("timeseries", "phase", "spectrum", "ice_temp"),
  var = c("T", "Td", "I"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  Object of class "hasselmann_sim" from
  [`hasselmann_sim`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md).

- type:

  Type of plot:

  "timeseries"

  :   Time series of selected variable (default)

  "phase"

  :   3D plotly visualization of T-Td-I space

  "spectrum"

  :   Power spectral density showing red noise character

  "ice_temp"

  :   Temperature-ice scatter with equilibrium curve

- var:

  Variable to plot for timeseries/spectrum: "T" (surface temperature,
  default), "Td" (deep ocean), or "I" (ice extent).

- n_points:

  Maximum number of points to plot; if NULL, plots all points. Useful
  for large simulations.

- ...:

  Additional arguments passed to internal plot functions.

## Value

Returns the plot object invisibly (ggplot2 or plotly depending on type).

## Details

The Hasselmann model produces red noise spectra characteristic of
climate variability, with the deep ocean acting as a low-pass filter.
Different plot types reveal different aspects of the stochastic climate
dynamics.

The spectrum plot shows power spectral density on log-log axes, with a
theoretical red noise reference (dashed line). The ice-temperature plot
shows the relationship between surface temperature and ice extent, with
the equilibrium curve I_eq(T) = max(0, -μT) shown as a dashed line.

The 3D phase plot (requires plotly) shows the trajectory in the full
T-Td-I state space, revealing the attractor structure.

## See also

[`hasselmann_sim`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md)
for simulation,
[`hasselmann_batch`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_batch.md)
for forcing scenarios.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate stochastic climate variability
sim <- hasselmann_sim(F_forcing = 0, n_steps = 200000, seed = 42)

# Inspect results
print(sim)
summary(sim)

# Different visualization types
plot(sim)                    # Time series of temperature
plot(sim, type = "phase")    # 3D phase space (plotly)
plot(sim, type = "spectrum") # Power spectrum
plot(sim, type = "ice_temp") # Ice-temperature relationship
} # }
```
