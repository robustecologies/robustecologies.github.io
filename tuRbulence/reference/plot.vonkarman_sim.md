# Plot a von Kármán turbulent flow simulation

Visualizes the von Kármán system using time series, phase space,
embedded attractor, or power spectral density representations.

## Usage

``` r
# S3 method for class 'vonkarman_sim'
plot(
  x,
  type = c("timeseries", "phase", "attractor", "psd"),
  var = c("theta", "x", "y", "z"),
  n_points = NULL,
  ...
)
```

## Arguments

- x:

  Object of class "vonkarman_sim" from
  [`vonkarman_sim`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md).

- type:

  Type of visualization:

  "timeseries"

  :   Time series of a state variable (default)

  "phase"

  :   2D phase space projection (x_m, y)

  "attractor"

  :   3D embedded attractor via peak embedding (plotly)

  "psd"

  :   Power spectral density

- var:

  State variable for time series or PSD: "theta" (default), "x", "y", or
  "z".

- n_points:

  Maximum points to display (subsampling for performance).

- ...:

  Additional arguments passed to internal plotting functions.

## Value

Returns the plot object invisibly (ggplot2 or plotly).

## Details

The "attractor" type uses peak embedding to reconstruct the strange
attractor in delay coordinates (\\\theta_m\\, \\\theta\_{m+1}\\,
\\\theta\_{m+2}\\), following the methodology of the original
experimental studies.

## See also

[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`print.vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_sim.md),
[`summary.vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_sim.md),
[`vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md),
[`vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_peaks.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Time series
plot(sim)
plot(sim, var = "y")

# Phase space
plot(sim, type = "phase")

# 3D embedded attractor
plot(sim, type = "attractor")

# Power spectral density
plot(sim, type = "psd")
} # }
```
