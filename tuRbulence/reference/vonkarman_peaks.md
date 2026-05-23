# Extract peaks from time series for attractor embedding

Implements the peak embedding procedure described in the paper,
extracting local maxima (or minima) with a minimum time separation
constraint.

## Usage

``` r
vonkarman_peaks(
  sim,
  t = NULL,
  min_separation = 0.05,
  extract_maxima = TRUE,
  smooth_window = 0
)
```

## Arguments

- sim:

  Object of class "vonkarman_sim" or numeric vector.

- t:

  Time vector (if sim is numeric).

- min_separation:

  Minimum time between consecutive peaks.

- extract_maxima:

  Logical; TRUE for maxima, FALSE for minima.

- smooth_window:

  Window size for moving average smoothing (0 = no smoothing).

## Value

Object of class "vonkarman_peaks" with extracted peaks.

## Details

Peak embedding creates a Poincaré-like section by extracting local
extrema from the time series. The number of peaks depends on the
oscillation frequency and the min_separation constraint. For typical
simulations with dt = 0.01 and n_steps = 100000 (total time = 1000),
expect roughly total_time / min_separation peaks (e.g., ~20000 peaks
with min_separation = 0.05).

To increase the number of embedded points: (1) increase n_steps in the
simulation, (2) decrease min_separation, or (3) use
create_trajectory_3d() for direct phase space visualization without peak
extraction.

## See also

[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`print.vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_peaks.md),
[`summary.vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_peaks.md),
[`plot.vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_peaks.md),
[`vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate chaotic attractor
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)

# Extract peaks for embedding (smaller min_separation = more points)
peaks <- vonkarman_peaks(sim, min_separation = 0.05)
print(peaks)  # Shows number of extracted peaks

# Create embedded attractor
attr <- vonkarman_attractor(peaks, embed_dim = 3)
plot(attr)
} # }
```
