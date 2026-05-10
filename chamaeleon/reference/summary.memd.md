# Summary method for memd objects

Compute detailed diagnostics for a multivariate empirical mode
decomposition, including energy distribution across IMFs, cumulative
energy, and frequency separation ratios. These diagnostics help assess
decomposition quality and interpret the multiscale structure.

## Usage

``` r
# S3 method for class 'memd'
summary(object, ...)
```

## Arguments

- object:

  Object of class `memd`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.memd` containing:

- n_imfs:

  Number of extracted intrinsic mode functions.

- n_samples:

  Number of time samples.

- n_channels:

  Number of channels (embedding dimension).

- mode_table:

  Data frame with frequency, timescale, energy, and cumulative energy
  for each mode.

- freq_separation:

  Vector of frequency ratios between adjacent modes.

- mean_freq_ratio:

  Mean frequency separation ratio.

- residual_energy_prop:

  Proportion of total energy in the residual.

- quality_assessment:

  Character string assessing decomposition quality.

## Details

The energy of each IMF is computed as the sum of squared amplitudes
across all channels. Well-separated modes typically show frequency
ratios between 1.5 and 3.0 (ideally near 2, corresponding to octave
separation).

**Quality indicators:**

- Good separation: mean frequency ratio 1.5-3.0, low residual energy

- Mode mixing: low frequency ratios (\<1.3) indicate modes with
  overlapping frequency content

- High residual: if residual contains \>30% of energy, the decomposition
  may be incomplete

## See also

[`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
for computing the decomposition;
[`print.memd`](https://robustecologies.github.io/chamaeleon/reference/print.memd.md),
[`plot.memd`](https://robustecologies.github.io/chamaeleon/reference/plot.memd.md).
