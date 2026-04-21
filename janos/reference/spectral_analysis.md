# Power spectral density estimation

Estimates the power spectral density of state variable time series on
the attractor using Welch's method, providing frequency-domain
characterization of system dynamics. This complements ACF-based period
estimation with a more reliable approach for distinguishing
quasi-periodic from chaotic dynamics.

## Usage

``` r
spectral_analysis(
  x,
  variables = NULL,
  spans = c(5, 5),
  detrend = TRUE,
  log_scale = TRUE
)
```

## Arguments

- x:

  A `dyn_sim` object.

- variables:

  Character vector of state variables to analyze. If NULL, uses the
  model's state names.

- spans:

  Odd integer vector for Daniell smoothing kernel widths (default:
  `c(5, 5)`). Larger values produce smoother spectra.

- detrend:

  Logical; remove linear trend before spectral estimation (default:
  TRUE).

- log_scale:

  Logical; compute log10 of PSD for display (default: TRUE).

## Value

An S3 object of class `spectral_analysis` containing:

- spectra:

  Named list of `spec` objects from
  [`stats::spec.pgram`](https://rdrr.io/r/stats/spec.pgram.html), one
  per variable

- dominant_frequencies:

  Data frame with dominant frequency (cycles per time unit), period
  (time units), and peak-to-background ratio for each variable

- parameters:

  Model parameters from the input object

- metadata:

  Analysis settings and diagnostics

## Details

The PSD is computed via
[`stats::spec.pgram`](https://rdrr.io/r/stats/spec.pgram.html) with
Daniell smoothing (modified Welch periodogram). For chaotic systems, the
spectrum is broadband with no sharp peaks; for periodic systems, the
spectrum concentrates at the fundamental frequency and its harmonics;
for quasi-periodic systems, the spectrum shows multiple incommensurate
peaks.

The spectral classification is driven by the normalised spectral entropy
\\H = -\sum_k p_k \log p_k / \log N_f\\, with \\p_k = \Phi_k / \sum_j
\Phi_j\\, the distribution of power across frequencies. \\H\\ is bounded
in \\\[0, 1\]\\: close to zero for a spectrum concentrated on a single
line (periodic), close to one for a flat broadband spectrum, and
intermediate for narrowband mixtures. Peak-to-background ratios and peak
frequencies are reported as descriptive diagnostics, but are not used as
classifiers because the median background level is not informative for
heavy-tailed (\\1/f^\alpha\\) chaotic spectra, where it can underflow by
many orders of magnitude. Note that a narrowband spectrum is consistent
with both a noisy limit cycle and a Rossler-type chaotic attractor with
a dominant fundamental; definitive chaos detection requires a positive
Lyapunov exponent via
[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
or a high `K` statistic via
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md).

## References

Welch, P. D. (1967). The use of fast Fourier transform for the
estimation of power spectra: a method based on time averaging over
short, modified periodograms. *IEEE Transactions on Audio and
Electroacoustics*, 15(2), 70-73.
[doi:10.1109/TAU.1967.1161901](https://doi.org/10.1109/TAU.1967.1161901)

## See also

[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md),
[`lyapunov_spectrum`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md),
[`zero_one_test`](https://robustecologies.github.io/janos/reference/zero_one_test.md),
[`correlation_dimension`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate Lotka-Volterra dynamics
model <- model_spec(
    rhs = list(
        prey ~ alpha * prey - beta * prey * predator,
        predator ~ delta * prey * predator - gamma * predator
    ),
    state_names = c("prey", "predator"),
    parms = list(alpha = 1.0, beta = 0.1, delta = 0.075, gamma = 1.5),
    init = c(prey = 40, predator = 9)
)
result <- dyn_sim(model, t_max = 100, discard_transient = 50)

# Compute spectral analysis
spec <- spectral_analysis(result)
print(spec)
summary(spec)

# Visualize power spectral density
plot(spec)
plot(spec, type = "combined")
} # }
```
