# Complete chameleon attractor analysis

Perform full multiscale analysis of a time series to characterize
potential chameleon attractors. This is the main entry point combining
Takens embedding, MEMD decomposition, and EVT-based metrics computation.

## Usage

``` r
chameleon_analysis(
  x,
  m = NULL,
  tau = NULL,
  n_directions = 64,
  q = 0.98,
  lowpass_freq = NULL,
  dt = 1,
  n_cores = 1,
  statistical_test = FALSE,
  n_surrogates = 199,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- x:

  Numeric vector. Univariate time series to analyze. Must have at least
  200 observations for reliable results. Longer series (n \> 2000)
  provide more robust EVT estimates.

- m:

  Integer. Embedding dimension for Takens reconstruction. If NULL
  (default), estimated automatically using the false nearest neighbors
  method. Typical values range from 2 to 10 depending on the system's
  complexity.

- tau:

  Integer. Time delay in samples for Takens embedding. If NULL
  (default), estimated from the autocorrelation function. The optimal
  tau corresponds to decorrelation of successive coordinates.

- n_directions:

  Integer. Number of projection directions for MEMD envelope computation
  (default 64). More directions improve accuracy but increase
  computation time approximately linearly.

- q:

  Numeric. Quantile for EVT threshold selection (default 0.98). Higher
  values focus on more extreme recurrences. Must be between 0.9 and
  0.999.

- lowpass_freq:

  Numeric. Optional low-pass filter cutoff frequency (Hz). Applied
  before analysis to remove high-frequency noise, following the
  preprocessing in Alberti et al. (2023). Set to NULL to skip filtering.

- dt:

  Numeric. Time step of the series in seconds (default 1). Used for
  frequency estimation in MEMD and scale-dependent metrics.

- n_cores:

  Integer. Number of CPU cores for parallel computation (default 1).
  Currently affects EVT metrics computation.

- statistical_test:

  Logical. Perform rigorous statistical testing using surrogate-based
  permutation tests (default FALSE). When TRUE, the detection is based
  on p-values from
  [`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
  rather than arbitrary thresholds. Recommended for publication-quality
  results.

- n_surrogates:

  Integer. Number of surrogates for statistical testing (default 199).
  Only used when `statistical_test = TRUE`.

- alpha:

  Numeric. Significance level for statistical testing (default 0.05).
  Only used when `statistical_test = TRUE`.

- verbose:

  Logical. Print progress information during analysis (default TRUE).
  Useful for monitoring long-running analyses.

## Value

A list of class `"chameleon_analysis"` containing:

- x_original:

  Numeric vector. Original input time series.

- x_filtered:

  Numeric vector or NULL. Low-pass filtered series (only present if
  `lowpass_freq` was specified).

- embedded:

  Matrix of class `"takens_embedding"`. Reconstructed phase-space
  trajectory with dimensions (n_points, m).

- memd:

  Object of class `"memd"`. MEMD decomposition containing multivariate
  intrinsic mode functions.

- evt_metrics:

  Object of class `"evt_metrics"`. Instantaneous dimension d(t) and
  inverse persistence theta(t) for the full attractor.

- scale_metrics:

  Object of class `"scale_metrics"`. Scale-dependent metrics D(t,f) and
  theta(t,f) across frequency bands.

- params:

  List of analysis parameters (m, tau, n_directions, q, lowpass_freq,
  dt).

- is_chameleon:

  List with detection results including `detected` (logical), `D_range`,
  `theta_range`, and coefficient of variation metrics.

- statistical_test:

  Object of class `"chameleon_test"` (only present if
  `statistical_test = TRUE`). Contains p-values, effect sizes, and
  rigorous statistical conclusions.

## Details

A chameleon attractor is a strange attractor whose geometric and
topological properties (dimension, persistence) vary dynamically across
time and frequency scales. This function implements the complete
methodology from Alberti et al. (2023) to detect and characterize such
behavior.

The analysis follows the seven-step workflow from Alberti et al. (2023):

1.  **Preprocessing**: Optional low-pass filtering to isolate
    dynamically relevant frequencies and reduce noise contamination.

2.  **Parameter estimation**: Automatic estimation of embedding
    dimension (m) and time delay (tau) if not provided.

3.  **Takens embedding**: Reconstruction of the attractor in
    m-dimensional phase space using delay coordinates.

4.  **MEMD decomposition**: Extraction of multivariate intrinsic mode
    functions representing dynamics at different timescales.

5.  **Full-attractor EVT**: Computation of instantaneous dimension and
    persistence for the complete embedded trajectory.

6.  **Scale-dependent metrics**: EVT analysis on partial sums of MIMFs
    to reveal frequency-dependent attractor properties.

7.  **Chameleon detection**: Statistical assessment of whether metrics
    vary significantly across scales.

**Interpretation of results:**

- **Scale-invariant attractor**: D(f) and theta(f) remain relatively
  constant across frequency scales. The attractor has uniform geometric
  properties regardless of the observation timescale.

- **Chameleon attractor**: D(f) and/or theta(f) vary significantly with
  frequency. The attractor appears qualitatively different when observed
  at different timescales, indicating complex multiscale dynamics.

**Computational considerations:**

- EVT metrics computation is O(n^2) in trajectory length

- MEMD is O(n \* log(n) \* n_directions)

- Total runtime scales primarily with n^2

## References

Alberti T, Daviaud F, Donner RV, Dubrulle B, Faranda D, Lucarini V
(2023). Chameleon attractors in turbulent flows. Chaos, Solitons and
Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

Faranda D, Messori G, Yiou P (2017). Dynamical proxies of North Atlantic
predictability and extremes. Scientific Reports 7:41278.
[doi:10.1038/srep41278](https://doi.org/10.1038/srep41278)

## See also

[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
for standalone statistical testing,
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for computing scale metrics.

## Examples

``` r
# Generate test data with multiple scales
set.seed(42)
t <- seq(0, 100, by = 0.01)
# Low frequency component
x_low <- sin(2*pi*0.05*t)
# High frequency component
x_high <- 0.3*sin(2*pi*2*t)
# Noise
x_noise <- 0.1*rnorm(length(t))
x <- x_low + x_high + x_noise

# Perform analysis
if (FALSE) { # \dontrun{
result <- chameleon_analysis(x, verbose = TRUE)
print(result)

# With rigorous statistical testing
result <- chameleon_analysis(x, statistical_test = TRUE, n_surrogates = 199)
print(result$statistical_test)
plot(result$statistical_test)
} # }
```
