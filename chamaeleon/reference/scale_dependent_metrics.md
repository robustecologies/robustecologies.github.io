# Compute scale-dependent EVT metrics

Compute instantaneous dimension D(t,f) and inverse persistence
theta(t,f) for scale-filtered signals. This is the core analysis for
detecting chameleon behavior, following the methodology of Alberti et
al. (2023).

## Usage

``` r
scale_dependent_metrics(
  x,
  memd_obj = NULL,
  freq_cutoffs = NULL,
  q = 0.98,
  n_directions = 64,
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- x:

  Numeric matrix. Embedded phase-space trajectory, typically the output
  of
  [`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md).
  Must have at least 100 rows and 2 columns.

- memd_obj:

  Optional object of class `"memd"`. Pre-computed MEMD decomposition
  from
  [`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md).
  If NULL (default), MEMD is computed internally with the specified
  parameters.

- freq_cutoffs:

  Numeric vector. Frequency cutoffs (Hz) at which to compute
  scale-dependent metrics. If NULL (default), uses the frequencies of
  the extracted MIMFs. Should be in descending order for proper
  interpretation.

- q:

  Numeric. Quantile for EVT threshold selection (default 0.98). Same
  interpretation as in
  [`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md).

- n_directions:

  Integer. Number of MEMD projection directions (default 64). Only used
  if `memd_obj` is NULL.

- n_cores:

  Integer. Number of CPU cores for parallel computation (default 1).

- verbose:

  Logical. Print progress information (default FALSE).

## Value

A list of class `"scale_metrics"` containing:

- D:

  Matrix of dimensions (n_times, n_freqs). Scale-dependent instantaneous
  dimension D(t,f). Each column corresponds to a frequency cutoff, each
  row to a time point.

- theta:

  Matrix of dimensions (n_times, n_freqs). Scale-dependent inverse
  persistence theta(t,f).

- frequencies:

  Numeric vector. Frequency cutoffs used (Hz).

- timescales:

  Numeric vector. Corresponding timescales (1/frequency).

- mean_D:

  Numeric vector. Time-averaged dimension at each frequency, \<D(f)\>.
  This is the key diagnostic for chameleon detection.

- mean_theta:

  Numeric vector. Time-averaged persistence at each frequency,
  \<theta(f)\>.

- sd_D:

  Numeric vector. Standard deviation of D(t,f) at each frequency.

- sd_theta:

  Numeric vector. Standard deviation of theta(t,f) at each frequency.

- memd:

  Object of class `"memd"`. The MEMD decomposition used.

- q:

  Numeric. Quantile used for EVT threshold.

S3 methods available: [`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

The key insight is that a chameleon attractor exhibits different
geometric and dynamical properties at different frequency scales. By
filtering the embedded signal to retain only modes above a frequency
cutoff, we can characterize how the attractor properties change with
observation scale.

**Methodology:**

For each frequency cutoff f\*, the algorithm:

1.  Reconstructs a scale-filtered signal using partial sums of MIMFs:
    \$\$\Theta^f\_\mu(t) = \sum\_{k: f_k \> f^\*} C\_{\mu,k}(t)\$\$
    where \\C\_{\mu,k}(t)\\ are the multivariate intrinsic mode
    functions with frequencies f_k.

2.  Computes EVT metrics (d, theta) on this filtered signal.

3.  Stores the results as D(t,f\*) and theta(t,f\*).

**Interpretation:**

- If \<D(f)\> is constant across frequencies, the attractor geometry is
  scale-invariant.

- If \<D(f)\> varies significantly (increasing or decreasing with f),
  the attractor exhibits chameleon behavior, appearing differently at
  different timescales.

- Similarly for \<theta(f)\>: constant implies uniform dynamics,
  variable implies scale-dependent persistence.

**Computational notes:** This function performs EVT analysis n_freqs
times, so runtime scales as O(n_freqs \* n^2) where n is the trajectory
length. Progress reporting is enabled with `verbose=TRUE`.

## References

Alberti T, Daviaud F, Donner RV, Dubrulle B, Faranda D, Lucarini V
(2023). Chameleon attractors in turbulent flows. Chaos, Solitons and
Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

## See also

[`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
for the MEMD decomposition,
[`evt_metrics`](https://robustecologies.github.io/chamaeleon/reference/evt_metrics.md)
for the underlying EVT computation,
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for the complete workflow.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
t <- seq(0, 50, by = 0.01)
x <- sin(2*pi*0.1*t) + 0.5*sin(2*pi*t) + 0.2*rnorm(length(t))
embedded <- takens_embed(x, m = 3, tau = 15)

# Compute scale-dependent metrics
scale_met <- scale_dependent_metrics(embedded, verbose = TRUE)
} # }
```
