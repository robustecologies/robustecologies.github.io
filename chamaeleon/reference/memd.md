# Multivariate empirical mode decomposition

Decompose a multivariate signal into multivariate intrinsic mode
functions (MIMFs) using the MEMD algorithm. This data-driven technique
separates oscillatory components at different characteristic timescales
while preserving mode alignment across all channels.

## Usage

``` r
memd(
  x,
  n_directions = 64,
  max_imf = NULL,
  max_sift = 100,
  sift_thresh = 0.05,
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- x:

  Numeric matrix. Multivariate signal with dimensions (n_samples,
  n_variables). In the context of chameleon analysis, this is typically
  the output from
  [`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md),
  where each column represents a delay coordinate.

- n_directions:

  Integer. Number of projection directions for envelope estimation
  (default 64). The MEMD algorithm projects the multivariate signal onto
  random directions on the unit hypersphere. More directions improve the
  accuracy of envelope estimation but increase computation time.
  Recommended range: 32-128.

- max_imf:

  Integer. Maximum number of IMFs to extract. Default is NULL, which
  extracts modes until the residual becomes monotonic. Typically yields
  log2(n) + 1 modes.

- max_sift:

  Integer. Maximum sifting iterations per IMF extraction (default 100).
  Sifting refines each mode until it satisfies IMF criteria.

- sift_thresh:

  Numeric. Stopping criterion threshold for the sifting process (default
  0.05). Based on the normalized mean squared difference between
  successive sifting iterations. Lower values yield more refined modes
  but require more iterations.

- n_cores:

  Integer. Number of CPU cores for parallel processing (default 1).
  Parallelization is applied to envelope computation.

- verbose:

  Logical. Print progress information (default FALSE).

## Value

A list of class `"memd"` containing:

- mimfs:

  List of N matrices. Each element is a MIMF with the same dimensions as
  the input (n_samples x n_variables). Ordered from highest to lowest
  frequency.

- residual:

  Matrix. The monotonic residual after all IMFs are extracted,
  representing the underlying trend.

- frequencies:

  Numeric vector of length N. Estimated mean frequency (Hz) of each
  MIMF, computed from zero-crossing counts.

- timescales:

  Numeric vector of length N. Corresponding mean timescales
  (1/frequency) of each MIMF.

- n_imfs:

  Integer. Number of extracted intrinsic mode functions.

S3 methods available: [`print()`](https://rdrr.io/r/base/print.html).

## Details

The signal is decomposed as: \$\$\Theta\_\mu(t) = \sum\_{k=1}^{N}
C\_{\mu,k}(t) + R\_\mu(t)\$\$

where each \\C\_{\mu,k}(t)\\ is a multivariate intrinsic mode function
representing oscillations at a characteristic frequency, and
\\R\_\mu(t)\\ is the monotonic residual (trend).

**Algorithm overview:**

The MEMD algorithm extends univariate EMD to handle multivariate signals
while ensuring mode alignment, i.e., corresponding modes in each channel
have the same center frequency. The key steps are:

1.  Generate uniformly distributed projection directions on the
    p-dimensional unit hypersphere.

2.  For each direction, project the multivariate signal to obtain a
    scalar projection.

3.  Find local maxima and minima of each projection.

4.  Interpolate through extrema using cubic splines to form upper and
    lower envelopes.

5.  Compute the mean envelope as the average across all directions.

6.  Subtract the mean envelope from the signal.

7.  Repeat (sifting) until the result satisfies IMF criteria: number of
    extrema and zero crossings differ by at most one.

8.  Subtract the extracted IMF from the original and repeat for the next
    mode.

**IMF properties:** Each extracted MIMF is nearly monocomponent, meaning
it contains oscillations in a narrow frequency band. The modes are
ordered from highest to lowest frequency after extraction.

**Connection to scale-dependent metrics:** In chameleon analysis, the
MIMFs enable construction of scale-filtered signals by summing modes
above a frequency cutoff. EVT metrics computed on these filtered signals
reveal how attractor properties change with observation scale.

**Computational complexity:** O(n \* log(n) \* n_directions \* n_imfs),
where n is signal length. The dominant cost is cubic spline
interpolation for each direction.

## References

Rehman N, Mandic DP (2010). Multivariate empirical mode decomposition.
Proceedings of the Royal Society A 466:1291-1302.
[doi:10.1098/rspa.2009.0502](https://doi.org/10.1098/rspa.2009.0502)

Huang NE et al. (1998). The empirical mode decomposition and the Hilbert
spectrum for nonlinear and non-stationary time series analysis.
Proceedings of the Royal Society A 454:903-995.
[doi:10.1098/rspa.1998.0193](https://doi.org/10.1098/rspa.1998.0193)

Alberti T et al. (2023). Chameleon attractors in turbulent flows. Chaos,
Solitons and Fractals 168:113195.
[doi:10.1016/j.chaos.2023.113195](https://doi.org/10.1016/j.chaos.2023.113195)

## See also

[`memd_partial_sums`](https://robustecologies.github.io/chamaeleon/reference/memd_partial_sums.md)
for constructing scale-filtered signals,
[`scale_dependent_metrics`](https://robustecologies.github.io/chamaeleon/reference/scale_dependent_metrics.md)
for the complete analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create embedded signal
set.seed(42)
t <- seq(0, 20, by = 0.01)
x <- sin(2*pi*0.5*t) + 0.5*sin(2*pi*2*t) + 0.2*rnorm(length(t))
embedded <- takens_embed(x, m = 3, tau = 10)

# Decompose using MEMD
decomp <- memd(embedded, verbose = TRUE)

# Inspect results
print(decomp)
summary(decomp)

# Visualize decomposition
plot(decomp, type = "modes")
plot(decomp, type = "spectrum")
plot(decomp, type = "reconstruction")
} # }
```
