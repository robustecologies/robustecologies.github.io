# Plot MEMD decomposition

Visualize the multivariate empirical mode decomposition results using
ggplot2. Three visualization types are available: intrinsic mode
functions, power spectra, and cumulative reconstruction.

## Usage

``` r
# S3 method for class 'memd'
plot(
  x,
  type = c("modes", "spectrum", "reconstruction", "all"),
  channel = 1,
  dt = 1,
  max_modes = NULL,
  ...
)
```

## Arguments

- x:

  Object of class `memd` returned by
  [`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md).

- type:

  Character. Type of visualization:

  "modes"

  :   (Default) Faceted time series of all MIMFs and residual, labeled
      with their characteristic frequencies.

  "spectrum"

  :   Power spectral density of each MIMF, showing the frequency content
      of each mode.

  "reconstruction"

  :   Cumulative partial sums showing how progressively adding
      lower-frequency modes reconstructs the original signal.

  "all"

  :   Combined panel showing modes (left), spectrum (top-right), and
      reconstruction (bottom-right) in a single figure.

- channel:

  Integer. Which channel (variable) of the multivariate signal to
  display. Default is 1 (first channel).

- dt:

  Numeric. Sampling interval in seconds for proper frequency axis
  labeling. Default is 1 (samples).

- max_modes:

  Integer. Maximum number of modes to display. Default is NULL (show
  all). Useful when many modes are extracted.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object that can be further customized or printed.

## Details

The MEMD decomposes a multivariate signal into intrinsic mode functions
(MIMFs) representing oscillations at different characteristic
timescales. These visualizations help assess the quality of the
decomposition and understand the multiscale structure of the signal.

**Mode visualization:** Each MIMF is shown as a separate facet, ordered
from highest to lowest frequency. The facet labels indicate the
estimated mean frequency. The residual (monotonic trend) is shown in the
bottom panel.

**Spectral visualization:** Power spectral density computed via
periodogram for each MIMF. Well-separated modes should show distinct
spectral peaks with minimal overlap.

**Reconstruction visualization:** Shows how the signal is progressively
reconstructed by cumulatively adding modes from highest to lowest
frequency. The final panel should closely match the original signal.

## References

Rehman N, Mandic DP (2010). Multivariate empirical mode decomposition.
Proceedings of the Royal Society A 466:1291-1302.
[doi:10.1098/rspa.2009.0502](https://doi.org/10.1098/rspa.2009.0502)

## See also

[`memd`](https://robustecologies.github.io/chamaeleon/reference/memd.md)
for computing the decomposition,
[`memd_partial_sums`](https://robustecologies.github.io/chamaeleon/reference/memd_partial_sums.md)
for scale-filtered signals.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test signal with multiple scales
set.seed(42)
t <- seq(0, 20, by = 0.01)
x <- sin(2*pi*0.5*t) + 0.3*sin(2*pi*2*t) + 0.1*rnorm(length(t))
embedded <- takens_embed(x, m = 3, tau = 10)

# Compute MEMD
decomp <- memd(embedded)

# Visualize modes
plot(decomp, type = "modes")

# Visualize spectra
plot(decomp, type = "spectrum", dt = 0.01)

# Visualize reconstruction
plot(decomp, type = "reconstruction")
} # }
```
