# Plot method for a noise_spec object

Two-panel diagnostic of a `noise_spec` object. The top panel shows a
sample realisation of length 4096 drawn under a fixed seed; for
multivariate `correlated` noise all components are overlaid. The bottom
panel shows either the log-log periodogram of the sample (for scalar
noise families) or the empirical cross-correlation matrix (for
multivariate correlated noise), providing a visual check of the spectral
or covariance structure.

## Usage

``` r
# S3 method for class 'noise_spec'
plot(x, type = c("default", "acf"), title = NULL, ...)
```

## Arguments

- x:

  A `noise_spec` object.

- type:

  Character, either `"default"` (realisation + spectrum) or `"acf"`
  (realisation + empirical autocorrelation).

- title:

  Optional plot title (overrides the default).

- ...:

  Unused, kept for S3 compatibility.

## Value

A ggplot object (patchwork of two panels for scalar families; single
ggplot with raster for multivariate correlated noise).

## Details

The realisation is generated once with `set.seed(1)`; every call
produces the same panel. The periodogram is computed with no tapering.
For Levy noise, amplitude-dependent outliers are expected and the
periodogram is not informative in the usual sense; the panel
nevertheless keeps the same structure so the comparison across noise
families is uniform.

## See also

[`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
[`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md),
[`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
[`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md)
. constructors;
[`print.noise_spec()`](https://robustecologies.github.io/janos/reference/print.noise_spec.md)
. compact header;
[`summary.noise_spec()`](https://robustecologies.github.io/janos/reference/summary.noise_spec.md)
. sample statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(colored_noise(beta = 1))
plot(fbm_noise(H = 0.7))
} # }
```
