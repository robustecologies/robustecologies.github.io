# Print method for a noise_spec object

Formats a `noise_spec` object as a compact one-screen header naming the
noise family (`"correlated"`, `"levy"`, `"fbm"`, `"colored"`), its
principal parameters (covariance for correlated, stability index and
skewness for Levy, Hurst exponent for fBm, spectral exponent and
amplitude for colored) and any short interpretation of the parameters
(for example the Cauchy case \\\alpha = 1\\ for Levy or the pink-noise
case \\\beta = 1\\ for colored).

## Usage

``` r
# S3 method for class 'noise_spec'
print(x, ...)
```

## Arguments

- x:

  A `noise_spec` object.

- ...:

  Unused, kept for S3 compatibility.

## Value

The input `x`, invisibly.

## See also

[`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
[`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md),
[`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
[`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md)
— constructors;
[`summary.noise_spec()`](https://robustecologies.github.io/janos/reference/summary.noise_spec.md)
— sample statistics;
[`plot.noise_spec()`](https://robustecologies.github.io/janos/reference/plot.noise_spec.md)
— realisation and spectrum panel.

## Examples

``` r
if (FALSE) { # \dontrun{
print(colored_noise(beta = 1, amplitude = 0.1))
} # }
```
