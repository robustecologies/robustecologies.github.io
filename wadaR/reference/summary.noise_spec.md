# Summary method for a noise_spec object

Reports the distributional and second-order structure of a `noise_spec`
object. For every family, a synthetic realisation of length 4096 is
generated with a fixed seed and its sample mean, variance,
autocorrelation at lag one and empirical spectral slope (when
meaningful) are reported. For `correlated` noise the cross-correlation
matrix is summarised by its spectral abscissa and condition number; for
`levy` noise the tail exponent is the theoretical `alpha`; for `fbm`
noise the Hurst exponent is re-estimated from the sample via the
variance-of-increments method.

## Usage

``` r
# S3 method for class 'noise_spec'
summary(object, ...)
```

## Arguments

- object:

  A `noise_spec` object.

- ...:

  Unused.

## Value

A list of class `summary.noise_spec`, invisibly.

## Details

The sample statistics are reproducible because the realisation is drawn
under `set.seed(1)`. They are diagnostics, not definitions: for Levy
distributions the variance is infinite in theory, so the reported sample
variance is informative only up to the truncation introduced by a finite
realisation.

## See also

[`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
[`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md),
[`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
[`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md)
— constructors;
[`print.noise_spec()`](https://robustecologies.github.io/janos/reference/print.noise_spec.md)
— compact header;
[`plot.noise_spec()`](https://robustecologies.github.io/janos/reference/plot.noise_spec.md)
— realisation and spectrum panel.

## Examples

``` r
if (FALSE) { # \dontrun{
summary(fbm_noise(H = 0.7))
} # }
```
