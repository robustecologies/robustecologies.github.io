# Print method for a summary.noise_spec object

Formats the sample-diagnostic summary of a `noise_spec` object as a
compact report: distributional parameters, sample mean and variance,
lag-one autocorrelation and empirical spectral slope.

## Usage

``` r
# S3 method for class 'summary.noise_spec'
print(x, ...)
```

## Arguments

- x:

  A `summary.noise_spec` object.

- ...:

  Unused.

## Value

The input `x`, invisibly.

## See also

[`summary.noise_spec()`](https://robustecologies.github.io/janos/reference/summary.noise_spec.md)
. summary constructor.

## Examples

``` r
if (FALSE) { # \dontrun{
print(summary(colored_noise(beta = 1)))
} # }
```
