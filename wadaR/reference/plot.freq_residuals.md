# Plot freq residuals

S3 method: apply
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to objects of
class `freq_residuals`.

## Usage

``` r
# S3 method for class 'freq_residuals'
plot(x, col = NULL, ...)
```

## Arguments

- x:

  An object of class `freq_residuals`.

- col:

  Optional character. Override the default residual colour.

- ...:

  Currently unused.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

Produces a four-panel residual diagnostic: residuals versus fitted
values, Q-Q plot against theoretical Normal quantiles, histogram of
residuals, and scale-location plot.

## See also

[`freq_residuals`](https://robustecologies.github.io/lucifer/reference/freq_residuals.md),
[`print.freq_residuals`](https://robustecologies.github.io/lucifer/reference/print.freq_residuals.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.freq_residuals
} # }
```
