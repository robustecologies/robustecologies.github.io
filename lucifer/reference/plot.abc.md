# Plot method for ABC objects

Produces diagnostic plots for an ABC fit.

## Usage

``` r
# S3 method for class 'abc'
plot(
  x,
  type = "posterior",
  ground_truth = NULL,
  true_values = NULL,
  col = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `abc`.

- type:

  Character: `"posterior"` (default), `"distances"`, `"temperature"`
  (SA-ABC only), `"intervals"`, or `"pairs"`.

- ground_truth:

  Optional named numeric vector of true parameter values. When provided
  and `type = "posterior"`, vertical dashed lines are drawn at the true
  values on each posterior histogram.

- true_values:

  Alias for `ground_truth` (backward compatibility).

- col:

  Optional character vector of hex color strings. When non-`NULL`,
  overrides the default RElab contrasting palette.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns the plot(s).

## Details

Produces diagnostic plots of an approximate Bayesian computation fit
produced by
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate
Bayesian computation in population genetics. *Genetics*, 162(4),
2025-2035.
[doi:10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025)

## See also

[`print.abc`](https://robustecologies.github.io/lucifer/reference/print.abc.md),
[`summary.abc`](https://robustecologies.github.io/lucifer/reference/summary.abc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving plot.abc
} # }
```
