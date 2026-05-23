# Plot a prescription object

Produces publication-quality visualizations of the algorithm
recommendation. Four plot types are available.

## Usage

``` r
# S3 method for class 'prescription'
plot(x, type = c("default", "scores", "profile", "landscape"), n.top = 15, ...)
```

## Arguments

- x:

  An object of class `prescription`.

- type:

  Character. Plot type: `"default"` (two-panel scores + radar),
  `"scores"` (same as default), `"profile"` (radar fingerprint only), or
  `"landscape"` (bubble chart of algorithms colored by MCMC subcategory
  with evaluation cost as size).

- n.top:

  Integer. Number of top methods to display. Default 15.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns the ggplot object(s).

## Details

Produces diagnostic plots of a Prescribe recommendation object produced
by
[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md).
Summary of the content is given below. Default output renders a
multi-panel graphic (trace, density, and autocorrelation where
applicable). The `PDF` argument captures the graphic to a file;
otherwise the current device is used. Font and colour choices follow
[`theme_relab`](https://robustecologies.github.io/lucifer/reference/theme_relab.md).

## References

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Buerkner, P.-C., & Modrak, M. (2020).
Bayesian workflow. arXiv:2011.01808.

## See also

[`Prescribe`](https://robustecologies.github.io/lucifer/reference/Prescribe.md),
[`print.prescription`](https://robustecologies.github.io/lucifer/reference/print.prescription.md),
[`summary.prescription`](https://robustecologies.github.io/lucifer/reference/summary.prescription.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rx <- Prescribe(Model, Data, Initial.Values)
plot(rx)
plot(rx, type = "landscape", n.top = 30)
plot(rx, type = "profile")
} # }
```
