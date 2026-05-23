# Plot samples from the output of MISS

This may be used to plot, or save plots of, samples in an object of
class `miss`. Plots include a trace plot, density plot, and
autocorrelation or ACF plot.

## Usage

``` r
# S3 method for class 'miss'
plot(x, PDF = FALSE, col = NULL, ...)
```

## Arguments

- x:

  This required argument is an object of class `miss`.

- PDF:

  This logical argument indicates whether or not the user wants Lucifer
  to save the plots as a .pdf file.

- col:

  Optional character vector of colors. When `NULL` (default), the RElab
  palette is used.

- ...:

  Additional arguments are unused.

## Value

Invisibly returns `x`. Called for its graphical side effect.

## Details

The plots are arranged in a \\3 \times 3\\ matrix. Each row represents
the predictive distribution of a missing value. The left column displays
trace plots, the middle column displays kernel density plots, and the
right column displays autocorrelation (ACF) plots.

Trace plots show the thinned history of the predictive distribution,
with its value in the y-axis moving by iteration across the x-axis.
Simulations of a predictive distribution with good properties do not
suggest a trend upward or downward as it progresses across the x-axis
(it should appear stationary), and it should mix well. A red, smoothed
line also appears to aid visual inspection.

Kernel density plots depict the marginal posterior distribution. There
is no distributional assumption about this density.

Autocorrelation plots show the autocorrelation or serial correlation
between sampled values at nearby iterations. The ideal autocorrelation
plot shows perfect correlation at zero lag, and quickly falls to zero
autocorrelation for all other lags.

## See also

[`MISS`](https://robustecologies.github.io/lucifer/reference/MISS.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### See the MISS function for an example.
} # }
```
