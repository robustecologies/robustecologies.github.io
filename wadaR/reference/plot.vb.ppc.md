# Plots of posterior predictive checks for vb

This may be used to plot, or save plots of, samples in an object of
class `vb.ppc`. A variety of plots is provided.

## Usage

``` r
# S3 method for class 'vb.ppc'
plot(x, Style = NULL, Data = NULL, Rows = NULL, PDF = FALSE, ...)
```

## Arguments

- x:

  This required argument is an object of class `vb.ppc`.

- Style:

  This optional argument specifies one of several styles of plots, and
  defaults to `NULL` (which is the same as `"Density"`). Styles of plots
  are indicated in quotes. Optional styles include `"Covariates"`,
  `"Covariates, Categorical DV"`, `"Density"`, `"DW"`,
  `"DW, Multivariate, C"`, `"ECDF"`, `"Fitted"`,
  `"Fitted, Multivariate, C"`, `"Fitted, Multivariate, R"`,
  `"Jarque-Bera"`, `"Jarque-Bera, Multivariate, C"`, `"Mardia"`,
  `"Predictive Quantiles"`, `"Residual Density"`,
  `"Residual Density, Multivariate, C"`,
  `"Residual Density, Multivariate, R"`, `"Residuals"`,
  `"Residuals, Multivariate, C"`, `"Residuals, Multivariate, R"`,
  `"Space-Time by Space"`, `"Space-Time by Time"`, `"Spatial"`,
  `"Spatial Uncertainty"`, `"Time-Series"`,
  `"Time-Series, Multivariate, C"`, and
  `"Time-Series, Multivariate, R"`.

- Data:

  This optional argument accepts the data set used when updating the
  model. Data is required only with certain plot styles.

- Rows:

  This optional argument is for a vector of row numbers that specify the
  records associated by row in the object of class `vb.ppc`. Only these
  rows are plotted. The default is to plot all rows. Some plots do not
  allow rows to be specified.

- PDF:

  This logical argument indicates whether or not the user wants Lucifer
  to save the plots as a .pdf file.

- ...:

  Additional arguments are unused.

## Value

See Details.

## Details

This function can be used to produce a variety of posterior predictive
plots, and the style of plot is selected with the `Style` argument. See
the full documentation of posterior predictive check plot styles in the
manual.

## References

Durbin, J., and Watson, G.S. (1950). "Testing for Serial Correlation in
Least Squares Regression, I." *Biometrika*, 37, p. 409–428.

Jarque, C.M. and Bera, A.K. (1980). "Efficient Tests for Normality,
Homoscedasticity and Serial Independence of Regression Residuals".
*Economics Letters*, 6(3), p. 255–259.

Mardia, K.V. (1970). "Measures of Multivariate Skewness and Kurtosis
with Applications". *Biometrika*, 57(3), p. 519–530.

## See also

[`predict.vb`](https://robustecologies.github.io/lucifer/reference/predict.vb.md)
and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### See the VariationalBayes function for an example.
} # }
```
