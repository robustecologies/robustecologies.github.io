# Centering and scaling

This function either centers and scales a continuous variable and
provides options for binary variables, or returns an untransformed
variable from a centered and scaled variable.

## Usage

``` r
CenterScale(x, Binary = "none", Inverse = FALSE, mu, sigma, Range, Min)
```

## Arguments

- x:

  This is a vector to be centered and scaled, or to be untransformed if
  `Inverse=TRUE`.

- Binary:

  This argument indicates how binary variables will be treated, and
  defaults to `"none"`, which keeps the original scale, or transforms
  the variable to the 0-1 range, if not already there. With `"center"`,
  it will center the binary variable by subtracting the mean. With
  `"center0"`, it centers the binary variable at zero, recoding a 0 to
  -0.5, and a 1 to 0.5. Finally, `"centerscale"` will center and scale
  the binary variable, subtracting the mean and dividing by two standard
  deviations.

- Inverse:

  Logical. If `TRUE`, then a centered and scaled variable `x` will be
  transformed to its original, un-centered and un-scaled state. This
  defaults to `FALSE`.

- mu:

  The mean of the original `x`. Required only when `Inverse=TRUE`.

- sigma:

  The standard deviation of the original `x`. Required only when
  `Inverse=TRUE`.

- Range:

  The range of the original `x`. Required only when `Inverse=TRUE` and
  `Binary="none"` or `Binary="center0"`.

- Min:

  The minimum of the original `x`. Required only when `Inverse=TRUE` and
  `Binary="none"` or `Binary="center0"`.

## Value

The `CenterScale` function returns a centered and scaled vector, or the
untransformed vector.

## Details

Gelman (2008) recommends centering and scaling continuous predictors to
facilitate MCMC convergence and enable comparisons between coefficients
of centered and scaled continuous predictors with coefficients of
untransformed binary predictors. A continuous predictor is centered and
scaled as follows: `x.cs <- (x - mean(x)) / (2*sd(x))`. This is an
improvement over the usual practice of standardizing predictors, which
is `x.z <- (x - mean(x)) / sd(x)`, where coefficients cannot be validly
compared between binary and continuous predictors.

In MCMC, such as in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
a centered and scaled predictor often results in a higher effective
sample size
([`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md)),
and therefore the chain mixes better. Centering and scaling is a method
of re-parameterization to improve mixing.

Griffin and Brown (2013) also assert that the user may not want to scale
predictors that are measured on the same scale, since scaling in this
case may increase noisy, low signals. In this case, centering (without
scaling) is recommended. To center a predictor, subtract its mean.

## References

Gelman, A. (2008). "Scaling Regression Inputs by Dividing by Two
Standard Devations". *Statistics in Medicine*, 27, p. 2865–2873.

Griffin, J.E. and Brown, P.J. (2013) "Some Priors for Sparse Regression
Modelling". *Bayesian Analysis*, 8(3), p. 691–702.

## See also

[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
and [`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## Examples

``` r
if (FALSE) { # \dontrun{
### See the lucifer function for an example in use.
library(lucifer)
x <- rnorm(100,10,1)
x.cs <- CenterScale(x)
x.orig <- CenterScale(x.cs, Inverse=TRUE, mu=mean(x), sigma=sd(x))
} # }
```
