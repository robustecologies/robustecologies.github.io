# Minnesota prior

The Minnesota prior, also called the Litterman prior, is a shrinkage
prior for autoregressive parameters in vector autoregressive (VAR)
models. There are many variations of the Minnesota prior. This Minnesota
prior is calculated as presented in Lutkepohl (2005, p. 225), and
returns one or more prior covariance matrices in an array.

## Usage

``` r
MinnesotaPrior(J, lags = c(1, 2), lambda = 1, theta = 0.5, sigma)
```

## Arguments

- J:

  The scalar number of time-series in the VAR.

- lags:

  An integer vector of lags of the autoregressive parameters. The lags
  are not required to be successive.

- lambda:

  A scalar, positive-only hyperparameter that controls how tightly the
  parameter of the first lag is concentrated around zero. A smaller
  value results in smaller diagonal variance. When equal to zero, the
  posterior equals the prior and data is not influential. When equal to
  infinity, no shrinkage occurs and posterior expectations are closest
  to estimates from ordinary least squares (OLS). It has been asserted
  that as the number, \\J\\, of time-series increases, this
  hyperparameter should decrease.

- theta:

  A scalar hyperparameter in the interval \[0,1\]. When one,
  off-diagonal elements have variance similar or equal to diagonal
  elements. When zero, off-diagonal elements have zero variance. A
  smaller value is associated with less off-diagonal variance.

- sigma:

  A vector of length \\J\\ of residual standard deviations of the
  dependent variables given the expectations.

## Value

A \\J \times J \times L\\ array for \\J\\ time-series and \\L\\ lags.

## Details

The Minnesota prior was introduced in Doan, Litterman, and Sims (1984)
as a shrinkage prior for autoregressive parameters in vector
autoregressive (VAR) models. The Minnesota prior was reviewed in
Litterman (1986), and numerous variations have been presented since.
This is the version of the Minnesota prior as described in Lutkepohl
(2005, p. 225) for stationary time-series.

Given one or more \\J \times J\\ matrices of autoregressive parameters
in a VAR model, the user specifies two tuning hyperparameters for the
Minnesota prior: `lambda` and `theta`. Each iteration of the numerical
approximation algorithm, the latest vector of residual standard
deviation parameters is supplied to the `MinnesotaPrior` function, which
then returns an array that contains one or more prior covariance
matrices for the autoregressive parameters. Multiple prior covariance
matrices are returned when multiple lags are specified. The tuning
hyperparameters, `lambda` and `theta`, can be estimated from the data
via hierarchical Bayes.

It is important to note that the Minnesota prior does not technically
return a covariance matrix, because the matrix is not symmetric, and
therefore not positive-definite. For this reason, a Minnesota prior
covariance matrix should not be supplied as a covariance matrix to a
multivariate normal distribution, such as with the
[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md)
function, though it would be accepted and then (incorrectly) converted
to a symmetric matrix. Instead,
[`dnormv`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Variance.md)
should be used for element-wise evaluation.

While the Minnesota prior is used to specify the prior covariance for
VAR autoregressive parameters, prior means are often all set to zero, or
sometimes the first lag is set to an identity matrix.

An example is provided in the Examples vignette.

## References

Doan, T., Litterman, R.B. and Sims, C.A. (1984). "Forecasting and
Conditional Projection using Realistic Prior Distributions".
*Econometric Reviews*, 3, p. 1–144.

Litterman, R.B. (1986). "Forecasting with Bayesian Vector
Autoregressions - Five Years of Experience". *Journal of Business &
Economic Statistics*, 4, p. 25–38.

Lutkepohl, H. (2005). "New Introduction to Multiple Time Series
Analysis". Springer, Germany.

## See also

[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md),
[`dnormv`](https://robustecologies.github.io/lucifer/reference/dist.Normal.Variance.md),
and
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving MinnesotaPrior
} # }
```
