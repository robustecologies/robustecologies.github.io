# Proposal covariance

This function returns the most recent covariance matrix or a list of
blocking covariance matrices from an object of class `demonoid`, the
most recent covariance matrix from `iterquad`, `laplace`, or `vb`, and a
number of covariance matrices of an object of class `pmc` equal to the
number of mixture components. The returned covariance matrix or matrices
are intended to be the initial proposal covariance matrix or matrices
for future updates. A variance vector from an object of class `demonoid`
is converted to a covariance matrix.

## Usage

``` r
as.covar(x)
```

## Arguments

- x:

  This is an object of class `demonoid`, `iterquad`, `laplace`, `pmc`,
  or `vb`.

## Value

The returned value is a matrix (or array in the case of PMC with
multiple mixture components) of the latest observed or proposal
covariance, which may now be used as an initial proposal covariance
matrix or matrices for a future update.

## Details

Unless it is known beforehand how many iterations are required for
iterative quadrature, Laplace Approximation, or Variational Bayes to
converge, MCMC to appear converged, or the normalized perplexity to
stabilize in PMC, multiple updates are necessary. An additional update,
however, should not begin with the same proposal covariance matrix or
matrices as the original update, because it will have to repeat the work
already accomplished. For this reason, the `as.covar` function may be
used at the end of an update to change the previous initial values to
the latest values.

The `as.covar` function is most helpful with objects of class `pmc` that
have multiple mixture components. For more information, see
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md).

## See also

[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model and data
N <- 25
y <- rnorm(N, 2, 0.5)
Data <- list(N = N, y = y, mon.names = "LP",
  parm.names = c("mu", "log.sigma"))
Model <- function(parm, Data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
  LP <- LL + dnormv(mu, 0, 1000, log = TRUE) +
    dhalfcauchy(sigma, 25, log = TRUE)
  yhat <- rnorm(Data$N, mu, sigma)
  Monitor <- LP
  return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
    yhat = yhat, parm = parm))
}

# Fit model
fit <- lucifer(Model, Data, Initial.Values = c(0, 0),
  Iterations = 1000, Status = 200, Thinning = 1,
  Algorithm = "CHARM", Specs = NULL)

# Extract covariance for next update
Covar <- as.covar(fit)
} # }
```
