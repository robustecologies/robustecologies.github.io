# Initial values

This function returns the most recent posterior samples from an object
of class `demonoid`, the posterior means of an object of class
`iterquad`, the posterior modes of an object of class `laplace` or `vb`,
the posterior means of an object of class `pmc` with one mixture
component, or the latest means of the importance sampling distribution
of an object of class `pmc` with multiple mixture components. The
returned values are intended to be the initial values for future
updates.

## Usage

``` r
as.initial.values(x)
```

## Arguments

- x:

  This is an object of class `demonoid`, `iterquad`, `laplace`, `pmc`,
  or `vb`.

## Value

The returned value is a vector (or matrix in the case of a multi-chain
`demonoid` or `pmc` with multiple mixture components) of the latest
values, which may now be used as initial values for a future update.

## Details

Unless it is known beforehand how many iterations are required for
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
or
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
to converge, MCMC in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
to appear converged, or the normalized perplexity to stabilize in
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md),
multiple updates are necessary. An additional update, however, should
not begin with the same initial values as the original update, because
it will have to repeat the work already accomplished. For this reason,
the `as.initial.values` function may be used at the end of an update to
change the previous initial values to the latest values.

## See also

[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md),
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

# Extract initial values for next update
Initial.Values <- as.initial.values(fit)
} # }
```
