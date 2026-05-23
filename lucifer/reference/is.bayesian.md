# Logical check of a Bayesian model

This function provides a logical test of whether or not a `Model`
specification function is Bayesian.

## Usage

``` r
is.bayesian(Model, Initial.Values, Data)
```

## Arguments

- Model:

  A model specification function. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

- Initial.Values:

  A vector of initial values, or current parameter values. For more
  information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

- Data:

  A list of data. For more information, see the
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
  function.

## Value

A logical value of `TRUE` when the model is Bayesian, and `FALSE`
otherwise.

## Details

This function tests whether or not a model is Bayesian by comparing the
first two returned arguments: the logarithm of the unnormalized joint
posterior density (`LP`) and deviance (`Dev`). The deviance (D) is

\$\$\mathrm{D} = -2 \mathrm{LL}\$\$

where LL is the log-likelihood. Consequently,

\$\$\mathrm{LL} = \mathrm{D} / -2\$\$

and LP is the sum of LL and prior probability densities. If LP = LL,
then the model is not Bayesian, because prior densities are absent.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define a Bayesian model (includes prior densities)
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

# Check if the model is Bayesian
is.bayesian(Model, Initial.Values = c(0, 0), Data) # TRUE
} # }
```
